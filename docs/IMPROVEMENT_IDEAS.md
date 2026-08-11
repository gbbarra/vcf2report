# Ideias de melhoria — validadas em uso real (WGS, Genoma_low_Sample)

Três aprendizados de uma sessão real rodando o vcf2report num VCF de WGS (~5M variantes) numa
máquina de 8GB de RAM. Os três já foram prototipados manualmente (fora do pipeline) e funcionaram —
aqui vai o que cada um precisa pra virar parte do produto.

---

## 1. Suporte nativo a WGS sem estourar RAM

**Problema:** rodar o pipeline padrão (`annotate_vcf.sh` + `run_headless.py`) direto num VCF de
genoma inteiro trava a máquina — confirmado duas vezes, em duas etapas diferentes, numa máquina de
8GB.

**Causa raiz dupla:**
- `scripts/annotate_vcf.sh` roda o SnpEff com `java -Xmx8g` **fixo** (todos os 3 branches de
  resolução do jar/PATH fixam esse valor — não tem flag de override). Numa máquina de 8GB, só o
  heap já é maior que a RAM inteira.
- `run_headless.py` → `pipeline.py` → `vcf/parse.py::parse_vcf()` carrega **todas** as variantes
  numa lista Python de uma vez (não é streaming/generator). ~5M objetos de variante estoura RAM
  independente das stores DuckDB/Parquet (essas são eficientes via scan colunar; o gargalo é a
  lista em memória).

**O que funcionou (prototipado manualmente nesta sessão):**
1. Escopar aos 25 cromossomos primários (chr1-22, X, Y, M) — cobrem ~98,6% das variantes num VCF
   de WGS típico; o resto (alt/decoy/unplaced/HLA, ~195 contigs) não compensa o custo.
2. Anotar por cromossomo: `bcftools view -r <chrom>` → `bcftools norm -m -any` → renomear
   chr↔Ensembl (a DB MANE é Ensembl-style) → SnpEff com heap reduzido (`-Xmx3g` funcionou; pico
   medido ~1,5GB de RSS até no maior cromossomo, chr1) → renomear de volta → concatenar tudo no
   final com `bcftools concat` (ordem cromossômica, não alfabética).
3. Classificar por cromossomo também: `run_headless.py` direto em cada `<chrom>.ann.vcf.gz`
   separadamente (pico medido ~515MB em chr21). Classificação ACMG é local por variante — não
   depende de contexto cross-cromossomo — então é seguro rodar assim.
4. Mesclar os `results.json` no final: somar contadores inteiros do funil de QC, concatenar as
   listas (`classifications`, cada bucket, `qc_rescued`, `clinvar_do_not_dismiss`), e em
   `seq_quality` fazer média ponderada por `n_variants` nas métricas que já são médias/medianas
   (nunca média simples das médias por cromossomo).

**Proposta pro produto:**
- Detectar automaticamente escala do input (contagem de variantes/tamanho do arquivo) e, acima de
  um limiar (ex: >500k variantes, ou heurística de RAM disponível vs. tamanho do VCF), ativar
  automaticamente um modo "chunked" que faz o que está descrito acima — sem o usuário precisar
  escrever scripts manuais.
- No mínimo, expor uma flag de heap configurável em `annotate_vcf.sh` (`SNPEFF_XMX` ou similar) em
  vez do `-Xmx8g` fixo, e documentar a estratégia de chunking em `docs/ANNOTATION.md` /
  `docs/ARCHITECTURE.md` como caminho recomendado pra WGS.
- Considerar tornar `parse_vcf()` um generator/streaming em vez de carregar tudo em lista, o que
  resolveria a raiz do problema do Estágio 6 sem precisar de chunking externo.

---

## 2. Doença associada ao gene no laudo

**Problema:** o laudo mostra o gene e a classificação ACMG, mas não diz **a qual doença aquele
gene está associado** — informação básica que hoje fica implícita (o revisor precisa saber de
cabeça ou pesquisar por fora).

**O que funcionou (prototipado manualmente):** cruzar o gene contra o arquivo oficial
`genes_to_phenotype.txt` do HPO/Monarch (https://github.com/obophenotype/human-phenotype-ontology
/releases, arquivo `genes_to_phenotype.txt` — colunas `gene_symbol`, `disease_id`,
`hpo_id`/`hpo_name`) pra pegar os `disease_id` (OMIM/Orphanet) associados ao gene, e cruzar esses
IDs contra `phenotype.hpoa` (mesmo release, colunas `database_id`↔`disease_name`) pra ter o nome
legível da doença.

**Proposta pro produto:** baixar/congelar esses dois arquivos HPO como uma store adicional (mesmo
padrão das stores gnomAD/AlphaMissense/ClinVar já existentes — Parquet local, sem rede em runtime),
e no card de cada candidato no laudo adicionar uma linha "Doença associada ao gene" com nome +
disease_id. Fácil de já ter pronto porque o pipeline já lê HPO pra fenótipo do paciente — é
reaproveitar a mesma fonte de dados olhando pra doença em vez de fenótipo.

---

## 3. Tipo de herança no laudo — por doença, não por gene

**Problema:** quando adicionamos a doença associada, a tentação é resumir a herança num único rótulo
por **gene**. Isso é enganoso: um mesmo gene frequentemente causa doenças diferentes com herança
diferente. Exemplo real encontrado nesta sessão — **JUP**: Naxos disease é **Autossômica
Recessiva**, mas Displasia Arritmogênica de Ventrículo Direito 12 (mesmo gene, doença diferente) é
**Autossômica Dominante**. Um rótulo agregado "AD/AR" no card esconderia justamente a informação
que mais importa pra interpretação.

**O que funcionou (prototipado manualmente):** dentro de `genes_to_phenotype.txt`, filtrar as
linhas por **gene E disease_id juntos** (não só gene) e checar quais termos HPO de herança aparecem
pra aquele par específico:
- `HP:0000006` → Autossômica Dominante
- `HP:0000007` → Autossômica Recessiva
- `HP:0001417` / `HP:0001419` / `HP:0001423` → Ligada ao X (genérica/recessiva/dominante)
- `HP:0001427` → Mitocondrial
- (lista completa de termos de herança na HPO em si, categoria "Mode of inheritance" HP:0000005)

**Proposta pro produto:** ao implementar o item 2 acima, estruturar o dado como
`{gene: [{disease_id, disease_name, inheritance: [...]}]}` — herança sempre amarrada à doença
específica, nunca um campo solto no nível do gene. No laudo, cada doença listada no card carrega
sua própria herança ao lado do nome.

---

## 4. Padronizar os campos de qualidade em todos os cards do laudo, não só nos sinalizados

**Problema:** no laudo, só os achados sinalizados manualmente (ex: QC-rescued por `MosaicLowAF`)
mostravam os campos de qualidade bruta da chamada — Zigosidade, Fração alélica, GQ, QUAL, Filter.
Os cards de candidatos "normais" (VUS, Pathogenic, etc.) só traziam campos derivados de anotação
(ClinVar, gnomAD AF, AlphaMissense) — sem dar pro revisor visibilidade da confiança/qualidade da
chamada em si, que é justamente o que diferencia um achado sólido (AB~50%, GQ alto, PASS) de um
duvidoso (AB baixo, GQ baixo, filtro suspeito) mesmo dentro do mesmo tier ACMG.

**O que funcionou (prototipado manualmente):** adicionar ao `.kv facts` de **todo** card os mesmos
5 campos que já existiam só nos achados sinalizados — Zigosidade, Fração alélica (formatada como
"XX,X% (n/N reads)", derivada de `allele_balance × depth`), GQ, QUAL, Filter — antes dos campos de
anotação (ClinVar/gnomAD/AlphaMissense). QUAL não vem em `results.json` (só depth/gq/allele_balance
/filter_status); precisou puxar do VCF anotado por posição pra cada variante.

**Proposta pro produto:**
- Incluir `qual` no dict `variant` do `results.json` (ao lado de `depth`/`gq`/`allele_balance`/
  `filter_status`), pra não precisar reconsultar o VCF depois.
- No template do laudo (`report_template.html`), tornar esse bloco de 5 campos de qualidade um
  padrão fixo em **todo** `.vcard`, independente de tier ou bucket — não só nos casos QC-rescued.
  Isso deixa o laudo auto-suficiente pra julgar confiança da chamada sem abrir o VCF original.

---

## 5. QC-rescue está restrito demais — só pega ClinVar Pathogenic/LP, perde VUS de alto valor

**Problema (achado ao vivo, não hipotético):** o usuário desconfiou de um resultado do laudo
("nenhum candidato pro fenótipo de aneurisma de aorta") e pediu pra checar de novo. Reinspeção
manual do VCF anotado (não só da lista de `candidates` já filtrada) achou **MYH11 c.2281T>C
(p.Tyr761His)** — ultra-raro (gnomAD AF=3e-06), **ClinVar VUS com 2 estrelas
especificamente pra "Familial thoracic aortic aneurysm and aortic dissection"** (a condição do
fenótipo informado!), MYH11 sendo gene estabelecido de FTAAD4. Essa variante nunca apareceu nos
956 candidatos porque tem GQ=14 (abaixo do `min_GQ=20` do pipeline) — e como não é ClinVar
**Pathogenic/Likely Pathogenic** (só VUS), não bateu no critério do mecanismo de `qc_rescued` que
salvou MBD4/KCNQ2 (esses sim P/LP) na mesma amostra. Resultado: um achado genuinamente relevante
ficou invisível tanto pro pipeline quanto pro laudo — quase virou um falso negativo relatado como
conclusão.

**O que isso revela:** o `qc_rescued` (`pipeline.py`) hoje só dispara pra
`clinvar_significance in {Pathogenic, Likely pathogenic}`. Isso deixa de fora exatamente o caso
mais comum na prática — variante rara, QC-dropped, **VUS bem revisado (2+ estrelas) numa condição
que bate com o fenótipo do paciente** — que é justamente o tipo de achado que mais precisa de
visibilidade pra decisão de repetir/confirmar a chamada.

**Proposta pro produto:**
- Ampliar o critério de `qc_rescued` pra também disparar quando: `clinvar_review_status` tem ≥2
  estrelas **E** (`clinvar_significance` é qualquer coisa não-Benign **OU** a condição do ClinVar
  bate com o HPO do paciente via o mesmo scorer do PP4).
- Mais geral ainda: rodar o rescue check pra **qualquer variante rara (abaixo do teto de PM2) que
  caiu no QC E está em um gene com HPO match ≥ cutoff**, independente de já estar no ClinVar —
  esse é o caso mais amplo (novo/não catalogado, mas no gene certo, no fenótipo certo).
- No laudo, dar a esses achados o mesmo tratamento visual dos QC-rescued existentes (card com
  alerta, campos de qualidade completos, nota de "confirmar antes de descartar") — não os deixar
  invisíveis só porque não bateram tier de candidato.

**Lição de processo (não só de produto):** antes de declarar "nenhum candidato" pra um gene/
fenótipo específico, sempre checar o VCF anotado bruto por esse gene, não só a lista final de
`candidates`/`classifications` já filtrada — a ausência na lista filtrada pode significar "não
existe variante" ou pode significar "existe mas foi descartada num estágio anterior", e são coisas
muito diferentes pra reportar.

---

*Contexto: os cinco itens foram validados manualmente numa sessão de análise da amostra
Genoma_low_Sample (WGS, BaseSpace project "Genoma_low"), com o resultado publicado em
`~/GBBGENOME/vcf2report_out/Genoma_low_Sample_laudo.html`. Os scripts/lógica usados ali são um bom
ponto de partida pra implementar isso no pipeline principal.*
