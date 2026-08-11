# CLAUDE.md

Instruções permanentes para o Claude neste repositório.
Leia este arquivo antes de qualquer alteração de código.

---

## 1. Regras invioláveis

Estas regras têm precedência sobre qualquer pedido meu feito no calor do momento.
Se um pedido meu conflitar com elas, **pare e me avise** em vez de obedecer.

1. **Nunca altere um teste existente para fazê-lo passar.**
   Se um teste falha, o código está errado — ou o teste está errado e você deve me
   perguntar antes de tocar nele. Editar a asserção para casar com o output é proibido.

2. **Nunca remova, comente, pule (`skip`/`xfail`) ou afrouxe um teste** sem
   autorização explícita minha, mensagem por mensagem.

3. **Nunca reduza um limiar de qualidade** (cobertura mínima, complexidade máxima,
   regras de lint) para fazer o pipeline passar.

4. **Não declare "pronto" ou "funcionando" sem ter executado os testes** e colado a
   saída real. Se você não rodou, diga que não rodou.

5. **Teste antes de código.** Para qualquer comportamento novo: escreva o teste,
   mostre-o falhando, só então implemente. Sem exceção para "mudança pequena".

6. **Nenhum `assert` decorativo.** São proibidos como asserção única:
   `assert result is not None`, `assert result`, `assert len(x) > 0`,
   `assert isinstance(x, dict)`. Toda asserção verifica um valor concreto esperado.

7. **Nada de mock do que está sendo testado.** Mock é para I/O externo (rede, banco,
   sistema de arquivos, API paga). Se você precisou mockar a lógica de negócio para o
   teste passar, o desenho está errado — me avise.

---

## 2. Fluxo de trabalho padrão

Para cada tarefa, siga esta ordem e **pare para confirmação** entre as etapas 2 e 3:

1. **Entender** — releia o código existente antes de escrever. Não presuma a API.
2. **Especificar** — descreva em uma frase o comportamento esperado e os casos de
   borda que vai cobrir. Espere meu OK.
3. **Teste vermelho** — escreva o(s) teste(s) e rode. Cole a saída mostrando a falha.
4. **Implementar** — o mínimo de código para passar. Nada de funcionalidade extra
   "que pode ser útil depois".
5. **Verde** — rode a suíte inteira, não só o teste novo. Cole a saída.
6. **Refatorar** — só com a suíte verde, e rodando novamente ao final.
7. **Gate** — rode lint, tipos e cobertura. Reporte os números.

---

## 3. Comandos do projeto

Ajustados à realidade deste repositório em 2026-08-04. Linhas marcadas **(pendente)**
dependem de ferramentas ainda não declaradas em `pyproject.toml` — ver §8.

```bash
# suíte completa (funciona hoje; ~5 min com os stores montados)
python3 -m pytest tests/ -q

# suíte + cobertura de branches                              (pendente: pytest-cov)
python3 -m pytest tests/ -q --cov=src --cov-branch --cov-report=term-missing --cov-fail-under=80

# um teste específico durante o desenvolvimento
python3 -m pytest tests/test_acmg_rules.py::test_nome -x -vv

# lint e formatação                                          (pendente: ruff)
python3 -m ruff check . && python3 -m ruff format --check .

# tipagem estática                                           (pendente: mypy)
python3 -m mypy src/

# testes de aceitação (BDD)                                  (pendente: pytest-bdd)
python3 -m pytest tests/features/ -v

# métricas de complexidade e manutenibilidade                (pendente: radon)
python3 -m radon cc -s -a src/ && python3 -m radon mi src/

# teste de mutação (lento — sob demanda)                     (pendente: mutmut)
python3 -m mutmut run --paths-to-mutate src/vcf2report/acmg/
python3 -m mutmut results

# gate dos stores de anotação (obrigatório antes de medir o benchmark)
python3 scripts/check_stores.py --gate

# benchmark de 200 exomas, comparando contra o baseline versionado
python3 scripts/run_benchmark.py --annotated <bench>/realistic_annotated \
    --bench <bench> --jobs 4 --withhold-clinvar \
    --out after.tsv --compare data/benchmark/hpo-spiked-exomes-baseline.tsv
```

---

## 4. Padrões por tipo de teste

### 4.1 Testes unitários
- Ficam em `tests/unit/`, espelhando a estrutura de `src/`.
- Nome descreve o comportamento, não a função:
  `test_retorna_erro_quando_campo_info_ausente`, não `test_parse_1`.
- Estrutura **Arrange / Act / Assert**, com linha em branco separando os blocos.
- Um comportamento por teste. Vários `assert` só se verificarem o mesmo comportamento.
- Sem rede, sem banco, sem I/O real. Use `tmp_path` do pytest para arquivos.
- Para toda função nova, cubra no mínimo: caso feliz, entrada vazia/nula, entrada
  malformada, e o limite (`boundary`) de qualquer comparação numérica.
- Erros esperados se testam com `pytest.raises(TipoDoErro, match="trecho da mensagem")`.
- Use `@pytest.mark.parametrize` em vez de copiar e colar variações do mesmo teste.

> **Exceção deste repositório:** existem testes que consultam fontes reais de propósito
> (`tests/test_gnomad_live.py` faz range-request no bucket público do gnomAD). Eles
> **não** são unitários e não devem ser mockados — a decisão de projeto é não criar
> artefato fictício. Ficam marcados `@pytest.mark.integration` e fora de `tests/unit/`.

### 4.2 Testes Gherkin / BDD
- Reservados para **regras de negócio que um revisor não-programador precisa validar**.
  Não use BDD para função utilitária — isso é teste unitário.
- `.feature` em `tests/features/`, escrito em português, com `# language: pt`.
- Steps em `tests/features/steps/`, usando `pytest-bdd`.
- O `.feature` descreve **o quê**, nunca **como**. Proibido citar nome de função,
  classe, endpoint ou estrutura de dados no cenário.
- Escreva o `.feature` primeiro, mostre para mim, e só depois implemente os steps.

```gherkin
# language: pt
Funcionalidade: Classificação de variantes
  Cenário: Critérios de patogenicidade forte e moderado
    Dado uma variante com o critério PVS1 atribuído
    E o critério PM2 atribuído
    Quando a classificação for calculada
    Então o resultado deve ser "Provavelmente Patogênica"
```

> **Nota de idioma:** os `.feature` são documentação interna e ficam em português.
> O **laudo gerado é sempre em inglês** — nenhum cenário deve exigir saída em português.

### 4.3 Cobertura
- Sempre com `--cov-branch`. Cobertura de linha sozinha esconde `if` não testado.
- Mínimo do projeto: **80%**. Módulos de lógica de decisão: **95%** (lista no §7).
- Cobertura é métrica de **ausência** — mostra o que não foi testado, não valida o que
  foi. Nunca use "100% de cobertura" como argumento de que está correto.
- Ao priorizar, ataque branches de decisão. Ignore boilerplate, `__init__.py`,
  `__repr__` e blocos `if TYPE_CHECKING`.

### 4.4 Teste de mutação
- Rode sob demanda (antes de release ou ao fechar um módulo), nunca no loop de
  desenvolvimento — é lento.
- Para cada mutante sobrevivente: mostre o diff do mutante e escreva o teste que o mata.
- Mutante sobrevivente em código de decisão clínica ou de cálculo é **bloqueante**.
- Se um mutante for genuinamente equivalente (não altera comportamento observável),
  explique por quê em vez de escrever teste artificial.

### 4.5 Métricas de qualidade
- Complexidade ciclomática máxima por função: **10**.
- Função acima de 50 linhas ou com mais de 5 parâmetros: proponha refatoração.
- Ao encontrar violação, **proponha o plano primeiro** e espere meu OK. Não refatore
  código que eu não pedi para tocar.

### 4.6 Quality gates (CI)
O pipeline deve falhar (exit code ≠ 0) se qualquer etapa não passar:

- [ ] `ruff check` sem erros
- [ ] `ruff format --check` sem diferenças
- [ ] `mypy` sem erros
- [x] `pytest` com toda a suíte verde
- [x] cobertura de branches ≥ 80%
- [ ] `bandit` sem achado de severidade alta
- [ ] nenhum segredo commitado (`detect-secrets` ou equivalente)

Espelhe os mesmos gates em `pre-commit` para pegar antes do push.

---

## 5. O que NÃO fazer

- ❌ Escrever teste depois do código "para bater a cobertura".
- ❌ Ajustar o teste até passar.
- ❌ `try/except: pass` para silenciar erro que o teste expôs.
- ❌ Adicionar dependência nova sem me perguntar.
- ❌ Alterar arquivo de configuração (`pyproject.toml`, CI, `.pre-commit-config.yaml`)
  no meio de uma tarefa de código, sem avisar.
- ❌ Reescrever módulo inteiro quando eu pedi uma correção pontual.
- ❌ Resumir a saída dos testes. Cole a saída real, inclusive as falhas.

---

## 6. Comunicação

- Responda em português.
- Se um requisito estiver ambíguo, **pergunte antes de implementar** — não escolha por mim.
- Se você discordar de uma decisão minha, diga, com o motivo técnico. Não obedeça em
  silêncio a algo que você acha errado.
- Ao terminar, reporte: testes que passaram/falharam, cobertura antes e depois, e o
  que ficou de fora.

---

## 7. Contexto do projeto

- **Domínio:** interpretação clínica de exoma. VCF single-proband GRCh38 + termos HPO →
  laudo ACMG/AMP 2015 auditável (com refinamentos ClinGen SVI). Roda **local e offline**.
  O laudo é **rascunho para revisão por profissional**, não dispositivo diagnóstico.

- **Stack:** Python 3.10–3.12 · pytest (+hypothesis) · DuckDB/Parquet para os stores de
  anotação · pysam (extra `tabix`) · Jinja2 opcional para o laudo · MCP opcional.
  Sem framework web, sem banco de dados, sem serviço externo em tempo de execução.

- **Estrutura de diretórios:**
  ```
  src/vcf2report/
    acmg/        criteria.py (28 critérios) · rules.py (Tabela 5) · engine.py
    annotate/    gnomad{,_parquet,_local,_remote}.py · clinvar*.py · alphamissense*.py
                 hpo.py · local_cohort.py · inheritance.py · cache.py
    report/      assemble.py (baldes + conclusão) · render.py · vus_triage.py · explore.py
    vcf/         parse.py · annparse.py · filter.py · qc.py · seqqc.py
    pipeline.py · models.py · config.py · stores.py · concordance.py · demo.py
  tests/         59 arquivos, ~741 testes (organizados por comportamento/defeito)
  scripts/       build/fetch dos stores · run_benchmark.py · sweep_cohort.py · check_stores.py
  data/          stores parquet (git-ignored) · example/ · benchmark/ · hpo/ · constraint/
  templates/     report.md.j2
  ```

- **Módulos críticos** (exigem cobertura ≥ 95% e mutação limpa):
  | módulo | por quê |
  |---|---|
  | `acmg/criteria.py` | decide cada um dos 28 critérios ACMG |
  | `acmg/rules.py` | combina critérios no tier final (Tabela 5 / pontos) |
  | `acmg/engine.py` | orquestra a classificação |
  | `report/assemble.py` | roteia achados nos baldes e escreve a conclusão |
  | `report/vus_triage.py` | decide o que é VUS provável-patogênica |
  | `vcf/filter.py` | decide o que sobrevive à triagem — um descarte aqui é invisível |
  | `annotate/gnomad_parquet.py` | decide se uma frequência é "ausente", "desconhecida" ou um valor |
  | `pipeline.py` | ordem das etapas e propagação de estado entre elas |

- **Restrições regulatórias / de dado sensível:**
  - Dado de paciente **nunca sai da máquina**. Nenhuma chamada de rede com conteúdo do VCF.
  - Todo laudo carrega o aviso **DRAFT — not for clinical use** e é **sempre em inglês**.
  - **Invariante de honestidade:** ausência de dado nunca pode ser reportada como evidência.
    "Não consultamos" ≠ "consultamos e não achou". Vale para PM2/BA1/BS1/BS2 e para toda
    coluna de evidência do laudo.
  - Licenças dos stores: gnomAD **ODbL-1.0** (atribuição + share-alike) · AlphaMissense
    **CC BY 4.0** (não redistribuído por este projeto) · ClinVar **domínio público**.
    O `_manifest.json` de cada store carimba a licença — não altere sem verificar a fonte.

---

## 8. Estado das ferramentas (2026-08-04)

Instalado e funcionando: `pytest`, `hypothesis`.

**Ainda não declarado em `pyproject.toml`** — medido neste ambiente para levantar a linha
de base, mas não adicionado como dependência (§5 proíbe sem autorização):
`ruff`, `mypy`, `pytest-cov`, `radon`, `bandit`, `pytest-bdd`, `detect-secrets`, `mutmut`.
