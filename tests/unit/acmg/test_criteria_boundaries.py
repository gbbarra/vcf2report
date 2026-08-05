"""As fronteiras numéricas de `acmg/criteria.py`, e o que decide exatamente EM CIMA delas.

O §4.1 do CLAUDE.md exige, para toda função nova, "o limite (`boundary`) de qualquer
comparação numérica". A primeira rodada de mutação sobre `criteria.py` mostrou que o módulo
inteiro estava sem essa camada: com **100% de cobertura de branches** — que é justamente o
que o §4.3 diz não valer como argumento —, dez mutantes sobreviveram à suíte COMPLETA, e
seis deles são um `>=` virando `>` ou um `<=` virando `<`.

Cada teste aqui foi confirmado dos dois lados: passa no fonte original e falha no mutante
correspondente. Sem essa segunda metade, um teste de fronteira que não toca a fronteira
passa despercebido.

Os mutantes que estes testes matam:

  _clinvar_reviewed  "practice guideline" -> caixa alta
  _pm1_signals       n_residues  >= HOTSPOT_MIN_RESIDUES   -> >
  _pm1_signals       enrichment  >= HOTSPOT_MIN_ENRICHMENT -> >
  _insilico_direction  revel >= REVEL_PATHOGENIC -> >
  _insilico_direction  cadd  >= CADD_PATHOGENIC  -> >
  _insilico_direction  revel <= REVEL_BENIGN     -> <
  _insilico_direction  cadd  <= CADD_BENIGN      -> <
  _benign_af         braz > faf -> >=
  _pm5_strength      n_other default 1 -> 2
  _pm5_strength      stars   default 0 -> 1

Três outros candidatos foram descartados aqui porque a suíte completa já os mata em
tests/test_audit_fixes.py (or->and em `_clinvar_reviewed` e em `_insilico_direction`, e a
caixa de "reviewed by expert"). Eram sobreviventes da SELEÇÃO do mutmut, não do código —
terceira rodada seguida em que isso acontece. Os três testes que os matam entraram na
seleção do setup.cfg por node ID.
"""

from __future__ import annotations

import pytest

from vcf2report.acmg import criteria as C
from vcf2report.acmg.engine import evaluate_criteria
from vcf2report.annotate import clinvar_residue as residue
from vcf2report.models import Annotation, Variant


def _v(consequence: str = "missense_variant") -> Variant:
    return Variant(
        chrom="1",
        pos=100,
        ref="C",
        alt="T",
        gene="TESTG",
        consequence=consequence,
        hgvs_p="p.Arg123Cys",
    )


def _crit(code: str, a: Annotation, v: Variant | None = None):
    return next(c for c in evaluate_criteria(v or _v(), a) if c.code == code)


# --------------------------------------------------------------- in-silico: os quatro limiares


@pytest.mark.parametrize(
    "campo, limiar, direcao_esperada",
    [
        ("revel", C.REVEL_PATHOGENIC, "pathogenic"),
        ("cadd_phred", C.CADD_PATHOGENIC, "pathogenic"),
        ("revel", C.REVEL_BENIGN, "benign"),
        ("cadd_phred", C.CADD_BENIGN, "benign"),
    ],
    ids=["revel-patogenico", "cadd-patogenico", "revel-benigno", "cadd-benigno"],
)
def test_um_preditor_exatamente_no_limiar_ja_conta(campo, limiar, direcao_esperada):
    """Os quatro limiares in-silico são inclusivos, e é o valor EXATO que prova isso.

    `_insilico_direction` compara `>= REVEL_PATHOGENIC` (0,70), `>= CADD_PATHOGENIC` (20,0),
    `<= REVEL_BENIGN` (0,15) e `<= CADD_BENIGN` (10,0). Trocar qualquer um dos quatro por
    estrito muda a resposta em um único ponto — o limiar — e a suíte inteira não notava:
    um REVEL de exatamente 0,70 deixaria de sustentar PP3, e um CADD de exatamente 10,0
    deixaria de sustentar BP4. É a diferença entre um in-silico contar e não contar, no
    valor que os autores dos escores escolheram como corte.
    """
    assert C._insilico_direction(Annotation(**{campo: limiar})) == direcao_esperada


@pytest.mark.parametrize(
    "campo, limiar, delta",
    [
        ("revel", C.REVEL_PATHOGENIC, -0.01),
        ("cadd_phred", C.CADD_PATHOGENIC, -0.1),
        ("revel", C.REVEL_BENIGN, +0.01),
        ("cadd_phred", C.CADD_BENIGN, +0.1),
    ],
    ids=["revel-patogenico", "cadd-patogenico", "revel-benigno", "cadd-benigno"],
)
def test_do_lado_de_dentro_do_limiar_nao_conta(campo, limiar, delta):
    """O outro lado da mesma fronteira: sem isto, o teste acima passaria com a comparação
    substituída por `is not None`, que aceita qualquer valor."""
    assert C._insilico_direction(Annotation(**{campo: limiar + delta})) is None


# ------------------------------------------------------------------- PM1: os dois limiares


def _hotspot(n_residues: int, enrichment: float) -> Annotation:
    return Annotation(
        clinvar_hotspot={
            "n_residues": n_residues,
            "n_changes": 9,
            "enrichment": enrichment,
        },
        gene_missense_tolerant=False,
    )


@pytest.mark.parametrize(
    "n_residues, dense_esperado",
    [
        (residue.HOTSPOT_MIN_RESIDUES - 1, False),
        (residue.HOTSPOT_MIN_RESIDUES, True),
    ],
    ids=["abaixo", "no-limiar"],
)
def test_o_hotspot_conta_a_partir_do_numero_minimo_de_residuos(
    n_residues, dense_esperado
):
    """`n_residues >= HOTSPOT_MIN_RESIDUES` (3). No limiar exato o hotspot já é denso, e é o
    que faz PM1 disparar — com `>` no lugar de `>=`, todo hotspot de exatamente 3 resíduos
    deixa de existir para o motor, silenciosamente."""
    s = C._pm1_signals(_v(), _hotspot(n_residues, 5.0))

    assert s["dense"] is dense_esperado
    assert s["fires"] is dense_esperado, "o veredicto de PM1 tem de seguir o sinal"


@pytest.mark.parametrize(
    "enrichment, enriched_esperado",
    [
        (residue.HOTSPOT_MIN_ENRICHMENT - 0.1, False),
        (residue.HOTSPOT_MIN_ENRICHMENT, True),
    ],
    ids=["abaixo", "no-limiar"],
)
def test_o_hotspot_conta_a_partir_do_enriquecimento_minimo(
    enrichment, enriched_esperado
):
    """`enrichment >= HOTSPOT_MIN_ENRICHMENT` (2,0) — a densidade tem de ser materialmente
    maior que a linha de base do próprio gene, e "materialmente" começa NO 2,0."""
    s = C._pm1_signals(_v(), _hotspot(9, enrichment))

    assert s["enriched"] is enriched_esperado
    assert s["fires"] is enriched_esperado


# ------------------------------------------------------- _benign_af: o empate entre as fontes


def test_coorte_local_igual_ao_faf95_nao_e_reportada_como_acima_dele():
    """`_benign_af` só troca a base quando a coorte local é ESTRITAMENTE maior.

    No empate o valor devolvido é o mesmo dos dois jeitos — 0,06 é 0,06 — então nenhum tier
    se move e nada na suíte reclamava. O que muda é a frase: com `>=` no lugar de `>`, o
    laudo passa a dizer "local cohort AF — **above** the gnomAD filtering AF" sobre um
    número que é igual, não maior. É a invariante de honestidade aplicada à coluna de
    evidência: a base citada tem de descrever a comparação que de fato aconteceu.
    """
    a = Annotation(gnomad_faf95=0.06, local_cohort_af=0.06)

    af, base = C._benign_af(a)
    ba1 = _crit("BA1", a)

    assert af == 0.06
    assert base == "gnomAD filtering AF (faf95, grpmax)"
    assert ba1.met is True, "o empate acima de 5% continua sendo BA1"
    assert "above the gnomAD filtering AF" not in ba1.reasoning
    assert ba1.evidence["basis"] == base


def test_coorte_local_acima_do_faf95_e_que_troca_a_base():
    """O outro lado: estritamente maior troca, e a frase passa a ser verdadeira."""
    af, base = C._benign_af(Annotation(gnomad_faf95=0.04, local_cohort_af=0.05))

    assert af == 0.05
    assert base == "local cohort AF — above the gnomAD filtering AF"


# ---------------------------------------------------- _pm5_strength: os defaults dos ausentes


@pytest.mark.parametrize(
    "match, forca_esperada",
    [
        ({"n_changes": 3}, "supporting"),
        ({"n_other": None, "stars": None}, "supporting"),
        ({"stars": 2}, "moderate"),
        ({"n_other": 2}, "moderate"),
        ({"n_other": 2, "stars": 2}, "strong"),
    ],
    ids=["ambos-ausentes", "ambos-nulos", "so-estrelas", "so-mudancas", "os-dois"],
)
def test_pm5_trata_campo_ausente_como_o_caso_mais_fraco(match, forca_esperada):
    """Um índice de resíduos que não informa `n_other` ou `stars` cai no default declarado
    no docstring: uma única outra mudança (`n_other=1`) apoiada em zero estrelas — o caso
    Supporting.

    Os dois defaults sobreviviam à suíte. Subir o de `n_other` para 2 promoveria a
    Moderate toda correspondência que não informa a contagem, e subir o de `stars` para 1
    faria o mesmo — em ambos os casos elevando a força de PM5 com base em dado que o índice
    não forneceu. É exatamente "ausência de dado virando evidência", na escolha de força.
    """
    assert C._pm5_strength(match) == forca_esperada


def test_pm5_sem_correspondencia_alguma_nao_tem_forca():
    """`None` e `{}` não são "o caso mais fraco" — são "não houve correspondência"."""
    assert C._pm5_strength(None) is None
    assert C._pm5_strength({}) is None


# ------------------------------------------------- _clinvar_reviewed: as três portas de 1★+


@pytest.mark.parametrize(
    "review_status, revisado",
    [
        ("criteria provided, single submitter", True),
        ("reviewed by expert panel", True),
        ("practice guideline", True),
        ("no assertion criteria provided", False),
        ("no assertion provided", False),
    ],
)
def test_as_tres_portas_de_revisao_do_clinvar_e_a_que_fica_de_fora(
    review_status, revisado
):
    """`_clinvar_reviewed` tem três aceitações independentes, e a de mais alto nível —
    "practice guideline", o 4★ do ClinVar — não era exercitada por teste nenhum.

    O gate decide se PP5 dispara e se PS1/PM5 se recolhem, então perdê-lo significa tratar
    a asserção mais bem revisada do ClinVar como se não tivesse critério algum. A linha de
    0★ que CONTÉM "criteria provided" fica de fora pelo `startswith`, e continua fora aqui.
    """
    assert (
        C._clinvar_reviewed(Annotation(clinvar_review_status=review_status)) is revisado
    )


def test_practice_guideline_chega_ate_o_pp5():
    """Pela API pública, para provar que o gate não é decoração interna."""
    a = Annotation(
        clinvar_significance="Pathogenic", clinvar_review_status="practice guideline"
    )

    assert _crit("PP5", a).met is True
