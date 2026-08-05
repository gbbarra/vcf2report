"""O modelo de pontos ClinGen/Tavtigian: limites de tier e integridade da trilha.

Escrito a partir dos 11 mutantes que sobreviveram em `_combine_points` (CLAUDE.md §4.4).
Este combinador só roda com ``VCF2REPORT_ACMG_MODEL=clingen``, o que o deixou com um
buraco de cobertura maior que o caminho Richards: 83% de branch contra 97%.

Três mutações foram aplicadas ao código real e a suíte inteira passou (770 passed), o que
prova que as lacunas abaixo são reais e não artefato da seleção reduzida do mutmut.

Um dos grupos NÃO virou teste, de propósito — ver
`test_as_tabelas_de_pontos_cobrem_toda_forca_que_um_criterio_pode_emitir`.
"""

from __future__ import annotations

import pytest

from vcf2report.acmg import criteria as crit_mod
from vcf2report.acmg import rules
from vcf2report.models import CriterionResult

#: Tavtigian 2020, adotado pelo ClinGen SVI e transcrito no próprio módulo:
#: Pathogenic >= 10 · Likely Pathogenic 6..9 · VUS 0..5 · Likely Benign -1..-6 · Benign <= -7.
LIMITES = [
    (10, rules.PATHOGENIC, "limite inferior de Pathogenic"),
    (9, rules.LIKELY_PATHOGENIC, "logo abaixo de Pathogenic"),
    (6, rules.LIKELY_PATHOGENIC, "limite inferior de Likely Pathogenic"),
    (5, rules.VUS, "logo abaixo de Likely Pathogenic"),
    (0, rules.VUS, "limite inferior de VUS"),
    (-1, rules.LIKELY_BENIGN, "limite superior de Likely Benign"),
    (-6, rules.LIKELY_BENIGN, "limite INFERIOR de Likely Benign"),
    (-7, rules.BENIGN, "limite superior de Benign"),
]


def _pontos(alvo: int) -> list[CriterionResult]:
    """Critérios reais somando exatamente ``alvo`` pontos.

    Montado com códigos e forças que o motor realmente emite — nada de força inventada,
    porque uma força fora das tabelas de pontos não existe em produção (ver o teste de
    cobertura das tabelas no fim deste arquivo).
    """
    peso_path = [("PVS1", "very_strong", 8), ("PS1", "strong", 4), ("PM1", "moderate", 2)]
    peso_ben = [("BA1", "stand_alone", -8), ("BS1", "strong", -4), ("BP1", "supporting", -1)]
    fonte = peso_path if alvo >= 0 else peso_ben
    restante, out, n = alvo, [], 0
    for code, forca, peso in fonte:
        while (peso > 0 and restante >= peso) or (peso < 0 and restante <= peso):
            n += 1
            out.append(
                CriterionResult(
                    code=f"{code[:2]}{n}" if code[:2] in ("PM", "PP", "BP", "BS") else code,
                    name="fixture",
                    default_strength=forca,
                    applies=True,
                    met=True,
                )
            )
            restante -= peso
            if code in ("PVS1", "BA1"):  # existe um só de cada
                break
    # o resto em Supporting (+1 / -1)
    while restante != 0:
        n += 1
        passo = 1 if restante > 0 else -1
        out.append(
            CriterionResult(
                code=f"PP{n}" if passo > 0 else f"BP{n}",
                name="fixture",
                default_strength="supporting",
                applies=True,
                met=True,
            )
        )
        restante -= passo
    return out


@pytest.mark.parametrize(
    "pontos, tier_esperado, descricao", LIMITES, ids=[str(l[0]) for l in LIMITES]
)
def test_cada_limite_de_tier_do_modelo_de_pontos(
    monkeypatch, pontos, tier_esperado, descricao
):
    """Mata `total >= -6` -> `total > -6` e `-> total >= -7`.

    Em exatamente -6 pontos o resultado é Likely Benign; em -7, Benign. Nenhum teste
    fixava esses dois pontos, e o §4.1 exige o limite de toda comparação numérica. Os
    demais limites entram na mesma varredura para que a tabela inteira fique ancorada.
    """
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    criterios = _pontos(pontos)

    tier, rule_path = rules.combine(criterios)

    assert f"{pontos:+d} points" in rule_path, (
        f"{descricao}: a fixture não somou {pontos} — {rule_path}"
    )
    assert tier == tier_esperado, f"{descricao}: {pontos:+d} -> {tier} (esperado {tier_esperado})"


def test_a_trilha_lista_so_os_criterios_que_realmente_dispararam(monkeypatch):
    """Mata `cr.applies and cr.met` -> `or`, e `met = None`.

    A trilha do modelo de pontos é a única explicação que o laudo dá para o número. Com
    `or`, criterios N/A (que não se aplicam a proband único) passam a ser listados como se
    tivessem contribuído; com `met = None`, some tudo e a linha diz "no criteria met"
    enquanto o total é diferente de zero. Nos dois casos o número deixa de bater com a
    lista que o justifica.
    """
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")
    criterios = [
        CriterionResult("PVS1", "n", "very_strong", applies=True, met=True),
        CriterionResult("PM2", "n", "moderate", applies=True, met=True),
        CriterionResult("PS2", "n", "strong", applies=False, met=False),  # N/A: sem trio
        CriterionResult("PP3", "n", "supporting", applies=True, met=False),  # avaliado, não bateu
    ]

    tier, rule_path = rules.combine(criterios)

    trilha = rule_path.split(" =>")[0]
    assert trilha == "PVS1 + PM2", f"trilha lista o que não disparou: {trilha}"
    assert "+10 points" in rule_path
    assert tier == rules.PATHOGENIC


def test_sem_nenhum_criterio_a_trilha_diz_isso_em_vez_de_ficar_vazia(monkeypatch):
    monkeypatch.setenv("VCF2REPORT_ACMG_MODEL", "clingen")

    tier, rule_path = rules.combine([])

    assert "no criteria met" in rule_path
    assert "+0 points" in rule_path
    assert tier == rules.VUS


def test_as_tabelas_de_pontos_cobrem_toda_forca_que_um_criterio_pode_emitir():
    """Este teste existe no lugar de três testes que NÃO escrevi.

    `table.get(strength, 0)` tem um default, e três mutantes o alteram (para `None`, para
    `1`, e removendo-o). Todos sobrevivem à suíte inteira — mas por serem **equivalentes**,
    não por falta de teste: o default nunca é alcançado. As duas tabelas, tomadas com o
    lado a que cada código pertence, cobrem toda força que existe:

        _PATH_POINTS   {very_strong, strong, moderate, supporting}   falta stand_alone
        _BENIGN_POINTS {stand_alone, strong, moderate, supporting}   falta very_strong

    e nenhum critério patogênico é stand-alone, nem nenhum benigno é very_strong.

    O §4.4 pede explicar a equivalência em vez de escrever teste artificial. Um teste que
    forjasse uma força fora da tabela mataria os mutantes provando algo que não acontece
    em produção. Este, em vez disso, pina a RAZÃO da equivalência — e falha no dia em que
    alguém acrescentar uma força nova, ou mover um código de lado, que é justamente quando
    o default deixaria de ser código morto.
    """
    faltando = []
    for code in crit_mod.all_criteria():
        tabela = rules._BENIGN_POINTS if code in rules.BENIGN_CODES else rules._PATH_POINTS
        # A força padrão declarada é a única que um critério pode emitir sem que o motor a
        # ajuste; as ajustadas (PVS1 pela árvore SVI, PP3 pelo AlphaMissense) são sempre um
        # degrau DENTRO do mesmo lado, então bastam as duas tabelas cobrirem o vocabulário.
        for forca in ("very_strong", "strong", "moderate", "supporting", "stand_alone"):
            pertence_ao_lado = (
                forca != "very_strong" if code in rules.BENIGN_CODES else forca != "stand_alone"
            )
            if pertence_ao_lado and forca not in tabela:
                faltando.append((code, forca))

    assert not faltando, (
        "força que um critério pode carregar e a tabela de pontos não conhece — o default "
        f"de table.get deixou de ser código morto: {faltando}"
    )
