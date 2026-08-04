"""Limites numéricos e precedência de força no combinador da Tabela 5.

Escritos a partir de mutantes que SOBREVIVERAM à suíte completa (CLAUDE.md §4.4). Cada
teste aqui nomeia, no docstring, a mutação que ele mata — porque um teste de limite sem
essa referência é indistinguível de um teste decorativo.

Verificação feita antes de escrever: as duas mutações abaixo foram aplicadas ao código
real e a suíte inteira passou (736 passed), o que prova que a lacuna é real e não
artefato da seleção reduzida usada pelo mutmut.
"""
from __future__ import annotations

import pytest

from vcf2report.acmg import rules
from vcf2report.models import CriterionResult


def _met(code: str, default_strength: str, applied_strength: str | None = None):
    """Um critério ATENDIDO, do jeito que o motor o entrega ao combinador."""
    return CriterionResult(
        code=code, name=f"{code} (fixture)", default_strength=default_strength,
        applies=True, met=True, applied_strength=applied_strength,
    )


def _supporting(n: int) -> list[CriterionResult]:
    """`n` critérios de suporte patogênico distintos (códigos reais, PP1..PP5)."""
    return [_met(f"PP{i}", "supporting") for i in range(1, n + 1)]


# --------------------------------------------------------- PATH-3: 1 Strong + Moderate/Supporting

@pytest.mark.parametrize("n_supporting, esperado", [
    (3, rules.LIKELY_PATHOGENIC),   # abaixo do limite
    (4, rules.PATHOGENIC),          # NO limite — é este ponto que o mutante move
    (5, rules.PATHOGENIC),          # acima
])
def test_path3_com_um_moderado_dispara_a_partir_de_quatro_supporting(n_supporting, esperado):
    """Mata `pp >= 4` -> `pp > 4` em PATH-3.

    Richards Tabela 5, PATH-3: 1 Strong + 1 Moderate + **≥4** Supporting é Patogênica.
    Com `pp > 4` o caso de exatamente 4 cai para Likely Pathogenic — uma variante rebaixada
    de Patogênica num limite que a diretriz define. Nenhum teste da suíte fixava pp == 4.
    """
    criteria = [_met("PS1", "strong"), _met("PM1", "moderate")] + _supporting(n_supporting)

    tier, rule_path = rules.combine(criteria)

    assert tier == esperado
    if esperado == rules.PATHOGENIC:
        assert rule_path.startswith("PS1 + PM1")
        assert "PATH-3" in rule_path


def test_path3_exige_o_moderado_e_nao_so_os_supporting():
    """Mata `pm == 1` -> `pm != 1` em PATH-3.

    Sem nenhum Moderado, 1 Strong + 4 Supporting **não** é Patogênica — a cláusula pede
    `pm == 1`. Com `pm != 1` o zero satisfaz a condição e a variante sobe indevidamente.
    """
    criteria = [_met("PS1", "strong")] + _supporting(4)

    tier, _ = rules.combine(criteria)

    assert tier == rules.LIKELY_PATHOGENIC


# --------------------------------------------- precedência da força APLICADA sobre a padrão

def test_a_forca_aplicada_prevalece_sobre_a_padrao_no_veto_por_evidencia_decisiva():
    """Mata `cr.applied_strength or cr.default_strength` -> `and` em _discarded_decisive.

    A árvore SVI rebaixa PVS1 para `strong` num nulo de último éxon. O veto por evidência
    decisiva precisa ler a força APLICADA; com `and`, um critério que tem `applied_strength`
    definido devolve a `default_strength`, e o veto passa a julgar uma força que não é a que
    o motor atribuiu.

    Montagem: BA1 (stand-alone benigno) vence, e do lado patogênico perdedor há um PVS1
    REBAIXADO a `moderate`. `moderate` não é decisivo, então NÃO há conflito e a chamada
    benigna deve permanecer. Sob o mutante, a força lida vira `very_strong` (a padrão), que
    é decisiva, e o resultado viraria VUS.
    """
    criteria = [
        _met("BA1", "stand_alone"),
        _met("PVS1", "very_strong", applied_strength="moderate"),
    ]

    tier, _ = rules.combine(criteria)

    assert tier == rules.BENIGN


def test_um_pvs1_nao_rebaixado_do_lado_perdedor_forca_o_conflito():
    """O par do teste acima: sem rebaixamento, PVS1 É decisivo e o conflito deve aparecer.

    Sem este, o teste anterior passaria por acidente caso o veto parasse de funcionar por
    completo — a asserção seria satisfeita pelo motivo errado.
    """
    criteria = [
        _met("BA1", "stand_alone"),
        _met("PVS1", "very_strong"),
    ]

    tier, rule_path = rules.combine(criteria)

    assert tier == rules.VUS
    assert "conflicting" in rule_path


# ----------------------------------------------- o rótulo deve nomear a classe que disparou

def test_o_rotulo_de_lp1_nomeia_moderate_quando_o_moderado_e_que_disparou():
    """Mata `pm >= 1` -> `pm > 1` na ESCOLHA DO RÓTULO de LP-1.

    O tier continua Likely Pathogenic sob o mutante — o que muda é a trilha de auditoria,
    que passaria a dizer "PVS1 + Supporting" quando quem disparou foi um Moderado. O laudo
    é auditável: a linha que explica a decisão precisa nomear a evidência correta.
    """
    criteria = [_met("PVS1", "very_strong"), _met("PM2", "moderate")]

    tier, rule_path = rules.combine(criteria)

    assert tier == rules.LIKELY_PATHOGENIC
    assert "LP-1 (PVS1 + Moderate)" in rule_path


def test_o_rotulo_de_lp1_nomeia_supporting_quando_nao_ha_moderado():
    """O outro lado do mesmo `if`: só com Supporting, o rótulo tem de dizer Supporting."""
    criteria = [_met("PVS1", "very_strong"), _met("PP3", "supporting")]

    tier, rule_path = rules.combine(criteria)

    assert tier == rules.LIKELY_PATHOGENIC
    assert "LP-1 (PVS1 + Supporting)" in rule_path
