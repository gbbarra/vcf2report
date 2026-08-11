"""Precedência entre as regras patogênicas e as prováveis-patogênicas.

`combine` avalia `_pathogenic_rule(c) or _likely_pathogenic_rule(c)`, nessa ordem. A
consequência é que combinações que satisfazem uma regra PATH nunca chegam às regras LP —
e isso torna trechos de `_likely_pathogenic_rule` inalcançáveis pela API pública.

Escrito a partir da terceira rodada de mutação (CLAUDE.md §4.4), que deixou 14 vivos em
`_likely_pathogenic_rule`. Triados:

  12  cosmética do rótulo (`XX...XX`, caixa) — equivalentes para efeito de auditoria,
      justificados em tests/unit/acmg/test_rule_path_invariant.py
   1  `1 <= pm <= 2` -> `1 <= pm < 2` — sobrevivente FALSO: a suíte completa o mata em
      tests/test_property.py::test_more_pathogenic_evidence_never_more_benign, uma
      propriedade Hypothesis de monotonicidade. O arquivo entrou na seleção do mutmut.
   1  `1 <= pm <= 2` -> `1 <= pm <= 3` — genuinamente EQUIVALENTE, e é o que este arquivo
      documenta.

Ou seja: **zero lacunas reais** em `_likely_pathogenic_rule`. O resultado honesto desta
rodada é a ausência de defeito, mais uma correção na ferramenta de medição.

A equivalência depende de uma premissa verificável: `_likely_pathogenic_rule` tem UM único
chamador em todo o repositório — o `or` de `combine` — e portanto só é avaliada quando
`_pathogenic_rule` já devolveu None. Se alguém passar a chamá-la direto na produção, o
mutante `<= 3` deixa de ser equivalente (devolveria "LP-2" onde a original devolve "LP-4"),
e esta triagem precisa ser refeita.
"""

from __future__ import annotations

import pytest

from vcf2report.acmg import rules
from vcf2report.models import CriterionResult


def _met(code: str, strength: str) -> CriterionResult:
    return CriterionResult(
        code=code,
        name=f"{code} (fixture)",
        default_strength=strength,
        applies=True,
        met=True,
    )


@pytest.mark.parametrize(
    "n_moderados, tier_esperado, regra",
    [
        (1, rules.LIKELY_PATHOGENIC, "LP-2"),
        (2, rules.LIKELY_PATHOGENIC, "LP-2"),
        (3, rules.PATHOGENIC, "PATH-3"),
    ],
    ids=["1-moderado", "2-moderados", "3-moderados"],
)
def test_um_strong_mais_moderados_sobe_para_pathogenic_no_terceiro(
    n_moderados, tier_esperado, regra
):
    """A faixa de LP-2 é `1 <= pm <= 2`, e o limite superior existe porque em pm==3 quem
    responde é PATH-3, não LP-2.

    É esta precedência que torna o mutante `1 <= pm <= 2` -> `1 <= pm <= 3` **equivalente**:
    a condição alargada nunca é avaliada com pm==3, porque `combine` já retornou PATH-3.
    Um teste que forçasse `_likely_pathogenic_rule` diretamente mataria o mutante provando
    algo que a API pública não faz — o teste artificial que o §4.4 proíbe.

    O que este teste pina, em vez disso, é a razão: a fronteira entre LP-2 e PATH-3 em
    pm==3. Se alguém inverter a ordem de avaliação em `combine`, ou mexer no gatilho de
    PATH-3, ele falha — e é exatamente aí que o mutante deixaria de ser equivalente.
    """
    criterios = [_met("PS1", "strong")] + [
        _met(f"PM{i}", "moderate") for i in range(1, n_moderados + 1)
    ]

    tier, rule_path = rules.combine(criterios)

    assert tier == tier_esperado
    assert regra in rule_path, f"pm={n_moderados}: esperava {regra}, veio {rule_path}"


def test_a_regra_patogenica_e_avaliada_antes_da_provavel_patogenica():
    """A precedência em si, sem depender de qual combinação específica está em jogo.

    `combine` faz `_pathogenic_rule(c) or _likely_pathogenic_rule(c)`. PVS1 + 1 Strong +
    1 Moderate satisfaz as DUAS: PATH-1 por `pvs>=1 and ps>=1`, e LP-1 por
    `pvs>=1 and pm>=1`. As duas primeiras asserções não são cerimônia — se alguém mexer
    nos gatilhos e a fixture deixar de satisfazer ambas, o teste vira tautologia e passa
    sem testar nada. Elas falham antes disso.
    """
    contagens = {"PVS": 1, "PS": 1, "PM": 1, "PP": 0, "BA": 0, "BS": 0, "BP": 0}

    assert rules._pathogenic_rule(contagens) is not None, (
        "a fixture não satisfaz PATH-1"
    )
    assert rules._likely_pathogenic_rule(contagens) is not None, (
        "nem LP-1 — fixture inútil"
    )

    tier, rule_path = rules.combine(
        [_met("PVS1", "very_strong"), _met("PS1", "strong"), _met("PM1", "moderate")]
    )

    assert tier == rules.PATHOGENIC
    assert "PATH-1" in rule_path and "LP-1" not in rule_path
