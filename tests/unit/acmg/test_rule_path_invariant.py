"""O rule path tem de identificar a regra que realmente disparou.

`test_integrity.py::test_acmg_combining_matrix` verifica o **tier** de cada combinação e
nunca o rótulo. Consequência medida por mutação: 54 mutantes que mexem no texto do rule
path sobrevivem à suíte inteira — inclusive um (`pm >= 1` -> `pm > 1` na escolha do rótulo
de LP-1) que faz o laudo dizer "PVS1 + Supporting" quando quem disparou foi um Moderado.
Tier certo, trilha de auditoria errada.

A resposta NÃO é fixar as 54 strings. O rótulo é prosa e vai ser reescrito; congelá-lo
letra a letra criaria 54 testes que o CLAUDE.md §1.1 torna intocáveis para proteger
ortografia. O que precisa ser invariante é mais fraco e mais útil:

  1. cada combinação emite o identificador da SUA regra, e de nenhuma outra;
  2. onde o rótulo é escolhido em tempo de execução, a escolha corresponde à evidência;
  3. a trilha antes do `=>` nomeia só os critérios que dispararam.

Os três sobrevivem a qualquer reescrita do texto entre parênteses.
"""

from __future__ import annotations

import pytest

from vcf2report.acmg import rules
from vcf2report.models import CriterionResult

#: Todos os identificadores de regra que o combinador pode emitir. Usado para provar que
#: uma regra não emite o rótulo de outra — sem esta lista, um teste que só busca o id
#: esperado passaria com um rótulo que cita duas regras ao mesmo tempo.
TODOS_OS_IDS = (
    "PATH-1",
    "PATH-2",
    "PATH-3",
    "LP-1",
    "LP-2",
    "LP-3",
    "LP-4",
    "LP-5",
    "LP-6",
    "BEN-1",
    "BEN-2",
    "LB-1",
    "LB-2",
)


def _met(code: str, strength: str) -> CriterionResult:
    return CriterionResult(
        code=code,
        name=f"{code} (fixture)",
        default_strength=strength,
        applies=True,
        met=True,
    )


def _avaliado_e_nao_atendido(code: str, strength: str) -> CriterionResult:
    """O critério RODOU e a resposta foi não. `applies=True, met=False` é o estado que os
    avaliadores emitem para quase todo critério em quase toda variante."""
    return CriterionResult(
        code=code,
        name=f"{code} (fixture)",
        default_strength=strength,
        applies=True,
        met=False,
    )


def _pvs():
    return [_met("PVS1", "very_strong")]


def _ps(n: int):
    return [_met(f"PS{i}", "strong") for i in range(1, n + 1)]


def _pm(n: int):
    return [_met(f"PM{i}", "moderate") for i in range(1, n + 1)]


def _pp(n: int):
    return [_met(f"PP{i}", "supporting") for i in range(1, n + 1)]


def _bs(n: int):
    return [_met(f"BS{i}", "strong") for i in range(1, n + 1)]


def _bp(n: int):
    return [_met(f"BP{i}", "supporting") for i in range(1, n + 1)]


#: (descrição, critérios, tier esperado, identificador de regra esperado).
#: Uma linha por regra do combinador — se uma regra nova for adicionada sem entrar aqui,
#: test_toda_regra_do_combinador_esta_coberta falha.
CASOS = [
    ("PVS1 + 1 Strong", _pvs() + _ps(1), rules.PATHOGENIC, "PATH-1"),
    ("2 Strong", _ps(2), rules.PATHOGENIC, "PATH-2"),
    ("1 Strong + 3 Moderate", _ps(1) + _pm(3), rules.PATHOGENIC, "PATH-3"),
    ("PVS1 + 1 Moderate", _pvs() + _pm(1), rules.LIKELY_PATHOGENIC, "LP-1"),
    ("1 Strong + 1 Moderate", _ps(1) + _pm(1), rules.LIKELY_PATHOGENIC, "LP-2"),
    ("1 Strong + 2 Supporting", _ps(1) + _pp(2), rules.LIKELY_PATHOGENIC, "LP-3"),
    ("3 Moderate", _pm(3), rules.LIKELY_PATHOGENIC, "LP-4"),
    ("2 Moderate + 2 Supporting", _pm(2) + _pp(2), rules.LIKELY_PATHOGENIC, "LP-5"),
    ("1 Moderate + 4 Supporting", _pm(1) + _pp(4), rules.LIKELY_PATHOGENIC, "LP-6"),
    ("BA1", [_met("BA1", "stand_alone")], rules.BENIGN, "BEN-1"),
    ("2 Strong benigno", _bs(2), rules.BENIGN, "BEN-2"),
    ("1 Strong + 1 Supporting benigno", _bs(1) + _bp(1), rules.LIKELY_BENIGN, "LB-1"),
    ("2 Supporting benigno", _bp(2), rules.LIKELY_BENIGN, "LB-2"),
]


@pytest.mark.parametrize(
    "descricao, criterios, tier_esperado, id_esperado",
    CASOS,
    ids=[c[3] for c in CASOS],
)
def test_a_regra_que_disparou_se_identifica_no_rule_path(
    descricao, criterios, tier_esperado, id_esperado
):
    """Invariante 1: o rótulo nomeia a regra que decidiu, e só ela.

    Isto é deliberadamente indiferente ao texto depois do identificador — reescrever
    "(>=3 Moderate)" para "(3 ou mais Moderados)" não quebra nada aqui. O que não pode
    mudar em silêncio é QUAL regra o laudo diz ter aplicado.
    """
    tier, rule_path = rules.combine(criterios)

    assert tier == tier_esperado, f"{descricao}: {tier} (esperado {tier_esperado})"
    assert id_esperado in rule_path, (
        f"{descricao}: rule path sem {id_esperado} -> {rule_path}"
    )

    outros = [i for i in TODOS_OS_IDS if i != id_esperado and i in rule_path]
    assert not outros, (
        f"{descricao}: rule path cita outra(s) regra(s) {outros} -> {rule_path}"
    )


def test_toda_regra_do_combinador_esta_coberta():
    """Uma regra nova sem caso aqui passaria despercebida, e o teste acima só verifica as
    que já existem. Ancorado na lista de identificadores, que é a mesma que o invariante 1
    usa para detectar contaminação entre rótulos."""
    cobertos = {caso[3] for caso in CASOS}

    assert cobertos == set(TODOS_OS_IDS), (
        f"regras sem caso de teste: {set(TODOS_OS_IDS) - cobertos}; "
        f"casos para regras inexistentes: {cobertos - set(TODOS_OS_IDS)}"
    )


def test_a_trilha_lista_so_os_criterios_que_dispararam():
    """Invariante 3: a trilha é a lista de EVIDÊNCIAS, não a lista de critérios avaliados.

    `combine` monta `met = [cr.code for cr in criteria if cr.applies and cr.met]`. Trocar
    esse `and` por `or` não move o tier nem o rótulo da regra — os dois vêm de `_counts`,
    que filtra de novo — e faz o laudo imprimir
    `PVS1 + PM1 + BA1 + BP4 => Likely Pathogenic [LP-1 ...]`. BA1 e BP4 foram avaliados e
    NÃO dispararam; apareceriam como se tivessem pesado, e ainda por cima como evidência
    benigna dentro de uma conclusão patogênica. É a invariante de honestidade do projeto
    ("ausência de dado nunca é evidência") aplicada à linha que o laudo cita como trilha de
    auditoria. Sobreviveu à suíte inteira até esta rodada de mutação.
    """
    criterios = (
        _pvs()
        + _pm(1)
        + [
            _avaliado_e_nao_atendido("BA1", "stand_alone"),
            _avaliado_e_nao_atendido("BP4", "supporting"),
        ]
    )

    tier, rule_path = rules.combine(criterios)
    trilha = rule_path.split("=>")[0]

    assert tier == rules.LIKELY_PATHOGENIC, f"a fixture mudou de tier: {rule_path}"
    assert "PVS1" in trilha and "PM1" in trilha, (
        f"trilha perdeu evidência real: {trilha}"
    )
    for code in ("BA1", "BP4"):
        assert code not in trilha, (
            f"{code} foi avaliado e não atendido, e está na trilha: {rule_path}"
        )


@pytest.mark.parametrize(
    "extra, classe_esperada, classe_proibida",
    [
        (_pm(1), "Moderate", "Supporting"),
        (_pp(1), "Supporting", "Moderate"),
    ],
    ids=["moderado", "supporting"],
)
def test_lp1_nomeia_a_classe_de_criterio_que_realmente_disparou(
    extra, classe_esperada, classe_proibida
):
    """Invariante 2: onde o rótulo é escolhido em tempo de execução, ele descreve a
    evidência real.

    LP-1 é a única regra com rótulo dinâmico — "PVS1 + Moderate" ou "PVS1 + Supporting",
    conforme o que acompanhou o PVS1. O mutante `pm >= 1` -> `pm > 1` nessa escolha mantém
    o tier e troca a explicação, que é o tipo de erro que uma seção chamada "auditável"
    não pode ter.
    """
    tier, rule_path = rules.combine(_pvs() + extra)

    assert tier == rules.LIKELY_PATHOGENIC
    assert "LP-1" in rule_path
    assert classe_esperada in rule_path
    assert classe_proibida not in rule_path
