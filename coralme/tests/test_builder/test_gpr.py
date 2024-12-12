import pytest
from coralme.builder.gpr import expand_gpr

@pytest.mark.parametrize("test_input,expected",
                         [
                             ("A or B",[['A'], ['B']]),
                             ("A or (B and C)", [['A'], ['B', 'C']]),
                             ("A and (B or C)", [['A', 'B'], ['A', 'C']]),
                             ("(B and C) or A", [['B', 'C'], ['A']]),
                             ("(A or B) and (B or C)", [['A', 'B'], ['A', 'C'], ['B', 'B'], ['B', 'C']]),
                             ("(A and B) and (B or C)", [['A', 'B', 'B'], ['A', 'B', 'C']]),
                             ("(A or B) and (B and C)", [['A', 'B', 'C'], ['B', 'B', 'C']]),
                             ("(A and (B and (C or D))) and (E and F)", [['A', 'B', 'C', 'E', 'F'], ['A', 'B', 'D', 'E', 'F']])
                             ])
def test_expand_gpr(test_input, expected):
    assert expand_gpr(test_input) == expected
