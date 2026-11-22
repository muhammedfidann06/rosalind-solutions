"""Tests für PMCH, MMCH, SSET und ASPC mit den Beispielen von rosalind.info."""

import pytest

from rosalind import core


def test_pmch():
    # Beispiel: AGCUAGUCAU -> 12  (#A = #U = 3 -> 3! = 6, #G = #C = 2 -> 2! = 2)
    assert core.perfect_matchings("AGCUAGUCAU") == 12
    with pytest.raises(ValueError):
        core.perfect_matchings("AAU")


def test_mmch():
    # Beispiel: AUGCUUC -> 6
    assert core.maximum_matchings("AUGCUUC") == 6


def test_sset():
    assert core.subset_count(3) == 8
    assert core.subset_count(30) == 2 ** 30 % 1_000_000


def test_aspc():
    # Beispiel: n = 6, m = 3 -> 42
    assert core.sum_of_combinations(6, 3) == 42
