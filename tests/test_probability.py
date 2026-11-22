"""Tests für INOD, SIGN, LIA und PROB mit den Beispielen von rosalind.info."""

import pytest

from rosalind import core


def test_inod():
    assert core.internal_nodes(4) == 2
    with pytest.raises(ValueError):
        core.internal_nodes(2)


def test_sign():
    perms = core.signed_permutations(2)
    assert len(perms) == 8  # 2! * 2^2
    assert (-1, 2) in perms and (2, -1) in perms
    assert len(set(perms)) == 8  # keine Duplikate


def test_lia_beispiel():
    assert core.independent_alleles(2, 1) == pytest.approx(0.684, abs=1e-3)


def test_prob_beispiel():
    result = core.random_string_log_probability(
        "ACGATACAA", [0.129, 0.287, 0.423, 0.476, 0.641, 0.742, 0.783]
    )
    expected = [-5.737, -5.217, -5.263, -5.360, -5.958, -6.628, -7.009]
    assert result == pytest.approx(expected, abs=1e-3)
