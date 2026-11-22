"""Tests für SEXL, AFRQ, EVAL und KMP mit den Beispielen von rosalind.info."""

import pytest

from rosalind import core


def test_sexl():
    assert core.sex_linked_carriers([0.1, 0.5, 0.8]) == pytest.approx([0.18, 0.5, 0.32])


def test_afrq():
    assert core.carrier_probability([0.1, 0.25, 0.5]) == pytest.approx([0.532, 0.75, 0.914], abs=1e-3)


def test_eval():
    result = core.expected_motif_count(10, "AG", [0.25, 0.5, 0.75])
    assert result == pytest.approx([0.422, 0.563, 0.422], abs=1e-3)


def test_kmp_failure_array():
    s = "CAGCATGGTATCACAGCAGAG"
    expected = [0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 3, 4, 5, 3, 0, 0]
    assert core.failure_array(s) == expected
    assert core.failure_array("") == []
