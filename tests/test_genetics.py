"""Tests für IEV und TRAN."""

import pytest

from rosalind import core


def test_iev_beispiel_von_rosalind():
    # Beispiel von rosalind.info: 1 0 0 1 0 1 -> 3.5
    assert core.expected_dominant_offspring([1, 0, 0, 1, 0, 1]) == pytest.approx(3.5)


def test_tran_kleines_beispiel():
    # A->G und C->T sind Transitionen, G->C und T->A sind Transversionen -> 2 / 2
    assert core.transition_transversion_ratio("ACGT", "GTCA") == pytest.approx(1.0)


def test_tran_fehler():
    with pytest.raises(ValueError):
        core.transition_transversion_ratio("AC", "A")
    with pytest.raises(ValueError):
        core.transition_transversion_ratio("AG", "GA")  # nur Transitionen
