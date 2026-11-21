"""Tests für KMER und PDST."""

import pytest

from rosalind import core


def test_kmer_zaehlt_ueberlappend():
    counts = core.kmer_composition("AAAA", k=2)
    assert len(counts) == 16          # 4^2 mögliche 2-mere
    assert counts[0] == 3             # "AA" kommt 3-mal vor (überlappend)
    assert sum(counts) == 3


def test_kmer_reihenfolge():
    counts = core.kmer_composition("ACGT", k=2)
    # AC ist Index 1, CG ist Index 6, GT ist Index 11 (lexikografisch über ACGT)
    assert [i for i, c in enumerate(counts) if c] == [1, 6, 11]


def test_pdst_beispiel_von_rosalind():
    seqs = ["TTTCCATTTA", "GATTCATTTC", "TTTCCATTTT", "GTTCCATTTA"]
    m = core.p_distance_matrix(seqs)
    expected = [
        [0.0, 0.4, 0.1, 0.1],
        [0.4, 0.0, 0.4, 0.3],
        [0.1, 0.4, 0.0, 0.2],
        [0.1, 0.3, 0.2, 0.0],
    ]
    for row, exp in zip(m, expected):
        assert row == pytest.approx(exp)
