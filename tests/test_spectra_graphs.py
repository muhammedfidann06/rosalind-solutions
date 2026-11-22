"""Tests für SPEC, SCSP, DBRU und CONV mit den Beispielen von rosalind.info."""

import pytest

from rosalind import core


def _is_subsequence(sub, seq):
    it = iter(seq)
    return all(x in it for x in sub)


def test_spec():
    masses = [3524.8542, 3710.9335, 3841.974, 3970.0326, 4057.0646]
    assert core.protein_from_spectrum(masses) == "WMQS"


def test_scsp():
    s, t = "ATCTGAT", "TGCATA"
    sup = core.shortest_common_supersequence(s, t)
    assert len(sup) == 9  # z. B. "ATGCATGAT"
    assert _is_subsequence(s, sup) and _is_subsequence(t, sup)


def test_dbru():
    reads = ["TGAT", "CATG", "TCAT", "ATGC", "CATC", "CATC"]
    expected = {
        ("ATC", "TCA"), ("ATG", "TGA"), ("ATG", "TGC"), ("CAT", "ATC"), ("CAT", "ATG"),
        ("GAT", "ATG"), ("GCA", "CAT"), ("TCA", "CAT"), ("TGA", "GAT"),
    }
    assert set(core.de_bruijn_edges(reads)) == expected


def test_conv():
    s1 = [186.07931, 287.12699, 548.20532, 580.18077, 681.22845, 706.27446, 782.27613, 968.35544, 968.35544]
    s2 = [101.04768, 158.06914, 202.09536, 318.09979, 419.14747, 463.17369]
    count, value = core.spectral_convolution(s1, s2)
    assert count == 3
    assert value == pytest.approx(85.03163, abs=1e-4)
