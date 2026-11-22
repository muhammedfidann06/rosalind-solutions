"""Tests für TREE, CORR, LONG und LEXV mit den Beispielen von rosalind.info."""

import pytest

from rosalind import core


def test_tree():
    edges = [(1, 2), (2, 8), (4, 10), (5, 9), (6, 10), (7, 9)]
    assert core.edges_to_complete_tree(10, edges) == 3


def test_corr_beispiel():
    reads = [
        "TCATC", "TTCAT", "TCATC", "TGAAA", "GAGGA", "TTTCA",
        "ATCAA", "TTGAT", "TTTCC",
    ]
    fixes = set(core.correct_reads(reads))
    assert fixes == {
        ("TTCAT", "TTGAT"),
        ("GAGGA", "GATGA"),
        ("TTTCC", "TTTCA"),
    }


def test_long_beispiel():
    reads = ["ATTAGACCTG", "CCTGCCGGAA", "AGACCTGCCG", "GCCGGAATAC"]
    assert core.shortest_superstring(reads) == "ATTAGACCTGCCGGAATAC"


def test_long_fehler():
    with pytest.raises(ValueError):
        core.shortest_superstring(["AAAA", "CCCC"])


def test_lexv_beispiel():
    words = core.lexicographic_variable_length(["D", "N", "A"], 3)
    assert words[:5] == ["D", "DD", "DDD", "DDN", "DDA"]
    assert words[-1] == "AAA"
    assert len(words) == 3 + 9 + 27
