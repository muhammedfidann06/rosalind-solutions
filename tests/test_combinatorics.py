"""Tests für PERM, LEXF und SSEQ mit den Beispielen von rosalind.info."""

from rosalind import core


def test_perm():
    perms = core.permutations_of(3)
    assert len(perms) == 6
    assert perms[0] == (1, 2, 3)
    assert perms[-1] == (3, 2, 1)


def test_lexf():
    words = core.lexicographic_kmers(["A", "C", "G", "T"], 2)
    assert len(words) == 16
    assert words[:5] == ["AA", "AC", "AG", "AT", "CA"]
    assert words[-1] == "TT"


def test_sseq():
    pos = core.spliced_motif("ACGTACGTGACG", "GTA")
    # Jede Lösung ist gültig; wir prüfen, dass die Positionen wirklich "GTA" ergeben.
    s = "ACGTACGTGACG"
    assert "".join(s[p - 1] for p in pos) == "GTA"
    assert pos == sorted(pos)
    assert core.spliced_motif("AAA", "C") == []
