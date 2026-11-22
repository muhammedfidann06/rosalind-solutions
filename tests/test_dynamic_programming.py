"""Tests für LGIS, LCSQ und EDIT."""

from rosalind import core


def _is_subsequence(sub, seq):
    it = iter(seq)
    return all(x in it for x in sub)


def test_lgis_beispiel_von_rosalind():
    seq = [5, 1, 4, 2, 3]
    inc = core.longest_increasing_subsequence(seq)
    dec = core.longest_decreasing_subsequence(seq)
    assert len(inc) == 3 and inc == sorted(inc) and _is_subsequence(inc, seq)
    assert len(dec) == 3 and dec == sorted(dec, reverse=True) and _is_subsequence(dec, seq)


def test_lgis_leer():
    assert core.longest_increasing_subsequence([]) == []


def test_lcsq_beispiel():
    s, t = "AACCTTGG", "ACACTGTGA"
    lcs = core.longest_common_subsequence(s, t)
    assert len(lcs) == 6  # z. B. "AACTGG"
    assert _is_subsequence(lcs, s) and _is_subsequence(lcs, t)


def test_edit_beispiel_von_rosalind():
    assert core.edit_distance("PLEASANTLY", "MEANLY") == 5
    assert core.edit_distance("", "ABC") == 3
    assert core.edit_distance("ABC", "ABC") == 0
