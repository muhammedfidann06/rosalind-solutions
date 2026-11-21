"""Tests für MRNA, SPLC und ORF mit den Beispielen von rosalind.info."""

from rosalind import core


def test_mrna():
    # "MA": M hat 1 Codon, A hat 4, Stopp hat 3 -> 1 * 4 * 3 = 12
    assert core.rna_count_for_protein("MA") == 12


def test_splc():
    fasta = """>Rosalind_10
ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG
>Rosalind_12
ATCGGTCGAA
>Rosalind_15
ATCGGTCGAGCGTGT
"""
    assert core.splice_and_translate(fasta) == "MVYIADKQHVASREAYGHMFKVCA"


def test_orf():
    dna = (
        "AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAG"
        "CCTGAATGATCCGAGTAGCATCTCAG"
    )
    assert core.open_reading_frames(dna) == sorted(
        ["MLLGSFRLIPKETLIQVAGSSPCNLS", "M", "MGMTPRLGLESLLE", "MTPRLGLESLLE"]
    )
