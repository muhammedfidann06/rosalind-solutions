# Tests mit den Beispielen von rosalind.info.
# Starten mit: python -m pytest
import pytest

from rosalind import core


def test_dna():
    s = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
    assert core.count_nucleotides(s) == (20, 12, 17, 21)


def test_rna():
    assert core.transcribe("GATGGAACTTGACTACGTAAATT") == "GAUGGAACUUGACUACGUAAAUU"


def test_revc():
    assert core.reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"


def test_fib():
    assert core.rabbits(5, 3) == 19


def test_fibd():
    assert core.mortal_rabbits(6, 3) == 4


def test_gc():
    fasta = """>Rosalind_6404
CCTGCGGAAGATCGGCACTAGAATAGCCAGAACCGTTTCTCTGAGGCTTCCGGCCTTCCC
TCCCACTAATAATTCTGAGG
>Rosalind_5959
CCATCGGTAGCGCATCCTTAGTCCAATTAAGTCCCTATCCAGGCGCTCCGCCGAAGGTCT
ATATCCATTTGTCAGCAGACACGC
>Rosalind_0808
CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGAC
TGGGAACCTGCGGGCAGTAGGTGGAAT"""
    name, value = core.highest_gc(fasta)
    assert name == "Rosalind_0808"
    assert value == pytest.approx(60.919540, abs=1e-3)


def test_hamm():
    assert core.hamming("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7


def test_iprb():
    assert core.dominant_probability(2, 2, 2) == pytest.approx(0.78333, abs=1e-4)


def test_prot():
    rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    assert core.translate(rna) == "MAMAPRTEINSTRING"


def test_subs():
    assert core.find_motif("GATATATGCATATACTT", "ATAT") == [2, 4, 10]


def test_cons():
    seqs = ["ATCCAGCT", "GGGCAACT", "ATGGATCT", "AAGCAACC",
            "TTGGAACT", "ATGCCATT", "ATGGCACT"]
    cons, profile = core.consensus(seqs)
    assert cons == "ATGCAACT"
    assert profile["A"] == [5, 1, 0, 0, 5, 5, 0, 0]


def test_grph():
    fasta = """>Rosalind_0498
AAATAAA
>Rosalind_2391
AAATTTT
>Rosalind_2323
TTTTCCC
>Rosalind_0442
AAATCCC
>Rosalind_5013
GGGTGGG"""
    edges = set(core.overlap_graph(fasta))
    assert edges == {
        ("Rosalind_0498", "Rosalind_2391"),
        ("Rosalind_0498", "Rosalind_0442"),
        ("Rosalind_2391", "Rosalind_2323"),
    }


def test_lcsm():
    result = core.longest_shared_motif(["GATTACA", "TAGACCA", "ATACA"])
    # Es gibt mehrere richtige Antworten mit Länge 2.
    assert len(result) == 2
    assert all(result in s for s in ["GATTACA", "TAGACCA", "ATACA"])


def test_prtm():
    assert core.protein_mass("SKADYEK") == pytest.approx(821.392, abs=1e-3)


def test_revp():
    result = core.reverse_palindromes("TCAATGCATGCGGGTCTATATGCAT")
    assert (4, 6) in result and (17, 4) in result and len(result) == 8
