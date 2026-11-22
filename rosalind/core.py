"""
Rosalind problem solutions (bioinformatics).
"""

from collections import Counter
from itertools import product
from typing import Dict, Iterable, List, Tuple

# Das ist die Codon-Tabelle. Drei Buchstaben (RNA) geben eine Aminosäure.
# "*" bedeutet: Hier ist das Protein zu Ende (Stopp).
_BASES = "UCAG"
_AMINO = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODON_TABLE: Dict[str, str] = {
    a + b + c: _AMINO[i]
    for i, (a, b, c) in enumerate(product(_BASES, repeat=3))
}

# Masse von jeder Aminosäure (monoisotopisch), für das Problem PRTM.
MONO_MASS: Dict[str, float] = {
    "A": 71.03711, "C": 103.00919, "D": 115.02694, "E": 129.04259,
    "F": 147.06841, "G": 57.02146, "H": 137.05891, "I": 113.08406,
    "K": 128.09496, "L": 113.08406, "M": 131.04049, "N": 114.04293,
    "P": 97.05276, "Q": 128.05858, "R": 156.10111, "S": 87.03203,
    "T": 101.04768, "V": 99.06841, "W": 186.07931, "Y": 163.06333,
}


# ---------------------------------------------------------------------------
# DNA – Nukleotide zählen
# ---------------------------------------------------------------------------
def count_nucleotides(dna: str) -> Tuple[int, int, int, int]:
    # Wir zählen jeden Buchstaben: A, C, G und T.
    c = Counter(dna.upper())
    return c["A"], c["C"], c["G"], c["T"]


# ---------------------------------------------------------------------------
# RNA – DNA in RNA umschreiben
# ---------------------------------------------------------------------------
def transcribe(dna: str) -> str:
    # In RNA gibt es kein T. Wir schreiben U statt T.
    return dna.upper().replace("T", "U")


# ---------------------------------------------------------------------------
# REVC – Komplementärer Strang
# ---------------------------------------------------------------------------
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(dna: str) -> str:
    # Zuerst tauschen wir die Basen (A<->T, C<->G).
    # Dann lesen wir den Text von hinten nach vorne.
    return dna.upper().translate(_COMPLEMENT)[::-1]


# ---------------------------------------------------------------------------
# FIB / FIBD – Kaninchen
# ---------------------------------------------------------------------------
def rabbits(n: int, k: int) -> int:
    # Jeden Monat bekommt ein erwachsenes Paar k neue Paare.
    # Wir merken uns nur die letzten zwei Monate.
    young, old = 1, 0
    for _ in range(n - 1):
        young, old = old * k, old + young
    return young + old


def mortal_rabbits(n: int, m: int) -> int:
    # Hier sterben die Kaninchen nach m Monaten.
    # ages[i] = Anzahl Paare mit Alter i.
    ages = [0] * m
    ages[0] = 1
    for _ in range(n - 1):
        babies = sum(ages[1:])          # Nur Erwachsene bekommen Babys.
        ages = [babies] + ages[:-1]     # Alle werden einen Monat älter.
    return sum(ages)


# ---------------------------------------------------------------------------
# GC – GC-Gehalt
# ---------------------------------------------------------------------------
def parse_fasta(text: str) -> Dict[str, str]:
    # Eine FASTA-Datei hat Zeilen mit ">" (Name) und dann die Sequenz.
    records: Dict[str, str] = {}
    name = None
    for line in text.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = ""
        elif name is not None:
            records[name] += line.upper()
    return records


def gc_content(dna: str) -> float:
    # Wie viel Prozent sind G oder C?
    if not dna:
        return 0.0
    dna = dna.upper()
    return 100.0 * (dna.count("G") + dna.count("C")) / len(dna)


def highest_gc(fasta: str) -> Tuple[str, float]:
    # Wir suchen die Sequenz mit dem höchsten GC-Wert.
    records = parse_fasta(fasta)
    best = max(records, key=lambda k: gc_content(records[k]))
    return best, gc_content(records[best])


# ---------------------------------------------------------------------------
# HAMM – Hamming-Distanz
# ---------------------------------------------------------------------------
def hamming(s: str, t: str) -> int:
    # Die zwei Texte müssen gleich lang sein.
    if len(s) != len(t):
        raise ValueError("Die Sequenzen sind nicht gleich lang.")
    # Wir zählen die Stellen, wo die Buchstaben anders sind.
    return sum(a != b for a, b in zip(s, t))


# ---------------------------------------------------------------------------
# IPRB – Mendel
# ---------------------------------------------------------------------------
def dominant_probability(k: int, m: int, n: int) -> float:
    # k = homozygot dominant, m = heterozygot, n = homozygot rezessiv.
    # Wir rechnen die Chance für ein rezessives Kind und nehmen 1 minus das.
    total = k + m + n
    pairs = total * (total - 1)
    recessive = (
        n * (n - 1) * 1.0          # nn x nn -> immer rezessiv
        + 2 * n * m * 0.5          # nn x Mm -> halb rezessiv
        + m * (m - 1) * 0.25       # Mm x Mm -> ein Viertel rezessiv
    )
    return 1 - recessive / pairs


# ---------------------------------------------------------------------------
# PROT – RNA in Protein übersetzen
# ---------------------------------------------------------------------------
def translate(rna: str) -> str:
    protein: List[str] = []
    # Wir lesen immer drei Buchstaben zusammen (ein Codon).
    for i in range(0, len(rna) - 2, 3):
        aa = CODON_TABLE[rna[i:i + 3].upper()]
        if aa == "*":
            break  # Stopp-Codon: Wir hören auf.
        protein.append(aa)
    return "".join(protein)


# ---------------------------------------------------------------------------
# SUBS – Motiv finden
# ---------------------------------------------------------------------------
def find_motif(s: str, t: str) -> List[int]:
    # Wir suchen alle Positionen von t in s.
    # Rosalind zählt ab 1, nicht ab 0. Darum "+ 1".
    return [i + 1 for i in range(len(s) - len(t) + 1) if s[i:i + len(t)] == t]


# ---------------------------------------------------------------------------
# CONS – Konsens-Sequenz
# ---------------------------------------------------------------------------
def consensus(seqs: Iterable[str]) -> Tuple[str, Dict[str, List[int]]]:
    seqs = list(seqs)
    length = len(seqs[0])
    # Für jede Position zählen wir A, C, G, T.
    profile = {b: [0] * length for b in "ACGT"}
    for s in seqs:
        for i, base in enumerate(s):
            profile[base][i] += 1
    # Der häufigste Buchstabe an jeder Stelle gewinnt.
    cons = "".join(
        max("ACGT", key=lambda b: profile[b][i]) for i in range(length)
    )
    return cons, profile


# ---------------------------------------------------------------------------
# GRPH – Overlap-Graph
# ---------------------------------------------------------------------------
def overlap_graph(fasta: str, k: int = 3) -> List[Tuple[str, str]]:
    records = parse_fasta(fasta)
    edges = []
    # Eine Kante gibt es, wenn das Ende von s gleich dem Anfang von t ist.
    for a, s in records.items():
        for b, t in records.items():
            if a != b and s[-k:] == t[:k]:
                edges.append((a, b))
    return edges


# ---------------------------------------------------------------------------
# LCSM – Längstes gemeinsames Motiv
# ---------------------------------------------------------------------------
def longest_shared_motif(seqs: Iterable[str]) -> str:
    seqs = sorted(seqs, key=len)
    shortest, others = seqs[0], seqs[1:]
    # Wir probieren zuerst lange Stücke, dann kürzere.
    for size in range(len(shortest), 0, -1):
        for start in range(len(shortest) - size + 1):
            cand = shortest[start:start + size]
            if all(cand in o for o in others):
                return cand
    return ""


# ---------------------------------------------------------------------------
# PRTM – Protein-Masse
# ---------------------------------------------------------------------------
def protein_mass(protein: str) -> float:
    # Wir addieren die Masse von jeder Aminosäure.
    return round(sum(MONO_MASS[aa] for aa in protein.upper()), 3)


# ---------------------------------------------------------------------------
# MRNA – Wie viele RNA-Sequenzen ergeben dieses Protein?
# ---------------------------------------------------------------------------
def rna_count_for_protein(protein: str, modulo: int = 1_000_000) -> int:
    # Für jede Aminosäure zählen wir, wie viele Codons sie haben.
    codons_per_aa = Counter(CODON_TABLE.values())
    count = codons_per_aa["*"]  # am Ende kommt ein Stopp-Codon (3 Möglichkeiten)
    for aa in protein.upper():
        # Wir rechnen sofort modulo, damit die Zahl nicht riesig wird.
        count = (count * codons_per_aa[aa]) % modulo
    return count


# ---------------------------------------------------------------------------
# SPLC – RNA-Spleißen: Introns entfernen, dann übersetzen
# ---------------------------------------------------------------------------
def splice_and_translate(fasta: str) -> str:
    # Die erste Sequenz ist die DNA, alle anderen sind Introns.
    seqs = list(parse_fasta(fasta).values())
    dna, introns = seqs[0], seqs[1:]
    for intron in introns:
        dna = dna.replace(intron, "")
    return translate(transcribe(dna))


# ---------------------------------------------------------------------------
# ORF – Offene Leserahmen in allen 6 Rahmen
# ---------------------------------------------------------------------------
def open_reading_frames(dna: str) -> List[str]:
    # Ein ORF beginnt mit AUG (Start) und endet mit einem Stopp-Codon.
    # Wir lesen beide Stränge, jeweils in 3 Rahmen = 6 Rahmen.
    proteins = set()
    for strand in (dna.upper(), reverse_complement(dna)):
        rna = transcribe(strand)
        for start in range(len(rna) - 2):
            if rna[start:start + 3] != "AUG":
                continue
            protein: List[str] = []
            for i in range(start, len(rna) - 2, 3):
                aa = CODON_TABLE[rna[i:i + 3]]
                if aa == "*":
                    proteins.add("".join(protein))  # nur mit Stopp-Codon gültig
                    break
                protein.append(aa)
    return sorted(proteins)


# ---------------------------------------------------------------------------
# PERM – Alle Permutationen von 1..n
# ---------------------------------------------------------------------------
def permutations_of(n: int) -> List[Tuple[int, ...]]:
    # itertools.permutations liefert schon die lexikografische Reihenfolge.
    from itertools import permutations
    return list(permutations(range(1, n + 1)))


# ---------------------------------------------------------------------------
# LEXF – Alle Wörter der Länge n aus einem Alphabet, sortiert
# ---------------------------------------------------------------------------
def lexicographic_kmers(alphabet: List[str], n: int) -> List[str]:
    # product() geht die Buchstaben in der gegebenen Reihenfolge durch.
    # Das Alphabet bestimmt also selbst, was "alphabetisch" heißt.
    return ["".join(p) for p in product(alphabet, repeat=n)]


# ---------------------------------------------------------------------------
# SSEQ – Gespleißtes Motiv: t als Teilfolge (nicht zusammenhängend) in s
# ---------------------------------------------------------------------------
def spliced_motif(s: str, t: str) -> List[int]:
    # Gierig von links: Für jeden Buchstaben von t nehmen wir das nächste passende in s.
    positions: List[int] = []
    i = 0
    for ch in t:
        while i < len(s) and s[i] != ch:
            i += 1
        if i == len(s):
            return []  # t ist keine Teilfolge von s
        positions.append(i + 1)  # Rosalind zählt ab 1
        i += 1
    return positions


# ---------------------------------------------------------------------------
# IEV – Erwartete Anzahl Nachkommen mit dominantem Merkmal
# ---------------------------------------------------------------------------
def expected_dominant_offspring(couples: List[int], children_per_couple: int = 2) -> float:
    # Paare in dieser Reihenfolge: AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa
    # Wahrscheinlichkeit für ein dominantes Kind bei jedem Paar:
    p_dominant = [1.0, 1.0, 1.0, 0.75, 0.5, 0.0]
    # Erwartungswert = Summe (Anzahl Paare * Kinder * Wahrscheinlichkeit)
    return sum(n * children_per_couple * p for n, p in zip(couples, p_dominant))


# ---------------------------------------------------------------------------
# TRAN – Verhältnis Transition / Transversion
# ---------------------------------------------------------------------------
_PURINES = set("AG")      # Purine: A, G
_PYRIMIDINES = set("CT")  # Pyrimidine: C, T


def transition_transversion_ratio(s: str, t: str) -> float:
    # Transition: Purin <-> Purin oder Pyrimidin <-> Pyrimidin (A<->G, C<->T)
    # Transversion: Purin <-> Pyrimidin
    if len(s) != len(t):
        raise ValueError("Die Sequenzen sind nicht gleich lang.")
    transitions = transversions = 0
    for a, b in zip(s.upper(), t.upper()):
        if a == b:
            continue
        same_group = ({a, b} <= _PURINES) or ({a, b} <= _PYRIMIDINES)
        if same_group:
            transitions += 1
        else:
            transversions += 1
    if transversions == 0:
        raise ValueError("Keine Transversion – das Verhältnis ist nicht definiert.")
    return transitions / transversions


# ---------------------------------------------------------------------------
# KMER – k-mer-Zusammensetzung (wie oft kommt jedes k-mer vor?)
# ---------------------------------------------------------------------------
def kmer_composition(dna: str, k: int = 4) -> List[int]:
    # Alle k-mere über ACGT, lexikografisch sortiert (wie LEXF).
    kmers = lexicographic_kmers(["A", "C", "G", "T"], k)
    counts = Counter(dna[i:i + k] for i in range(len(dna) - k + 1))
    return [counts[km] for km in kmers]


# ---------------------------------------------------------------------------
# PDST – Matrix der p-Distanzen
# ---------------------------------------------------------------------------
def p_distance_matrix(seqs: List[str]) -> List[List[float]]:
    # p-Distanz = Anteil der Stellen, an denen sich zwei Sequenzen unterscheiden.
    n = len(seqs)
    return [[hamming(seqs[i], seqs[j]) / len(seqs[i]) for j in range(n)] for i in range(n)]


# ---------------------------------------------------------------------------
# INOD – Innere Knoten eines ungewurzelten Binärbaums
# ---------------------------------------------------------------------------
def internal_nodes(leaves: int) -> int:
    # Jeder innere Knoten hat Grad 3. Aus "Kanten = Knoten - 1" folgt: innere = Blätter - 2.
    if leaves < 3:
        raise ValueError("Ein ungewurzelter Binärbaum braucht mindestens 3 Blätter.")
    return leaves - 2


# ---------------------------------------------------------------------------
# SIGN – Vorzeichenbehaftete Permutationen
# ---------------------------------------------------------------------------
def signed_permutations(n: int) -> List[Tuple[int, ...]]:
    # Jede Permutation, und jede Zahl darf + oder - sein -> n! * 2^n Möglichkeiten.
    from itertools import permutations
    out: List[Tuple[int, ...]] = []
    for perm in permutations(range(1, n + 1)):
        for signs in product((1, -1), repeat=n):
            out.append(tuple(s * p for s, p in zip(signs, perm)))
    return out


# ---------------------------------------------------------------------------
# LIA – Mendels zweites Gesetz (unabhängige Allele)
# ---------------------------------------------------------------------------
def independent_alleles(k: int, n: int) -> float:
    # Mit einem AaBb-Partner ist jedes Kind mit Wahrscheinlichkeit 1/4 AaBb – egal welche Eltern.
    # In Generation k gibt es 2^k Organismen. Gesucht: P(mindestens n davon sind AaBb).
    from math import comb
    total = 2 ** k
    p = 0.25
    return 1 - sum(comb(total, i) * p ** i * (1 - p) ** (total - i) for i in range(n))


# ---------------------------------------------------------------------------
# PROB – Wahrscheinlichkeit einer Sequenz bei gegebenem GC-Gehalt (log10)
# ---------------------------------------------------------------------------
def random_string_log_probability(dna: str, gc_contents: List[float]) -> List[float]:
    # Bei GC-Gehalt x: P(G) = P(C) = x/2, P(A) = P(T) = (1-x)/2.
    # Wir addieren Logarithmen statt Wahrscheinlichkeiten zu multiplizieren (kein Unterlauf).
    from math import log10
    gc = sum(1 for ch in dna.upper() if ch in "GC")
    at = len(dna) - gc
    return [gc * log10(x / 2) + at * log10((1 - x) / 2) for x in gc_contents]


# ---------------------------------------------------------------------------
# REVP – Restriktionsstellen (Palindrome)
# ---------------------------------------------------------------------------
def reverse_palindromes(dna: str, lo: int = 4, hi: int = 12) -> List[Tuple[int, int]]:
    out = []
    # Ein "reverses Palindrom" ist gleich wie sein Komplement von hinten.
    for i in range(len(dna)):
        for size in range(lo, hi + 1):
            piece = dna[i:i + size]
            if len(piece) < size:
                break
            if piece == reverse_complement(piece):
                out.append((i + 1, size))
    return out
