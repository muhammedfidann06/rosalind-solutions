"""
Kleines Kommandozeilen-Tool.

Beispiel:
    python solve.py revc data/rosalind_revc.txt
"""
import argparse
from pathlib import Path

from rosalind import core


def main() -> None:
    # Hier lesen wir die Argumente: Problem-ID und Datei.
    p = argparse.ArgumentParser(description="Rosalind solver")
    p.add_argument("problem", help="z.B. dna, rna, revc, gc, hamm, prot ...")
    p.add_argument("file", type=Path)
    args = p.parse_args()

    text = args.file.read_text().strip()
    lines = text.splitlines()
    pid = args.problem.lower()

    # Für jedes Problem rufen wir die passende Funktion auf.
    if pid == "dna":
        print(*core.count_nucleotides(text))
    elif pid == "rna":
        print(core.transcribe(text))
    elif pid == "revc":
        print(core.reverse_complement(text))
    elif pid == "fib":
        n, k = map(int, text.split())
        print(core.rabbits(n, k))
    elif pid == "fibd":
        n, m = map(int, text.split())
        print(core.mortal_rabbits(n, m))
    elif pid == "gc":
        name, value = core.highest_gc(text)
        print(name)
        print(f"{value:.6f}")
    elif pid == "hamm":
        print(core.hamming(lines[0], lines[1]))
    elif pid == "iprb":
        k, m, n = map(int, text.split())
        print(f"{core.dominant_probability(k, m, n):.5f}")
    elif pid == "prot":
        print(core.translate(text))
    elif pid == "subs":
        print(*core.find_motif(lines[0], lines[1]))
    elif pid == "cons":
        cons, profile = core.consensus(core.parse_fasta(text).values())
        print(cons)
        for b in "ACGT":
            print(f"{b}:", *profile[b])
    elif pid == "grph":
        for a, b in core.overlap_graph(text):
            print(a, b)
    elif pid == "lcsm":
        print(core.longest_shared_motif(core.parse_fasta(text).values()))
    elif pid == "prtm":
        print(core.protein_mass(text))
    elif pid == "revp":
        seq = "".join(core.parse_fasta(text).values()) or text
        for pos, size in core.reverse_palindromes(seq):
            print(pos, size)
    else:
        # Das Problem kennen wir (noch) nicht.
        raise SystemExit(f"Unbekanntes Problem: {pid}")


if __name__ == "__main__":
    main()
