# 🧬 Rosalind Solutions

Solutions to bioinformatics problems from [rosalind.info](https://rosalind.info), written in clean, tested Python.
Code comments are written in **German (A2 level)** – I use this project to practice the language, too.

## Solved problems

| ID   | Topic                                | Function                    |
|------|--------------------------------------|-----------------------------|
| DNA  | Counting nucleotides                 | `count_nucleotides`         |
| RNA  | Transcription                        | `transcribe`                |
| REVC | Reverse complement                   | `reverse_complement`        |
| FIB  | Recurrence relations (rabbits)       | `rabbits`                   |
| FIBD | Mortal rabbits (dynamic programming) | `mortal_rabbits`            |
| GC   | GC content + FASTA parsing           | `highest_gc`                |
| HAMM | Hamming distance                     | `hamming`                   |
| IPRB | Mendel's first law (probability)     | `dominant_probability`      |
| PROT | RNA → protein translation            | `translate`                 |
| SUBS | Motif search                         | `find_motif`                |
| CONS | Consensus & profile matrix           | `consensus`                 |
| GRPH | Overlap graphs                       | `overlap_graph`             |
| LCSM | Longest shared motif                 | `longest_shared_motif`      |
| PRTM | Protein mass                         | `protein_mass`              |
| REVP | Restriction sites (palindromes)      | `reverse_palindromes`       |

## Usage

```bash
pip install pytest
python -m pytest            # runs all sample datasets
python solve.py revc data/rosalind_revc.txt
```

## Skills
Python · string algorithms · dynamic programming · probability · graph basics · unit testing
