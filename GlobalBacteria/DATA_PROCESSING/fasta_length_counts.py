#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path


def open_maybe_gzip(file_path):
    """
    Otevře běžný nebo .gz soubor v textovém režimu.
    """
    file_path = str(file_path)
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt", encoding="utf-8")
    return open(file_path, "r", encoding="utf-8")


def read_fasta_length_counts(fasta_path):
    """
    Načte FASTA / FASTA.GZ a vrátí slovník:
    {delka_sekvence: pocet_sekvenci}
    """
    length_counts = {}

    current_name = None
    current_seq_len = 0

    with open_maybe_gzip(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    length_counts[current_seq_len] = length_counts.get(current_seq_len, 0) + 1

                current_name = line[1:].split()[0]
                current_seq_len = 0
            else:
                current_seq_len += len(line)

        if current_name is not None:
            length_counts[current_seq_len] = length_counts.get(current_seq_len, 0) + 1

    return length_counts


def write_tsv(length_counts, output_path):
    """
    Zapíše TSV ve formátu:
    length    count
    """
    with open(output_path, "w", encoding="utf-8") as out:
        out.write("length\tcount\n")
        for length in sorted(length_counts):
            out.write(f"{length}\t{length_counts[length]}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Spočítá počty sekvencí podle délky z FASTA / FASTA.GZ a uloží je do TSV."
    )
    parser.add_argument("fasta", help="Vstupní FASTA nebo FASTA.GZ soubor")
    parser.add_argument(
        "-o", "--output",
        default="length_counts.tsv",
        help="Výstupní TSV soubor"
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise FileNotFoundError(f"Soubor neexistuje: {fasta_path}")

    length_counts = read_fasta_length_counts(fasta_path)

    if not length_counts:
        raise ValueError("Ve FASTA souboru nebyly nalezeny žádné sekvence.")

    write_tsv(length_counts, args.output)
    print(f"Uloženo do: {args.output}")


if __name__ == "__main__":
    main()
