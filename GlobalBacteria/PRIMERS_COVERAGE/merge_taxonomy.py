#!/usr/bin/env python3
import glob
import os
import csv

files = sorted(glob.glob("*_taxonomy_coverage.txt"))

if not files:
    raise SystemExit("Nebyly nalezeny žádné *_taxonomy_coverage.txt soubory.")

data = {}
all_taxonomies = set()

for fname in files:
    with open(fname, "r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue

            parts = line.rsplit(None, 1)
            if len(parts) != 2:
                print(f"Přeskakuji nevalidní řádek v {fname}: {line}")
                continue

            taxonomy, count = parts
            try:
                count = int(count)
            except ValueError:
                print(f"Přeskakuji nevalidní count v {fname}: {line}")
                continue

            all_taxonomies.add(taxonomy)
            if taxonomy not in data:
                data[taxonomy] = {}
            data[taxonomy][os.path.basename(fname)] = count

output_file = "merged_taxonomy_coverage.tsv"

with open(output_file, "w", encoding="utf-8", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    header = ["taxonomy"] + [os.path.basename(f) for f in files]
    writer.writerow(header)

    for taxonomy in sorted(all_taxonomies):
        row = [taxonomy]
        for fname in files:
            row.append(data.get(taxonomy, {}).get(os.path.basename(fname), 0))
        writer.writerow(row)

print(f"Hotovo. Výstup uložen do: {output_file}")
