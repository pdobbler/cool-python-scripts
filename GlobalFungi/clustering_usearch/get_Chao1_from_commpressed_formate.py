import csv
import sys
from collections import defaultdict

def compute_chao1(otu_abundances):
    S_obs = len(otu_abundances)
    f1 = sum(1 for v in otu_abundances.values() if v == 1)
    f2 = sum(1 for v in otu_abundances.values() if v == 2)
    if f2 == 0:
        chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1)) if f1 > 1 else S_obs
    else:
        chao1 = S_obs + (f1 ** 2) / (2 * f2)
    return S_obs, f1, f2, chao1

def main(input_file, output_file):
    sample_data = defaultdict(lambda: defaultdict(int))

    with open(input_file, "r") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            sample = row["sample"]
            otus = row["OTUs"].split(";")
            abundances = list(map(int, row["abundance"].split(";")))
            for otu, count in zip(otus, abundances):
                sample_data[sample][otu] += count

    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["sample", "S_obs", "F1", "F2", "Chao1"])
        for sample, otu_abundances in sample_data.items():
            S_obs, F1, F2, chao1 = compute_chao1(otu_abundances)
            writer.writerow([sample, S_obs, F1, F2, round(chao1, 2)])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python chao1_calc.py input_file.tsv output_file.tsv")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
