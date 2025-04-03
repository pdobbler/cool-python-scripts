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

    with open(input_file, "r") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            if line.strip() == "":
                continue
            parts = line.strip().split("\t")
            sample = parts[0]
            otus = parts[1].split(";")
            abundances = list(map(int, parts[2].split(";")))
            for otu, count in zip(otus, abundances):
                sample_data[sample][otu] += count

    with open(output_file, "w") as f:
        f.write("sample\tS_obs\tF1\tF2\tChao1\n")
        for sample, otu_abundances in sample_data.items():
            S_obs, F1, F2, chao1 = compute_chao1(otu_abundances)
            f.write(f"{sample}\t{S_obs}\t{F1}\t{F2}\t{round(chao1, 2)}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python chao1_calc.py input_file.tsv output_file.tsv")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
