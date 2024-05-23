import primer3
import sys

original_primers = sys.argv[1]
out_info = sys.argv[2]

def check_self_dimer(primer_sequence):
    # Check for self-dimer using the updated calc_hairpin function
    dimer_result = primer3.calc_hairpin(primer_sequence)
    if dimer_result.structure_found:
        return {
            "dG": dimer_result.dg,  # Free energy
            "Tm": dimer_result.tm,  # Melting temperature
            "structure_found": True,
        }
    else:
        return {
            "structure_found": False,
        }

# 528F_H001 TGTATAGCGGTAATTCCAGCTCCAA
out_file = open(out_info, "w")
out_file.write("name\tseq\tdG\tTm\n")
for line in open(original_primers):
    vals = line.rstrip().split('\t')
    if vals[0] != "":
        line = vals[0] + "\t" + vals[1]
        seq = vals[1]
        result = check_self_dimer(seq)
        if result["structure_found"]:
            #print(f"Self-dimer found with ΔG: {result['dG']} and Tm: {result['Tm']}")
            line += "\t" + str(result['dG']) + "\t" + str(result['Tm'])
        else:
            #print("No self-dimer found.")
            line += "\tOK\tOK"
        out_file.write(line + "\n")

out_file.close()


