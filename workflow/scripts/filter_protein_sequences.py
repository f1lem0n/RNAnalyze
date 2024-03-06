from Bio import SeqIO

if __name__ == "__main__":
    with open(snakemake.input[0], "r") as f:
        records = list(SeqIO.parse(f, "fasta"))
    filtered = []
    for r in records:
        if set(r.seq) <= set("AUGC"):
            continue
        elif set(r.seq) <= set("ATGC"):
            continue
        else:
            filtered.append(r)
    with open(snakemake.output[0], "w") as f:
        SeqIO.write(filtered, f, "fasta")
