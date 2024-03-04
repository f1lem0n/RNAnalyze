from Bio import Entrez, SeqIO
from pathlib import Path
from subprocess import run
from time import sleep


def batch_download(input: Path, output: Path) -> None:
    with open(input, "r") as f:
        blast_table = f.readlines()
    for line in blast_table:
        accession_number = line.split("\t")[1]
        homolog = line.split("\t")[0].split("_")[0].lower()
        start = int(line.split("\t")[8])
        stop = int(line.split("\t")[9])
        if homolog in str(output):
            download_sequence(accession_number, start, stop, output)
    combine_sequences(output)


def download_sequence(
    accession_number: str, start: int, stop: int, output: Path
) -> None:
    with Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        id=accession_number,
    ) as handle:
        record = SeqIO.read(handle, "fasta")
        record.seq = record.seq[start - 1 : stop]
        SeqIO.write(record, output / f"{accession_number}.fasta", "fasta")
    sleep(1)


def combine_sequences(output: Path) -> None:
    run(
        f"cat {output}/*.fasta > {output}/combined.fasta 2> {snakemake.log} "
        f"|| touch {output}/combined.fasta",
        shell=True,
    )


def main():
    Entrez.email = snakemake.config["email"]
    batch_download(Path(snakemake.input[0]), Path(snakemake.output[0]))


if __name__ == "__main__":
    main()
