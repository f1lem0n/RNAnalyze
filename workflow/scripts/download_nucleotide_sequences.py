from Bio import Entrez, SeqIO
from pathlib import Path
from subprocess import run
from time import sleep

def parse_accession_numbers(input: Path) -> list:
    return list(set(run(
        "awk '{print $2}' " + str(input),
        shell=True,
        capture_output=True,
        text=True
    ).stdout.strip().splitlines()))

def download_sequence(accession_number: str, output: Path) -> None:
    print(f"Downloading {accession_number}...")
    with Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        id=accession_number,
    ) as handle:
        record = SeqIO.read(handle, "fasta")
        SeqIO.write(record, output / f"{accession_number}.fasta", "fasta")
    sleep(1)

def combine_sequences(output: Path) -> None:
    print("Combining sequences...")
    run(
        f"cat {output}/*.fasta > {output}/combined.fasta",
        shell=True
    )

def main():
    Entrez.email = snakemake.config["email"]
    input = Path(snakemake.input[0])
    output = Path(snakemake.output.dir)
    accession_numbers = parse_accession_numbers(input)
    for accession_number in accession_numbers:
        download_sequence(accession_number, output)
    combine_sequences(output)


if __name__ == "__main__":
    main()
