from Bio import Entrez, SeqIO
from pathlib import Path
from subprocess import run
from time import sleep


def parse_accession_numbers(input: Path) -> list:
    segregated_accessions = {}
    with open(input, "r") as f:
        table = f.readlines()
    for row in table:
        pdb_id, ncbi_accession = row.split("\t")[:2]
        pdb_id = pdb_id.split("_")[0].lower()
        if pdb_id not in segregated_accessions:
            segregated_accessions[pdb_id] = [ncbi_accession]
        else:
            segregated_accessions[pdb_id].append(ncbi_accession)
    return segregated_accessions


def download_sequence(accession_number: str, output: Path) -> None:
    print(f"Downloading {accession_number} --> {output}")
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
        shell=True,
    )


def batch_download(segregated_accessions: dict, output: Path) -> None:
    # Batch download sequences
    for homolog, accessions in segregated_accessions.items():
        if homolog not in str(output):
            continue
        for a in set(accessions):
            download_sequence(a, output)
        combine_sequences(output)

    # if there are no sequences for a homolog, create placeholder files
    # this is crucial for the pipeline to not fail on subsequent rules
    for homolog, _ in segregated_accessions.items():
        if homolog not in list(output.iterdir()):
            (output / homolog).mkdir()
            with open(output / "combined.fasta", "w") as f:
                f.write("")


def main():
    Entrez.email = snakemake.config["email"]
    input = Path(snakemake.input[0])
    output = Path(snakemake.output.dir)
    segregated_accessions = parse_accession_numbers(input)
    batch_download(segregated_accessions, output)


if __name__ == "__main__":
    main()
