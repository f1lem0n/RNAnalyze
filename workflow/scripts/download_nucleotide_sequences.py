from sys import argv
from pathlib import Path
from Bio import Entrez, SeqIO
from time import sleep


def main():
    accession_numbers = list(set(argv[2:]))
    output = Path(argv[1])
    for accession_number in accession_numbers:
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


if __name__ == "__main__":
    main()
