from os import walk
from Bio import SeqIO


def split_sequence(
    sequence: str, threshold: int, chunk_size: int
) -> list | str:
    """
    Split a sequence into chunks of a given size if sequence is
    longer than threshold, else return unchanged sequence
    """
    if threshold < 0 or chunk_size < 0:
        return [sequence]
    if not len(sequence) >= threshold:
        return [sequence]
    n = len(sequence) // chunk_size
    chunks = [
        sequence[i * chunk_size : (i + 1) * chunk_size] for i in range(n)
    ]
    chunks.append(
        sequence[n * chunk_size : n * chunk_size + len(sequence) % chunk_size]
    )
    return chunks


def get_sequence_filepaths(directory: str) -> list:
    """
    Get all filepaths of sequence files in a directory
    """
    sequence_filepaths = []
    for root, _, files in walk(directory):
        for file in files:
            if file.endswith(".fasta") and "combined" not in file:
                sequence_filepaths.append(f"{root}/{file}")
    return sequence_filepaths


def write_chunks_to_files(chunks: list, rec, directory: str) -> None:
    for idx, chunk in enumerate(chunks):
        temp_rec = SeqIO.SeqRecord(seq=chunk)
        temp_rec.id = f"{rec.id}_chunk{idx}"
        temp_rec.name = f"{rec.name}_chunk{idx}"
        temp_rec.description = f"{rec.description}"
        with open(f"{directory}/{temp_rec.id}.fasta", "w") as handle:
            SeqIO.write(temp_rec, handle, "fasta")


def main():
    input_directory = snakemake.input[0]
    output_directory = snakemake.output.dir
    sequence_filepaths = get_sequence_filepaths(input_directory)
    for sf in sequence_filepaths:
        with open(sf, "r") as handle:
            rec = list(SeqIO.parse(handle, "fasta"))[0]
        chunks = split_sequence(
            rec.seq,
            snakemake.params.threshold,
            snakemake.params.chunk_size,
        )
        write_chunks_to_files(chunks, rec, output_directory)
    with open(f"{output_directory}/.done", "w") as handle:
        handle.write("")


if __name__ == "__main__":
    main()
