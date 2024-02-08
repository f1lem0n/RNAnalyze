from subprocess import run
from zipfile import ZipFile
from os import walk, path


def run_structural_homology_analysis(pdb_id):
    run(
        f"cd structural_homology && ./get_homologs {pdb_id}",
        shell=True,
        check=True,
    )
    run(f"cd structural_homology && ./align {pdb_id}", shell=True, check=True)


def run_sequence_homology_analysis(pdb_id):
    run(
        f"cd sequence_homology && ./get_homologs {pdb_id}",
        shell=True,
        check=True,
    )
    # run(f"cd sequence_homology && ./align {pdb_id}", shell=True, check=True)


def get_all_file_paths(directory):
    file_paths = []

    for root, _, files in walk(directory):
        for filename in files:
            filepath = path.join(root, filename)
            file_paths.append(filepath)

    return file_paths


def move_zip_output():
    run("mkdir -p output/structural_homology", shell=True, check=True)
    run("mkdir -p output/sequence_homology", shell=True, check=True)
    run(
        "mv structural_homology/data/* output/structural_homology/.",
        shell=True,
        check=True,
    )
    run(
        "mv sequence_homology/data/* output/sequence_homology/.",
        shell=True,
        check=True,
    )
    with ZipFile("output.zip", "w") as zipf:
        for file in get_all_file_paths("output"):
            zipf.write(file)


def main():
    run("echo $PATH", shell=True, check=True)
    pdb_id = input("Enter PDB ID: ").lower()
    run_structural_homology_analysis(pdb_id)
    run_sequence_homology_analysis(pdb_id)
    move_zip_output()


if __name__ == "__main__":
    main()
