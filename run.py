from subprocess import run
from os import environ


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


def main():
    pdb_id = environ["PDB_ID"]
    run_structural_homology_analysis(pdb_id)
    # run_sequence_homology_analysis(pdb_id)


if __name__ == "__main__":
    main()
