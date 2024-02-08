from sys import argv
from subprocess import run


def run_structural_homology_analysis(pdb_id):
    run(f"./structural_homology/get_homologs {pdb_id}", shell=True, check=True)
    run(f"./structural_homology/align {pdb_id}", shell=True, check=True)


def main():
    pdb_id = argv[1]


if __name__ == "__main__":
    main()
