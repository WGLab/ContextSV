import csv
import sys


def read_tsv(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        return [row for row in reader]

def calculate_mendelian_error(father_genotype, mother_genotype, child_genotype):
    # Generate all possible child genotypes
    child_genotypes = set()
    for allele1 in father_genotype.split('/'):
        for allele2 in mother_genotype.split('/'):
            child_genotypes.add('/'.join(sorted([allele1, allele2])))
    
    # Check if the child genotype is valid
    return 0 if child_genotype in child_genotypes else 1


def main(father_file, mother_file, child_file):
    father_records = read_tsv(father_file)
    mother_records = read_tsv(mother_file)
    child_records = read_tsv(child_file)

    if len(father_records) != len(mother_records) or len(father_records) != len(child_records):
        raise ValueError("All files must have the same number of records")

    total_records = len(father_records)
    error_count = 0

    for i in range(total_records):
        father_genotype = father_records[i][5]
        mother_genotype = mother_records[i][5]
        child_genotype = child_records[i][5]

        error_count += calculate_mendelian_error(father_genotype, mother_genotype, child_genotype)

    error_rate = error_count / total_records
    print(f"Mendelian Inheritance Error Rate: {error_rate:.2%} for {total_records} SVs")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python mendelian_inheritance.py <father_tsv> <mother_tsv> <child_tsv>")
        sys.exit(1)

    father_file = sys.argv[1]
    mother_file = sys.argv[2]
    child_file = sys.argv[3]

    main(father_file, mother_file, child_file)
