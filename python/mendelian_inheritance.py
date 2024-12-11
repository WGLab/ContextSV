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

    # Print the parent and child genotypes if invalid
    if child_genotype not in child_genotypes:
        print(f"ME: Father: {father_genotype}, Mother: {mother_genotype}, Child: {child_genotype}")
    
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

    sv_type_dict = {}
    sv_type_error_dict = {}

    for i in range(total_records):
        father_genotype = father_records[i][5]
        mother_genotype = mother_records[i][5]
        child_genotype = child_records[i][5]
        child_sv_type = child_records[i][2]
        sv_type_dict[child_sv_type] = sv_type_dict.get(child_sv_type, 0) + 1

        # Print SV size if error occurs
        error_value = calculate_mendelian_error(father_genotype, mother_genotype, child_genotype)
        if error_value == 1:
            # print(f"SV size: {father_records[i][2]}")
            sv_type_error_dict[child_sv_type] = sv_type_error_dict.get(child_sv_type, 0) + 1

        error_count += error_value
        # error_count += calculate_mendelian_error(father_genotype, mother_genotype, child_genotype)

    if total_records == 0:
        error_rate = 0
        print("No records found")
    else:
        error_rate = error_count / total_records

    print(f"Mendelian Inheritance Error Rate: {error_rate:.2%} for {total_records} shared trio SVs")

    print("SV Type Distribution:")
    for sv_type, count in sv_type_dict.items():
        error_count = sv_type_error_dict.get(sv_type, 0)
        error_rate = error_count / count
        print(f"{sv_type}: {error_rate:.2%} ({error_count}/{count})")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python mendelian_inheritance.py <father_tsv> <mother_tsv> <child_tsv>")
        sys.exit(1)

    father_file = sys.argv[1]
    mother_file = sys.argv[2]
    child_file = sys.argv[3]

    main(father_file, mother_file, child_file)
