"""
mendelian_error.py: Compute the Mendelian error rate from the VCF files of a
father, mother, and son.

Usage:
    mendelian_error.py <father> <mother> <son>

Arguments:
    <father>    Path to the father's VCF file.
    <mother>    Path to the mother's VCF file.
    <son>       Path to the son's VCF file.

Output:
    The Mendelian error rate (proportion of variants with Mendelian errors).

Example:
    mendelian_error.py father.vcf mother.vcf son.vcf
"""

import sys
import logging
import numpy as np
import pandas as pd

def get_genotype(sample):
    """
    Parse the genotype (GT) field from the SAMPLE column of a VCF file.
    """
    genotype = sample.split(':')[0]

    if genotype == './.':
        return None
    else:
        return genotype


def compute_mendelian_error_rates(father_file, mother_file, child_file):
    """
    Compute the Mendelian error rate from the VCF files of a father, mother,
    and child.
    """
    # Read the VCF files into pandas dataframes
    father_df = pd.read_csv(father_file, sep='\t', comment='#', header=None, \
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                            dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                   'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str}, nrows=10000)
    
    mother_df = pd.read_csv(mother_file, sep='\t', comment='#', header=None, \
                            names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                                dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                    'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str}, nrows=10000)
    
    child_df = pd.read_csv(child_file, sep='\t', comment='#', header=None, \
                            names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], \
                                dtype={'CHROM': str, 'POS': np.int64, 'ID': str, 'REF': str, 'ALT': str, 'QUAL': str, \
                                    'FILTER': str, 'INFO': str, 'FORMAT': str, 'SAMPLE': str}, nrows=10000)

    # Parse the genotype (GT) fields and compute the Mendelian error rates
    total_variants = len(child_df)
    mendelian_errors = 0

    for i in range(total_variants):
        # Loop through the father's variants and compare with the mother's and
        # child's variants

        # Get the current variant's location
        chrom = child_df['CHROM'][i]
        pos = child_df['POS'][i]
        svlen = child_df['INFO'][i].split(';')[0].split('=')[1]

        #print(f"Chrom: {chrom}, Pos: {pos}, SVLEN: {svlen}")

        # Find the same variant in the mother's and father's VCF files
        mother_df = mother_df[(mother_df['CHROM'] == chrom) & (mother_df['POS'] == pos)]
        father_df = father_df[(father_df['CHROM'] == chrom) & (father_df['POS'] == pos)]

        # Check if the variant is present in the mother's and child's VCF files
        if mother_df.empty or father_df.empty:
            #logging.warning("Variant not found in mother's or child's VCF file")
            continue
        else:
            print("Variant found in mother's and father's VCF file at %s:%d" % (chrom, pos))

        # Get the samples
        child_sample = child_df['SAMPLE'][i]
        mother_sample = mother_df['SAMPLE'].values[0]
        father_sample = father_df['SAMPLE'].values[0]

        # Get the genotypes
        father_genotype = get_genotype(father_sample)
        mother_genotype = get_genotype(mother_sample)
        child_genotype = get_genotype(child_sample)
        
        # Skip if any of the genotypes are missing
        if father_genotype is None or mother_genotype is None or child_genotype is None:
            logging.warning("Missing genotype(s) for variant at %s:%d", chrom, pos)
            continue

        # Print the genotypes
        print(f"Father: {father_genotype}, Mother: {mother_genotype}, Child: {child_genotype}")

        # Mendelian error: Child's genotype is inconsistent with inheritance of
        # exactly one allele from each parent.
        # Scenario 1: Father and mother have the same genotype, but the child's
        # genotype is different.
        # Scenario 2: Father and mother have different genotypes, but the
        # child's genotype is the same as one of the parents'.
        # See Smolka et al. (2022) for more details (preprint for Sniffles2):
        # https://www.biorxiv.org/content/10.1101/2022.04.04.487055v2.full
        
        # Scenario 1
        if father_genotype == mother_genotype and father_genotype != son_genotype:
            mendelian_errors += 1

        # Scenario 2
        if father_genotype != mother_genotype and (father_genotype == son_genotype or mother_genotype == son_genotype):
            mendelian_errors += 1

    mendelian_error_rate = mendelian_errors / total_variants
    
    return mendelian_error_rate


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logging.info("Running mendelian_error.py")
    if len(sys.argv) != 4:
        logging.error("Incorrect number of arguments")
        sys.exit(__doc__)

    father_file = sys.argv[1]
    mother_file = sys.argv[2]
    child_file = sys.argv[3]
    me_rate = compute_mendelian_error_rates(father_file, mother_file, child_file)
    logging.info("Mendelian error rate: %.4f", me_rate)
