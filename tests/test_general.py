"""
test_general.py: Test the general module.
"""

import os
import sys
import pytest

# Import lib from the parent directory
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from lib import contextsv

# Get the path to the test data directory.
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

# Get the path to the test data files.
TEST_BAM_FILE = os.path.join(TEST_DATA_DIR, 'test.bam')
TEST_REF_FILE = ""  # Empty string for now.
TEST_SNPS_FILE = os.path.join(TEST_DATA_DIR, 'snps.vcf.gz')
TEST_PFB_FILE = ""

# Set the output directory.
TEST_OUTDIR = os.path.join(os.path.dirname(__file__), 'output')

# Test the main function.
def test_run():
    
    # Set input parameters.
    input_data = contextsv.InputData()
    input_data.setShortReadBam(TEST_BAM_FILE)
    input_data.setLongReadBam(TEST_BAM_FILE)
    input_data.setRefGenome(TEST_REF_FILE)
    input_data.setSNPFilepath(TEST_SNPS_FILE)
    input_data.setRegion("chr3:60380533-60390533")
    input_data.setThreadCount(1)
    input_data.setMeanChromosomeCoverage("chr3:39.561,chr6:39.4096")
    input_data.setAlleleFreqFilepaths("")
    input_data.setHMMFilepath("")
    input_data.setOutputDir(TEST_OUTDIR)

    # Run the analysis.
    contextsv.run(input_data)

    # Check that the output file exists.
    output_file = os.path.join(TEST_OUTDIR, 'cnv_data.tsv')
    assert os.path.exists(output_file)

    # Check that the output file is not empty.
    assert os.path.getsize(output_file) > 0

    # Check that the output file has the correct number of lines.
    with open(output_file, 'r') as f:
        assert len(f.readlines()) == 23

    # Check that the output file has the correct header.
    with open(output_file, 'r') as f:
        assert f.readline().strip() == "chromosome\tposition\tb_allele_freq\tlog2_ratio\tcnv_state\tpopulation_freq"

    # Check that the output file has the correct SNP values in the last line
    with open(output_file, 'r') as f:
        last_line = f.readlines()[-1].strip('\n')
        print(last_line)
        actual_line="chr3\t60389325\t0.590909\t-0.048852\t6\t0.01"
        print(actual_line)
        assert last_line == actual_line

