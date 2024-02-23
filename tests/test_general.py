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
TEST_REF_FILE = ""
TEST_SNPS_FILE = os.path.join(TEST_DATA_DIR, 'snps.vcf.gz')
TEST_PFB_FILE = ""

# Set the output directory.
TEST_OUTDIR = os.path.join(os.path.dirname(__file__), 'output')

# Check if running unit tests on GitHub Actions
local_dir = os.path.expanduser("~/github/ContextSV")
if os.getcwd() == local_dir:
    test_dir_path = os.path.join(local_dir, 'tests/data')
else:
    test_dir_path = os.path.abspath(str("TestData"))

TEST_REF_FILE = os.path.join(test_dir_path, 'GRCh37_21.fa')
TEST_BAM_FILE = os.path.join(test_dir_path, 'ONT_Kit14_21.bam')
TEST_SNPS_FILE = os.path.join(test_dir_path, 'snps_21.vcf.gz')


# Test the main function.
def test_run():
    
    # Set input parameters.
    input_data = contextsv.InputData()
    input_data.setShortReadBam(TEST_BAM_FILE)
    input_data.setLongReadBam(TEST_BAM_FILE)
    input_data.setRefGenome(TEST_REF_FILE)
    input_data.setSNPFilepath(TEST_SNPS_FILE)
    input_data.setRegion("21:14486099-14515105")
    input_data.setThreadCount(1)
    input_data.setMeanChromosomeCoverage("21:80.6292")
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
        assert len(f.readlines()) == 55

    # Check that the output file has the correct header.
    with open(output_file, 'r') as f:
        assert f.readline().strip() == "chromosome\tposition\tsnp\tb_allele_freq\tlog2_ratio\tcnv_state\tpopulation_freq"

    # Check that the output file has the correct SNP values in the last line
    with open(output_file, 'r', encoding='utf-8') as f:
        last_line = f.readlines()[-1].strip('\n')
        print("The last line of the output file is: ")
        print(last_line)
        actual_line="21\t14508888\t0\t0.5\t0\t6\t0.01"
        print(actual_line)
        assert last_line == actual_line
