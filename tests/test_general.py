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

# Set the output directory.
TEST_OUTDIR = os.path.join(os.path.dirname(__file__), 'output')


def test_run():
    """Test the run function."""
    # Run the program.
    contextsv.run(
        TEST_BAM_FILE,
        TEST_REF_FILE,
        TEST_SNPS_FILE,
        TEST_OUTDIR,
        "chr3:60380533-60390533"
    )

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
        assert f.readline().strip() == "chromosome\tposition\tb_allele_freq\tlog2_ratio\tcnv_state"

    # Check that the output file has the correct SNP values (excluding predicted
    # state) in the last line
    with open(output_file, 'r') as f:
        last_line = f.readlines()[-1]

        print("Debugging:")
        print(last_line)
        print(last_line[:-2])

        # CNV2:
        #chr3	60389325	0.590909	-0.048852
        assert last_line[:-2] == "chr3\t60389325\t0.590909\t-0.048852\t"

        # CNV3:
        #assert last_line[:-2] == "60389325,0.590909,-0.0433203"
