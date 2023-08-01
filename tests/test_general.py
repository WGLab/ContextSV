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
TEST_SNPS_FILE = os.path.join(TEST_DATA_DIR, 'snps.vcf.gz')

# Set the output directory.
TEST_OUTDIR = os.path.join(os.path.dirname(__file__), 'output')


def test_run():
    """Test the run function."""
    # Run the program.
    contextsv.run(
        TEST_BAM_FILE,
        TEST_SNPS_FILE,
        TEST_OUTDIR,
        "chr3:60380533-60390533"
    )

    # Check that the output file exists.
    output_file = os.path.join(TEST_OUTDIR, 'snp_lrr_baf.csv')
    assert os.path.exists(output_file)

    # Check that the output file is not empty.
    assert os.path.getsize(output_file) > 0

    # Check that the output file has the correct number of lines.
    with open(output_file, 'r') as f:
        assert len(f.readlines()) == 23

    # Check that the output file has the correct header.
    with open(output_file, 'r') as f:
        assert f.readline().strip() == "position,baf,log2_ratio,cnv_state"

    # Check that the output file has the correct last line.
    with open(output_file, 'r') as f:
        last_line = f.readlines()[-1].strip()
        assert last_line == "60389325,0.590909,-0.048852,1"