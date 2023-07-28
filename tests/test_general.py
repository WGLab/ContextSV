"""
test_general.py: Test the general module.
"""

import os
import pytest

from lib import contextsv

# Get the path to the test data directory.
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')

# Get the path to the test data files.
TEST_BAM_FILE = os.path.join(TEST_DATA_DIR, 'test.bam')
TEST_REF_FILE = os.path.join(TEST_DATA_DIR, 'hg19.fa')
TEST_SNPS_FILE = os.path.join(TEST_DATA_DIR, 'snps.vcf.gz')



