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


def test_run():
    """Test the run function of the contextsv module with a small dataset."""
    
    # Set input parameters.
    input_data = contextsv.InputData()
    input_data.setShortReadBam(TEST_BAM_FILE)
    input_data.setLongReadBam(TEST_BAM_FILE)
    input_data.setRefGenome(TEST_REF_FILE)
    input_data.setSNPFilepath(TEST_SNPS_FILE)
    input_data.setChromosome("21")
    input_data.setRegion("14486099-14515105")
    input_data.setThreadCount(1)
    input_data.setAlleleFreqFilepaths("")
    input_data.setHMMFilepath("")
    input_data.setOutputDir(TEST_OUTDIR)
    input_data.saveCNVData(True)
    
    # Run the analysis.
    contextsv.run(input_data)

    # Check that the output file exists.
    output_file = os.path.join(TEST_OUTDIR, 'output.vcf')
    assert os.path.exists(output_file)

    # Check that the VCF file is not empty.
    assert os.path.getsize(output_file) > 0

    # Check that the VCF file has the correct number of lines.
    with open(output_file, 'r', encoding='utf-8') as f:
        assert len(f.readlines()) == 25

    # Check that the VCF file has the correct header, and the correct
    # VCF CHROM, POS, and INFO fields in the next 2 lines.
    header_line = 18
    with open(output_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            if i == header_line:
                assert line.strip() == "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
            elif i == header_line + 1:
                # Get the first, second, and eighth fields from the line
                fields = line.strip().split('\t')
                assert fields[0] == "21"
                assert fields[1] == "14458394"
                assert fields[7] == "END=14458394;SVTYPE=INS;SVLEN=1341;SUPPORT=1;SVMETHOD=CONTEXTSVv0.1;ALN=CIGARINS,;CLIPSUP=0;REPTYPE=NA;HMM=0.000000"
            elif i == header_line + 2:
                fields = line.strip().split('\t')
                assert fields[0] == "21"
                assert fields[1] == "14458394"
                assert fields[7] == "END=14458394;SVTYPE=INS;SVLEN=1344;SUPPORT=1;SVMETHOD=CONTEXTSVv0.1;ALN=CIGARINS,;CLIPSUP=0;REPTYPE=NA;HMM=0.000000"
                break
            