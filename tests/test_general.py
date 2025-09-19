"""
test_general.py: Test the general module.
"""

import os
import subprocess
import sys
import pytest


# Set the output directory.
TEST_OUTDIR = os.path.join(os.path.dirname(__file__), 'output')

# Check if running unit tests on GitHub Actions
local_dir = os.path.expanduser("~/github/ContextSV")
if os.getcwd() == local_dir:
    TEST_DATA_DIR = os.path.join(local_dir, 'tests/data')
else:
    TEST_DATA_DIR = os.path.abspath(str("tests/data"))

print("Current working directory:", os.getcwd())
print("Test data directory:", TEST_DATA_DIR)
print("Contents of test data directory:", os.listdir(TEST_DATA_DIR) if os.path.exists(TEST_DATA_DIR) else "Directory does not exist")

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
BAM_FILE = os.path.join(TEST_DATA_DIR, 'chr3_test.bam')
REF_FILE = os.path.join(TEST_DATA_DIR, 'GRCh38_noalts_chr3.fa')
SNPS_FILE = os.path.join(TEST_DATA_DIR, 'chr3_test.snps.vcf.gz')
GAP_FILE = os.path.join(TEST_DATA_DIR, 'Gaps-HG38-UCSC-chr3.bed')
HMM_FILE = os.path.join(TEST_DATA_DIR, 'wgs_test.hmm')

# Create a PFB file pointing to the SNP file in the test data directory
PFB_FILE = os.path.join(TEST_DATA_DIR, 'chr3_pfb.txt')
# Remove any existing PFB file to avoid conflicts
if os.path.exists(PFB_FILE):
    os.remove(PFB_FILE)

PFB_SNP_FILE = os.path.join(TEST_DATA_DIR, 'chr3_gnomad_snps.vcf.gz')
with open(PFB_FILE, 'w', encoding='utf-8') as pf:
    pf.write(f"3={PFB_SNP_FILE}\n")

def test_run_help():
    """Ensure the binary runs with --help and exits cleanly."""
    result = subprocess.run(
        ["./build/contextsv", "--help"],
        capture_output=True,
        text=True,
        check=True
    )

    # Print output for debugging purposes
    print("STDOUT:", result.stdout)

    # Check process exited successfully
    assert result.returncode == 0, f"Non-zero exit: {result.returncode}\n{result.stderr}"
    assert result.stdout, "No output from --help"
    assert "Usage:" in result.stdout, "Help text missing 'Usage:'"
    assert "Options:" in result.stdout, "Help text missing 'Options:'"

def test_run_version():
    """Ensure the binary runs with --version and exits cleanly."""
    result = subprocess.run(
        ["./build/contextsv", "--version"],
        capture_output=True,
        text=True,
        check=True
    )

    # Print output for debugging purposes
    print("STDOUT:", result.stdout)

    # Check process exited successfully
    assert result.returncode == 0, f"Non-zero exit: {result.returncode}\n{result.stderr}"
    assert result.stdout, "No output from --version"
    assert "ContextSV version" in result.stdout, "Version text missing 'ContextSV version'"

def test_run_basic():
    """Run ContextSV with basic required parameters."""
    out_dir = os.path.join(TEST_OUTDIR, 'test_run_basic')
    os.makedirs(out_dir, exist_ok=True)

    result = subprocess.run(
        ["./build/contextsv",
         "--bam", BAM_FILE,
         "--ref", REF_FILE,
         "--snp", SNPS_FILE,
         "--outdir", out_dir,
         "--hmm", HMM_FILE,
         "--eth", "nfe",
         "--pfb", PFB_FILE,
         "--sample-size", "20",
         "--min-cnv", "2000",
         "--eps", "0.1",
         "--min-pts-pct", "0.1",
         "--assembly-gaps", GAP_FILE,
         "--chr", "chr3",
         "--save-cnv",
         "--debug"
        ],
        capture_output=True,
        text=True,
        check=True
    )

    # Print output for debugging purposes
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    # Check process exited successfully
    assert result.returncode == 0, f"Non-zero exit: {result.returncode}\n{result.stderr}"
    assert "ContextSV finished successfully!" in result.stdout, "Did not complete successfully"

    # Check for expected output files
    expected_files = [
        os.path.join(out_dir, 'CNVCalls.json'),
        os.path.join(out_dir, 'output.vcf')
    ]
    for ef in expected_files:
        assert os.path.isfile(ef), f"Expected output file not found: {ef}"

    # Find the large duplication in the VCF output
    # chr3	61149366	.	N	<DUP>	.	PASS	END=61925600;SVTYPE=DUP;SVLEN=776235;SVMETHOD=ContextSVv1.0.0-1-g4bd038c;ALN=SPLIT,HMM;HMM=-2533.541937;SUPPORT=63;CLUSTER=23;ALNOFFSET=0;CN=6	GT:DP	1/1:63
    vcf_file = os.path.join(out_dir, 'output.vcf')
    found_dup = False
    with open(vcf_file, 'r') as vf:
        for line in vf:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom, pos, id, ref, alt, qual, filter, info = fields[:8]

            # Print each VCF record for debugging
            print(f"VCF Record: {chrom}, {pos}, {id}, {ref}, {alt}, {qual}, {filter}, {info}")

            if (chrom == 'chr3' and pos == '61149366' and alt == '<DUP>'):
                found_dup = True
                info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
                assert info_dict.get('SVTYPE') == 'DUP', "SVTYPE is not DUP"
                assert int(info_dict.get('SVLEN', 0)) == 776235, f"SVLEN is not 776235, got {info_dict.get('SVLEN')}"
                assert int(info_dict.get('CN', 0)) == 6, f"CN is not 6, got {info_dict.get('CN')}"
                break
    assert found_dup, "Expected duplication not found in VCF output"

    # Find the large duplication in the CNVCalls.json output
    # [
        # {
        # "chromosome": "chr3",
        # "start": 61149366,
        # "end": 61925600,
        # "sv_type": "DUP",
        # "likelihood": -2533.54,
        # "size": 776235,

    json_file = os.path.join(out_dir, 'CNVCalls.json')
    found_dup_json = False
    import json
    with open(json_file, 'r') as jf:
        cnv_data = json.load(jf)
        for entry in cnv_data:
            if (entry.get('chromosome') == 'chr3' and
                entry.get('start') == 61149366 and
                entry.get('end') == 61925600 and
                entry.get('sv_type') == 'DUP'):
                found_dup_json = True
                assert abs(entry.get('likelihood', 0) + 2533.54) < 0.1, f"Likelihood is not -2533.54, got {entry.get('likelihood')}"
                assert entry.get('size') == 776235, f"Size is not 776235, got {entry.get('size')}"
                break
    assert found_dup_json, "Expected duplication not found in CNVCalls.json output" 
