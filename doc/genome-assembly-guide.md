#!/bin/bash

# Hybrid Optical Map-Guided Assembly

## Download Test Data
Input data is from GIAB: Oxford Nanopore ultralong (guppy-V3.2.4_2020-01-22)
https://github.com/genome-in-a-bottle/giab_data_indexes

## De novo Assembly
We will perform de novo assembly of the hg002 data using miniasm (instructions from https://github.com/lh3/miniasm)

Also see Optical Kermit's optical map guided assembly script, which uses the
same tools and approach (https://github.com/Denopia/kermit-optical-maps/blob/master/genome-assembly.sh)

### Set up the software environment:

```
export PATH="$software_dir"/minimap2:$PATH
export PATH="$software_dir"/miniasm:$PATH
```

### Use minimap2 to self-align the hg002 data:

```echo "Running minimap2..."
minimap2 -x ava-ont -t 12 input_data/HG002_hs37d5_ONT_GIAB.fastq input_data/HG002_hs37d5_ONT_GIAB.fastq | gzip -1 > output/HG002_hs37d5_ONT_GIAB.paf.gz
```

### Use miniasm to assemble the self-aligned data:

```
echo "Running miniasm..."
miniasm -f inputs/HG002_hs37d5_ONT_GIAB.fastq outputs/HG002_hs37d5_ONT_GIAB.paf.gz > outputs/HG002_hs37d5_ONT_GIAB.gfa
```

### Convert GFA file to FASTA for downstream analysis 
(instructions from https://www.biostars.org/p/169516/)

```
echo "Converting GFA to FASTA..."
awk '/^S/{print ">"$2"\n"$3}' outputs/HG002_hs37d5_ONT_GIAB.gfa > outputs/HG002_hs37d5_ONT_GIAB.fasta
echo "Complete."
```

This yields an assembly which can be used for input into BioNano's hybrid scaffolding pipeline to generate a more accurate hybrid assembly using the BioNano data.
