# DX_UPD

Software to parse VCF files of family (full trio only with child, father, mother) and create IGV track with paternal or maternally inherited loci.

## Usage
```bash
usage: make_UPD_igv.py [-h] [-c] [--mindepth MINDEPTH] [--suffix SUFFIX]
                       [--maxlocus MAXLOCUS] [--includehet] [--includenormal]
                       ped_file run_id sample_id input_files [input_files ...]

positional arguments:
  ped_file             PED file
  run_id               run ID prefix
  sample_id            sample id to be processed
  input_files          input files (space separated)

optional arguments:
  -h, --help           show this help message and exit
  -c, --compressed     VCF input is compressed (.gz)
  --mindepth MINDEPTH  Threshold for minimum depth (DP) of SNV (default = 15)
  --suffix SUFFIX      suffix of VCF file to be searched (default = .vcf)
  --maxlocus MAXLOCUS  maximum size of locus to be printed. This reduces large
                       blocks in regions with low informativity (default =
                       50000)
  --includehet         Include (possible) heterodisomy SNVs
  --includenormal      Include normal inherited SNVs
``` 

## Installation
To run the UPD scripts we need to create a virtual python environment

### Make virtual python environment (note this could be UMCU-HPC specific)
```bash
python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
All scripts are tested using Python 3.6.8
