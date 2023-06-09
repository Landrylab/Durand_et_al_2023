# Scripts used in Durand et al, 2023

The repo contains all scripts used in the study Durand et al (2023). Raw data are provided (except pictures taken by the S & P robot) so that virtually all analyses and graphs can theoretically be reproduced.

## S & P analyses

Pictures taken by the S & P robot are not provided (too many heavy files). Below, we do provide the script used to process them and obtain the provided csv files.

1. Pictures were cropped and converted into inverted levels of gray using [this script](image_conversion.py)
2. Pictures were analyzed using `pyphe-quantify` from the [`pyphe` toolbox](https://github.com/Bahler-Lab/pyphe)
```
# --grid needs to be adapted depending on the size of the array (384, 1536)
# Command line is launched from the folder containing pictures to be analyzed
pyphe-quantify batch --grid auto_384 --pattern "*.JPG" --s 0.1
```
3. CSV files generated by pyphe-quantify were analyzed using [this notebook](robotpics_analysis_manuscript_edition.ipynb)

## Analyses of growth assays in plate readers
- [Growth of *FUR1* mutants](20230324_growthcurves_FUR1.ipynb)
- [Growth of mutants in liquid medium](growth_curves_TECAN384.ipynb)
- [Medium conditioning assay](20230331_medium_conditioning.ipynb)
- [Growth of LL13-040 strain in defined media](20230401_growthcurves_SC-SD.ipynb)
- [5-FC dose-response](20230331_dose-response.ipynb)
- [Growth in minimal medium with varying nitrogen sources](20230411_growthcurves_cytosine.ipynb)

## Analyses of cytometry data

- [Rhodamine accumulation assay](rhodamine.ipynb)
- [Fluorescence reporter assay](cytometry-FCY1.ipynb)

## Analyses of WGS data
- [Snakemake pipeline to analyze demultiplexed fastq files](Snakefile)
- [SNP calling and filtering with 'samtools'](README_bcftools.ipynb)
- [SNP calling and filtering with 'gatk'](README_gatk.ipynb)
- [Mutations detected by samtools/gatk/Sanger](FUR1_variants_v2.ipynb)
- [Coverage](coverage.ipynb)
