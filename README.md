# gridss_filter
**gridss_filter** is a framework designed to assist in the prioritization of germline structural variants called by the [GRIDSS](https://github.com/PapenfussLab/gridss) software.

## Prerequisites 
- GRIDSS v.2.13.2
- RepeatMasker v.4.1.5
- R v.4.4.2

Input BAM files must be processed with GRIDSS and annotated using `gridss_annotate_vcf_repeatmasker`.

---

## Installation

All required R packages will be automatically installed if missing:

```
packages <- c("dplyr", "optparse", "stringr", "vcfR", "GenomicRanges", "yaml", "furrr")
````

#### merge_gridss_vcfs.R ###

This script reads all VCFs annotated with gridss_annotate_vcf_repeatmasker and extracts all variants with PASS filter. It generates two dataframes:
- Variants with a single breakend
- Variants with two breakends, merged by ID and MATEID

````r
Rscript merge_gridss_vcfs.R --txt <paths_vcfs.txt> --output <output_folder> [--cores N]
```` 

--txt (-t):
Path to the txt file containing the following information: 

| path  | sample | run | genes.interest |
|-----------|-----------|-----------|-----------|
| path_sample1_vcf_repeatmasker.vcf | sample1  | run1  | gene1, gene2, gene3, gene4|
| path_sample2_vcf_repeatmasker.vcf | sample2  | run1  | gene1, gene3, gene5|

Sample(n) must be replaced by the exact sample name. 
gene(n) must be replaced by the exact gene name. 


--output (-o):
Folder to store the results.


--cores (-c):
The number of cores to use. 




