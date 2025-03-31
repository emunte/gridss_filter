# gridss_filter
gridss_filter is a framework to assist in the priorization of germline_variants called by GRIDSS software (https://github.com/PapenfussLab/gridss).

### Prerequisites ###
- GRIDSS v.2.13.2
- RepeatMasker v.4.1.5
- R v.4.4.2

Your bams of interest must be analyzed using gridss and annotated using gridss_annotate_vcf_repeatmasker.


### Running ####

#### merge_gridss_vcfs.R ###

 ````
Rscript merge_gridss_vcfs.R --txt <paths_vcfs.txt> --output <path_folder_output> [--cores n]
```` 

path	sample	run	genes.interest

| path  | sample | run | genes.interest |
|-----------|-----------|-----------|-----------|
| path_sample1_vcf_repeatmasker.vcf | sample1  | run 1  | gene1, gene2, gene3, gene4|
| path_sample2_vcf_repeatmasker.vcf | sample2  | run 1  | gene1, gene3, gene5, gene6|
