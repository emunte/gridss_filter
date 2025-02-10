################################################################################
## Filter gridss for germline variants
################################################################################


# 1 . Load libraries----

#Function to install libraries if necessary
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressMessages(library(pkg, character.only = TRUE, quietly = TRUE))
}

#List of libraries
packages <- c("dplyr", "optparse", "stringr", "GenomicRanges", "yaml", "DiagrammeR", "DiagrammeRsvg", "rsvg")


#Apply
invisible(sapply(packages, install_if_missing))



# 2. Build options list

option_list <- list(
  optparse::make_option(c("-i", "--input"), type="character", default="",
                        help="Folder with merge results", metavar="character"),
  optparse::make_option(c("-d", "--bedFile"), type="character", default="",
                       help="path to the bed file", metavar="character"),
  optparse::make_option(c("-p", "--params"), type="character", default="params.yaml",
                        help="Yaml file with all the specifications", metavar="integer"),
  optparse::make_option(c("-o", "--output"), type="character", default="",
                        help="Path to the output directory. Results will be stored in this folder", metavar="character")
)

opt_parser <- optparse::OptionParser(option_list=option_list);
# Load params
args <- optparse::parse_args(opt_parser)
input <- args$txt
output <- args$output
bedFile <- args$bedFile %>% rtracklayer::import()
params <- yaml::yaml.load(args$params)



# 3. Load data

#Get the two dataframes (merge_gridss_vcf output)
all.one.break <- list.files(input, pattern="merge_gridss_one_break.txt", recursive=TRUE, full.names=TRUE) %>%
  read.delim(header=TRUE, sep="\t") %>%
  dplyr::mutate(end=POS)


all.two.break <- list.files(input, pattern="merge_gridss_two_break_same_chrom.txt", recursive=TRUE, full.names=TRUE) %>%
  read.delim(header=TRUE, sep="\t") %>%  #reorder
  dplyr::mutate(
    chrom1=ifelse(POS.y < POS.x, CHROM.y, CHROM.x),
    chrom2=ifelse(POS.y < POS.x, CHROM.x, CHROM.y),
    start= ifelse(POS.x<POS.y, POS.x, POS.y),
    end= ifelse(POS.y>POS.x, POS.y, POS.x),
    alt1= ifelse(POS.x<POS.y, ALT.x, ALT.y),
    alt2= ifelse(POS.y>POS.x, ALT.y, ALT.x))



# 4. Pipeline----
#register number of initial variants---
filters.variants <- data.frame(type=c("one_break", "two_break"),totals=c(nrow(all.one.break), nrow(all.two.break)) )



## A) Frequency filters----
#Delete variants found in >10 samples

#one-break
count.obs <-  all.one.break %>%
  dplyr::count(CHROM, POS, ALT, name = "n") #count same variant (CHROM; POS; ALT)
count.obs2 <-  all.one.break %>%
  dplyr::count(CHROM, POS,  name = "n2") #count similar variants (CHROM; POS)

overlaps.one.sample.count <- all.one.break %>%  #merge results
  dplyr::left_join(count.obs2, by = c("CHROM", "POS")) %>%
  dplyr::left_join(count.obs, by=c("CHROM", "POS", "ALT"))%>%
  dplyr::mutate(start = POS) %>%
  dplyr::filter(n <= frequency.filter) %>% #filter
  dplyr::relocate(CHROM, start, end, run, sample, genes.interest, .before = POS)

#two-breaks
count.obs.two <- all.two.break %>%
  dplyr::count(chrom1,chrom2, start, end, alt1, alt2, name = "n") #count same variant (CHROM; POS; ALT)
count.obs.two2 <- all.two.break %>%
  dplyr::count(chrom1, chrom2, start, end, name = "n2") #count similar variants (CHROM; POS)

overlaps.two.sample.count <- all.two.break %>%
  dplyr::left_join(count.obs.two2, by = c("chrom1", "chrom2", "start", "end")) %>%
  dplyr::left_join(count.obs.two, by=c("chrom1", "chrom2", "start", "end", "alt1", "alt2")) %>%
  dplyr::filter(n <= frequency.filter) %>% #filter
  dplyr::relocate(run, sample, genes.interest, chrom1,chrom2, start, end, .before = EVENT)

#Store how many variants remain
filters.variants <- filters.variants %>%
  dplyr::mutate(nsamples=c(nrow(overlaps.one.sample.count), nrow(overlaps.two.sample.count)))


## B)  Genomic ranges, variants with the breakpoint inside the ROI-----

#One-break
# Convert `overlaps.one.sample.count` to a GRanges object
all.one.break.GR <- regioneR::toGRanges(overlaps.one.sample.count)
# Find overlaps between `bedFile` and `all.one.break.GR`
# This identifies which rows from `bedFile` and `all.one.break.GR` overlap
b <- GenomicRanges::findOverlaps(bedFile, all.one.break.GR)
# Compute the exact intersection regions between overlapping entries
overlaps.gr <- GenomicRanges::pintersect(all.one.break.GR[subjectHits(b)],bedFile[queryHits(b)])
# Extract gene names from `bedFile` metadata and assign them to the overlapping regions
overlaps.gr$gene <- elementMetadata(bedFile[queryHits(b)])$name
# Convert GRanges object to a data frame for easier manipulation
overlaps.one <- overlaps.gr %>% as.data.frame(row.names=NULL) %>%
  dplyr::rename("chr"="seqnames") %>%
  select(-start, -end, -width, -strand) %>% #we make this step because if not the variants are splited according to the bed file
  unique()


# Two-breaks
# We try to see if any of the two breakpoints (POS.x or POS.y) is inside the region of interest (ROI)
# Convert `all.two.break` into two separate GRanges objects for each breakpoint (x and y)
# Convert the first breakpoint (CHROM.x, POS.x) into a GRanges object
all.two.break.GR.x <-  all.two.break %>%
  dplyr::mutate(POS.x2 = POS.x) %>%
  dplyr::select(CHROM.x, POS.x, POS.x2, EVENT, MATEID.x, MATEID.y, sample) %>%
  regioneR::toGRanges()
# Convert the second breakpoint (CHROM.y, POS.y) into a GRanges object
all.two.break.GR.y <-  all.two.break %>%
  dplyr::mutate(POS.y2 = POS.y) %>%
  dplyr::select(CHROM.y, POS.y, POS.y2, EVENT, MATEID.x, MATEID.y, sample) %>%
  regioneR::toGRanges()

# Find overlaps between the breakpoints and the ROI (bedFile)
two.b.x <- GenomicRanges::findOverlaps(bedFile, all.two.break.GR.x)
two.b.y <- GenomicRanges::findOverlaps(bedFile, all.two.break.GR.y)
# Get the exact overlapping regions
overlaps.gr.two.x <- GenomicRanges::pintersect(all.two.break.GR.x[subjectHits(two.b.x)],bedFile[queryHits(two.b.x)])
overlaps.gr.two.y <- GenomicRanges::pintersect(all.two.break.GR.y[subjectHits(two.b.y)],bedFile[queryHits(two.b.y)])

# Extract gene names from `bedFile` and assign them to the overlapping regions
overlaps.gr.two.x$gene.x <- elementMetadata(bedFile[queryHits(two.b.x)])$name
overlaps.gr.two.y$gene.y <- elementMetadata(bedFile[queryHits(two.b.y)])$name

# Convert GRanges to data frame and keep only unique entries
overlaps.two.x <- overlaps.gr.two.x %>%
  as.data.frame(row.names = 1:length(overlaps.gr.two.x)) %>%
  unique()
overlaps.two.y <- overlaps.gr.two.y %>%
  as.data.frame(row.names = 1:length(overlaps.gr.two.y)) %>%
  unique()

# Merge results from both breakpoints (x and y) to keep both gene associations
# We merge by EVENT, MATEID.x, MATEID.y, and sample to maintain variant integrity
all.events <- merge(overlaps.two.x, overlaps.two.y, by=c("EVENT", "MATEID.x", "MATEID.y", "sample"), all.x=TRUE, all.y=TRUE) %>%
  dplyr::select(EVENT, MATEID.x, MATEID.y, sample, gene.x, gene.y) %>%
  unique()

# Merge the event information with additional variant data (`overlaps.two.sam
overlaps.two<- merge(all.events,  overlaps.two.sample.count , by=c(c("EVENT", "MATEID.x", "MATEID.y", "sample")))



# B.2) Only consider genes of interest (if you want to analyze all genes, put ALL_BEDFILE in the colunm genes.interest.

# Function to extract gene names from a column containing "ENSTxxx|GENE"
extract_gene_name <- function(gene_column) {
  stringr::str_extract(gene_column, "(?<=\\|)[^|]+")  # Agafa el text després de "|"
}
# Function to filter events based on genes of interest  # Extract gene names from `gene.x` and `gene.y` columns
filter_genes_of_interest <- function(df) {
  if ("gene.x" %in% colnames(df) & "gene.y" %in% colnames(df)) {
    # Two-breakend case: Extract gene names from `gene.x` and `gene.y`
    df <- df %>%
      mutate(
        gene.x.name = extract_gene_name(gene.x),
        gene.y.name = extract_gene_name(gene.y)
      )
  } else if ("gene" %in% colnames(df)) {
    # One-breakend case: Extract gene name from `gene`
    df <- df %>%
      mutate(gene.name = extract_gene_name(gene))
  } else {
    stop("Error: No gene columns found (expected 'gene.x' and 'gene.y' for two breakends or 'gene' for one breakend).")
  }

  # Filter data independently for each sample
  df.filtered <- df %>%
    group_by(sample) %>%  # Group by sample to apply different filters per sample
    filter(
      unique(genes.interest) == "ALL_BEDFILE" |
        any(across(any_of(c("gene.x.name", "gene.y.name", "gene.name")),
                                                            ~ . %in% unlist(strsplit(unique(genes.interest), ",\\s*"))
      ))
    ) %>%
    ungroup()  # Remove grouping
  return(df.filtered)  # Return the filtered dataset
}

ov.one.cascade <- filter_genes_of_interest(overlaps.one)
ov.two.cascade <- filter_genes_of_interest(overlaps.two)

filters.variants <- filters.variants %>%
  dplyr::mutate(phenotype=c(nrow(ov.one.cascade), nrow(ov.two.cascade)))

#C) Filter by gt_AF ----
AF.one <- ov.one.cascade %>%
  dplyr::filter(gt_AF >= args$gt_AF)

AF.two <- ov.two.cascade %>%
  dplyr::filter(gt_AF.x >= args$gt_AF &
                  gt_AF.y >= args$gt_AF)

filters.variants <- filters.variants %>%
  dplyr::mutate(gt_AF=c(nrow(AF.one), nrow(AF.two)))


#D) Delete variants in repetitive regions
ov.one.no.repeats <- AF.one %>%
  dplyr::filter(!(INSRMRC %in% c("Low_complexity", "Simple_repeat") ))
ov.two.no.repeats <- AF.two %>%
  dplyr::filter(!(INSRMRC.x %in% c("Low_complexity", "Simple_repeat") & INSRMRC.y %in% c("Low_complexity", "Simple_repeat") ))

filters.variants <- filters.variants %>%
  dplyr::mutate(simple=c(nrow(ov.one.no.repeats), nrow(ov.two.no.repeats)))

#E) Delete highly similar variants (CHR, POS)
ov.one.hot.spot <- ov.one.no.repeats %>%
  dplyr::filter(n2<=highly.similar)

ov.two.hot.spot <- ov.two.no.repeats %>%
  dplyr::filter(n2<=highly.similar)

filters.variants <- filters.variants %>%
  dplyr::mutate(similar=c(nrow(ov.one.hot.spot), nrow(ov.two.hot.spot)))


# F) For one break get only MEI variants
one.b.definitive <- ov.one.hot.spot %>%
  dplyr::filter(INSRMRC %in% c("SINE", "LINE"))

#E) For two breakends only consider variants that are bigger than X (if located at the same chr)

count_ATCG <- function(seq_vector) {
  extracted <- str_extract_all(seq_vector, "[ATCG]+")  # Extreu només A, T, C, G
  lengths <- sapply(extracted, nchar)  # Compta quants nucleòtids hi ha en cada entrada
  return(lengths)
}
two.b.definitive <- ov.two.hot.spot %>%
  dplyr::filter(CHROM.x!=CHROM.y | CHROM.x==CHROM.y & (end-start >= minimumLength | (end-start ==1 & count_ATCG(ALT.x) >= minimumLegth)))



#Plot diagram
workflow_filtering <- function(what, filters, one.def, two.def){
wf<- grViz(paste("
digraph flow {
  graph [layout = dot, rankdir = TB]

  node [shape = rectangle, style = filled, fillcolor = white, fontname = Helvetica, fontsize = 12]

  # Nodes
  Total_Variants [label = 'Total", what, "breakend", filters.variants$totals," variants', shape = box, style = filled, fillcolor = lightgrey]

  Filter1 [label = 'Present in \n ≤", args$frequency, "samples?', shape = diamond, fillcolor = white]
  Pass1 [label ='",  filters.variants$nsample,"variants', fillcolor = lightblue]
  Fail1 [label = '", filters.variants$totals-filters.variants$totals, "variants', fillcolor = lightcoral]

  Filter2 [label = 'Phenotype genes?', shape = diamond, fillcolor = white]
  Pass2 [label = '",  filters.variants$phenotype,"variants', fillcolor = lightblue]
  Fail2 [label = '",  filters.variants$nsample-filters.variants$phenotype," variants', fillcolor = lightcoral]

  Filter3 [label = 'VAF ≥", args$gt_AF, "', shape = diamond, fillcolor = white]
  Pass3 [label = '",  filters.variants$gt_AF,"variants', fillcolor = lightblue]
  Fail3 [label = '",  filters.variants$phenotype- filters.variants$gt_AF,"variants', fillcolor = lightcoral]

  Filter4 [label = 'outside a simple \n repeat  or low \n complexity region', shape = diamond, fillcolor = white]
  Pass4 [label = '",  filters.variants$simple,"variants', fillcolor = lightblue]
  Fail4 [label = '",  filters.variants$gt_AF- filters.variants$simple,"variants', fillcolor = lightcoral]

  Filter5 [label = 'highly similar \n variant ≤", args$similar, " ', shape = diamond, fillcolor = white]
  Pass5 [label = '",  filters.variants$similar, "variants', fillcolor = lightblue]
  Fail5 [label = '",  filters.variants$simple- filters.variants$similar,"variants', fillcolor = lightcoral]",

  if(what=="one"){
    paste("Filter6 [label = 'Transposable \n elements', shape = diamond, fillcolor = white]
    Pass6 [label = '",  nrow(one.b.definitive),"variants', fillcolor = lightblue]
    Fail6 [label = '",  filters.variants$similar- nrow(one.b.definitive),"variants', fillcolor = lightcoral]")
  }else{
    paste("Filter6 [label = 'minimum length \n ≥", args$minimumLength, "bp', shape = diamond, fillcolor = white]
    Pass6 [label = '",  nrow(two.b.definitive),"variants', fillcolor = lightblue]
    Fail6 [label = '",  filters.variants$similar- nrow(two.b.definitive),"variants', fillcolor = lightcoral]")
  },
  "
  # Edges
  Total_Variants -> Filter1
  Filter1 -> Pass1 [label = 'Yes']
  Filter1 -> Fail1 [label = 'No']

  Pass1 -> Filter2
  Filter2 -> Pass2 [label = 'Yes']
  Filter2 -> Fail2 [label = 'No']

  Pass2 -> Filter3
  Filter3 -> Pass3 [label = 'Yes']
  Filter3 -> Fail3 [label = 'No']

  Pass3 -> Filter4
  Filter4 -> Pass4 [label = 'Yes']
  Filter4 -> Fail4 [label = 'No']

  Pass4 -> Filter5
  Filter5 -> Pass5 [label = 'Yes']
  Filter5 -> Fail5 [label = 'No']

  Pass5 -> Filter6
  Filter6 -> Pass6 [label = 'Yes']
  Filter6 -> Fail6 [label = 'No']
}"))
return(wf)
}

plots.dir <- file.path(args$output, "plots")
dir.create(plots.dir)
one.wf <- workflow_filtering("one", filters.variants[1,], one.b.definitive, two.b.definitive)
two.wf <- workflow_filtering("two", filters.variants[2,], one.b.definitive, two.b.definitive)


# convert to SVG
svg.one <- DiagrammeRsvg::export_svg(one.wf)
svg.two <- DiagrammeRsvg::export_svg(two.wf)

# Save as PNG
rsvg_png(charToRaw(svg.one), file = file.path(plots.dir, "workflow_diagram_one.png"))
rsvg_png(charToRaw(svg.two), file = file.path(plots.dir, "workflow_diagram_two.png"))





