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
packages <- c("dplyr", "optparse", "stringr", "GenomicRanges", "yaml")

#Apply
invisible(sapply(packages, install_if_missing))



# 2. Build options list

option_list <- list(
  optparse::make_option(c("-i", "--input"), type="character", default="",
                        help="Folder with merge results", metavar="character"),
  optparse::make_option(c("-d", "--bedFile"), type="character", default="",
                       help="path to the bed file", metavar="character"),
  optparse::make_option(c("-f", "--frequency"), type="character", default=10,
                        help="Maximum allowed occurrences before a variant is excluded", metavar="integer"),
  optparse::make_option(c("-o", "--output"), type="character", default="",
                        help="Path to the output directory. Results will be stored in this folder", metavar="character")
)

opt_parser <- optparse::OptionParser(option_list=option_list);
# Load params
args <- optparse::parse_args(opt_parser)
input <- yaml::yaml.load(args$txt)
output <- yaml::yaml.load(args$output)
bedFile <- yaml::yaml.load(args$bedFile) %>% rtracklayer::import()
frequency.filter <- yaml::yaml.load(args$filter.frequency)


# 3. Load data

#Get the two dataframes (merge_gridss_vcf output)
all.one.break <- list.files(input, pattern="merge_gridss_one_break.txt", recursive=TRUE, full.names=TRUE) %>%
  read.delim(header=TRUE, sep="\t") %>%
  dplyr::mutate(end=POS)


all.two.break <- list.files(input, pattern="merge_gridss_two_break.txt", recursive=TRUE, full.names=TRUE) %>%
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
  dplyr::relocate(run, sample, genes.interest, CHROM, start, end, .before = POS)

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


#Two-breaks





