################################################################################
## Merge gridss results in two dataframes
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
packages <- c("dplyr", "optparse", "stringr", "vcfR", "GenomicRanges", "yaml")

#Apply
invisible(sapply(packages, install_if_missing))


# 2. Build options list
option_list <- list(
  optparse::make_option(c("-t", "--txt"), type="character", default="",
              help="Path to txt with vcf path and sample information. See README", metavar="character"),
  optparse::make_option(c("-o", "--output"), type="character", default="",
                        help="Path to the output folder. Results will be stored in this folder", metavar="character")
)

opt_parser <- optparse::OptionParser(option_list=option_list);
# Load params
args <- optparse::parse_args(opt_parser)
samples.txt <- yaml::yaml.load(args$txt)
output <- yaml::yaml.load(args$output)

# 3. Pipeline

  samples.df <- read.delim(samples.txt) #load the txt file
  #create dataframes
  files.pass.all <- data.frame()
  files.et.b.all <- data.frame()
  for (i in 1:nrow(samples.df)){
    print(i)
    #create a data.frame to store all the variants
    file.pass <- data.frame()
    #read the vcf file
    vcf <-vcfR::read.vcfR(samples.df$path[i])
    if(vcf@fix %>% nrow()>0){
      #convert vcf file to dataframe
      vcf.tidy <- vcfR::vcfR2tidy(vcf, single_frame = TRUE)
      file.et <- as.data.frame(vcf.tidy$dat) %>%
        dplyr::mutate(sample = samples.df$sample[i],
                      run = samples.df$run[i],
                      genes.interest= samples.df$genes.interest[i])

      #merge variants with the samme mateID and filter only PASS variants
      if(nrow(file.et)>0){
        file.et.h <- file.et %>% dplyr::filter(stringr::str_detect(ID, "[0-9]h"))
        file.et.o <- file.et %>% dplyr::filter(stringr::str_detect(ID, "[0-9]o"))
        files.et.all <- merge(file.et.o, file.et.h, by="EVENT") %>%
          dplyr::filter(FILTER.x=="PASS"|FILTER.y=="PASS")

      #variants with b are single-breakend (we will store them in a dataframe)
        file.et.b <- file.et %>% dplyr::filter(stringr::str_detect(ID, "[0-9]b$"), FILTER=="PASS")
        files.pass.all <- rbind(files.pass.all, files.et.all)
        files.et.b.all <- rbind(files.et.b.all, file.et.b)
      }
    }
  }
  #write results
  merge.dir <- file.path(output, "merge_gridss_vcfs")
  dir.create(merge.dir)
  write.table(files.pass.all, file.path(merge.dir, paste0(Sys.Date(), "merge_gridss_two_break.txt")), row.names=FALSE)
  write.txt(files.et.b.all , file.path(merge.dir, paste0(Sys.Date(), "merge_gridss_one_break.csv")), row.names=FALSE)






