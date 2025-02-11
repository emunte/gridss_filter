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
packages <- c("dplyr", "optparse", "stringr", "vcfR", "GenomicRanges", "yaml", "furrr")

#Apply
invisible(sapply(packages, install_if_missing))


# 2. Build options list----
option_list <- list(
  optparse::make_option(c("-t", "--txt"), type="character", default="",
              help="Path to txt with vcf path and sample information. See README", metavar="character"),
  optparse::make_option(c("-c", "--cores"), type="integer", default=1,
                        help="Number of cores to use", metavar="character"),
  optparse::make_option(c("-o", "--output"), type="character", default="",
                        help="Path to the output folder. Results will be stored in this folder", metavar="character")
)

opt_parser <- optparse::OptionParser(option_list=option_list);


# 3. Load params----
args <- optparse::parse_args(opt_parser)
samples.txt <- args$txt
output <- args$output
cores <- args$cores


# 3. Parallelization----
start <- print(Sys.time()) #start time
future::plan(multisession, workers = cores)

# 4. Functions---

#Process the vcf file provided by Gridss repeat Masker
process_vcf <- function(i, samples.df){
  print(i)
  file.pass <- data.frame()
  vcf <- tryCatch({
    vcfR::read.vcfR(samples.df$path[i])
  }, error = function(e) return(NULL))  # if an error occurs continue

  if (!is.null(vcf) && nrow(vcf@fix) > 0) {
    vcf.tidy <- vcfR::vcfR2tidy(vcf, single_frame = TRUE) #convert into a dataframe
    colnames(vcf.tidy$dat) <- gsub("^REF...61", "REF_num", colnames(vcf.tidy$dat)) #avoid repeated names
    colnames(vcf.tidy$dat) <- gsub("^REF...4", "REF", colnames(vcf.tidy$dat))
    #colnames(vcf.tidy$dat) <- make.unique(colnames(vcf.tidy$dat))
    
    file.et <- as.data.frame(vcf.tidy$dat) %>%
      mutate(sample = samples.df$sample[i],
             run = samples.df$run[i],
             genes.interest = samples.df$genes.interest[i])

    if (nrow(file.et) > 0) {
      file.et.h <- file.et %>% filter(str_detect(ID, "[0-9]h"))
      file.et.o <- file.et %>% filter(str_detect(ID, "[0-9]o"))

      files.et.all <- merge(file.et.o, file.et.h, by = c("EVENT", "sample", "run", "genes.interest")) %>% #merge two breakend variants
        filter(FILTER.x == "PASS" | FILTER.y == "PASS")

      file.et.b <- file.et %>%
        filter(str_detect(ID, "[0-9]b$"), FILTER == "PASS")  #one breakend variants
      gc()
      return(list(two_breaks = files.et.all, one_breaks = file.et.b))
    }
  }
  gc()
  return(NULL)

}


# 5 Run the pipeline----

#read the txt file
samples.df <- read.delim(samples.txt) 

#Execute the function
results <- furrr::future_map(1:nrow(samples.df), process_vcf, samples.df, .progress = TRUE, .options = furrr_options(seed = TRUE))

# Stop parallelization
future::plan(sequential)  
#Save memory
gc() 

#bind the results
files.pass.all <- bind_rows(lapply(results, `[[`, "two_breaks"))
files.et.b.all <- bind_rows(lapply(results, `[[`, "one_breaks"))

#Create the otuput dir
merge.dir <- file.path(output, paste0(Sys.Date(), "_merge_gridss_vcfs"))
dir.create(merge.dir, showWarnings = FALSE)

#print the results in txt files
write.table(files.pass.all, file.path(merge.dir, paste0(Sys.Date(), "_merge_gridss_two_break.txt")), row.names = FALSE, sep = "\t")
write.table(files.et.b.all, file.path(merge.dir, paste0(Sys.Date(), "_merge_gridss_one_break.txt")), row.names = FALSE, sep = "\t")


print(paste("Executed in", Sys.time()- start, "seconds"))


