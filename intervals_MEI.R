
intervals_INSRMRT <- function(df, merge.vars){
# 1) Prepare the 'df' data frame:
df <- all.one.break %>%
  # a) Filter out rows where INSRMRC is missing
  #dplyr::filter(!is.na(INSRMRC)) %>%

  # b) Sort by CHROM and POS
  #arrange(CHROM, POS) %>%

  # c) Group by CHROM and INSRMRT
  dplyr::group_by(CHROM, INSRMRC) %>%

  # d) Create intervals around POS (±4)
  mutate(
    lower = POS - merge.vars ,
    upper = POS + merge.vars
  ) %>%

  # e) Remove grouping
  ungroup()




# 2) Convert the data frame 'df' to a GRanges object

#   - seqnames = df$CHROM
#   - ranges = IRanges(start = df$lower, end = df$upper)
#   - transposon = df$INSRMRT
#   - row.id = sequential number for each row
gr <- GRanges(
  seqnames   = df$CHROM,
  ranges     = IRanges(start = df$lower, end = df$upper),
  transposon = df$INSRMRC,
  row.id     = seq_len(nrow(df))
)


# 3) Create a grouping label for each combination of (CHROM, INSRMRT)
df$group <- paste(df$CHROM, df$INSRMRC, sep = "_")

#    This splits 'gr' into a list (GRangesList), one entry per (CHROM, INSRMRT)
gr.list  <- split(gr, df$group)
# class(gr.list) => "GRangesList"


# 4) For each sub-GRanges, reduce overlapping intervals and assign cluster.id
res.list <- lapply(gr.list, function(gr.sub) {

  # a) Merge overlapping intervals
  red  <- reduce(gr.sub)

  # b) Determine which merged interval each original interval belongs to
  hits <- findOverlaps(gr.sub, red)

  # c) subjectHits(hits) gives the index of the "reduced" interval in 'red'
  #    This becomes the cluster ID for each interval.
  gr.sub$cluster.id <- subjectHits(hits)

  # d) Return the modified GRanges with a new column 'cluster.id'
  gr.sub
})


# 5) Combine the list of GRanges ('res.list') into a single GRanges object
#    The 'do.call(c, ...)' approach concatenates all GRanges into one
gr.combined <- do.call(c, unname(res.list))
# class(gr.combined) => "GRanges"


# 6) Convert the unified GRanges 'gr.combined' back to a data frame
df.clusters <- data.frame(
  CHROM       = as.character(seqnames(gr.combined)),  # chromosome
  start       = start(gr.combined),                   # interval start
  end         = end(gr.combined),                     # interval end
  cluster.id  = gr.combined$cluster.id,               # cluster ID
  transposon  = gr.combined$transposon,               # type of transposon
  row.id      = gr.combined$row.id                    # row ID to join back
)


# 7) Summarize: count how many variants in each cluster
df.summary <- df.clusters %>%
  group_by(CHROM, transposon, cluster.id) %>%
  summarise(
    n.variants = n(),  # number of variants in this cluster
    .groups = "drop"
  )


# 8) Merge 'n.variants' back into df.clusters, so each row knows how many
#    variants share its cluster.id
df.clusters <- df.clusters %>%
  left_join(
    df.summary,
    by = c("CHROM", "transposon", "cluster.id")
  )


# 9) Finally, join cluster information back to the original data frame 'df'
df.enriquit <- df %>%
  left_join(
    df.clusters %>%
      select(row.id, cluster.id, n.variants),
    by = "row.id"
  )

return(df.enriquit)
}

# 'df.enriquit' now contains all original columns plus:
#   - cluster.id: which cluster each row belongs to
#   - n.variants: how many variants are in that same cluster

