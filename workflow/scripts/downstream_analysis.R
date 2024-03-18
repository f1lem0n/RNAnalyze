# stdout/stderr logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log)

# libraries
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
})

# # read the data
tblastn_out <- read.table(
  file = snakemake@input[[1]],
  sep = "\t",
)
colnames(tblastn_out) <- c(
  "query",
  "subject",
  "identity",
  "length",
  "mismatches",
  "gaps",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore"
)

### Downstream analysis ###

unique_queries <- unique(tblastn_out$query)
unique_subjects <- unique(tblastn_out$subject)
binary_mtx <- outer(unique_queries, unique_subjects, Vectorize(function(x, y) {
  ifelse(any(tblastn_out$query == x & tblastn_out$subject == y), 1, 0)
}))
binary_df <- data.frame(binary_mtx, row.names = unique_queries)
colnames(binary_df) <- unique_subjects
binary_melt <- melt(as.matrix(binary_df))
colnames(binary_melt) <- c("query", "subject", "hit")
binary_melt$hit <- as.logical(binary_melt$hit)
qs_heatmap <- ggplot(
  binary_melt,
  aes(x = subject, y = query, fill = hit)
) +
  geom_tile() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  )

ggsave(
  filename = paste0(snakemake@output$dir, "/qs_mosaic.png"),
  plot = qs_heatmap,
  width = 10,
  height = 10
)

# create .done file
file.create(snakemake@output[[1]])
