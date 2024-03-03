# stdout/stderr logging
log <- file(snakemake@log[[1]], open = "wt")
sink(log)

# libraries
suppressMessages({
  library(dplyr)
})

# the rest of the script
file.create(snakemake@output[[1]])
