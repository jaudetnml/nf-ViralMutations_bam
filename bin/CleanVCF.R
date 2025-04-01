#!/usr/bin/env Rscript

library(tidyverse)

makeGT <- function(alt) {
  print(alt)
  n_alts <- str_count(alt, ",") + 1
  print(n_alts)
  if (n_alts > 1) {
    GTs <- 1:n_alts
    return(paste(GTs, collapse = "/"))
  } else {
    return("1/1")
  }
}

DepthFromCounts <- function(counts) {
  counts_split <- str_split(counts, ",")
  counts_numeric <- lapply(counts_split, as.numeric)
  depths <- sapply(counts_numeric, sum)
  return(depths)
}

args <- commandArgs(trailingOnly = TRUE)

filename <- args[[1]]
MinFreq <- as.numeric(args[[2]])
MinDepth <- as.numeric(args[[3]])

raw_vcf <- read_tsv(filename,
  comment = "#",
  col_names = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", NA, NA, "Dat"),
  col_types = "cicccnc--c"
)

raw_singles <- filter(raw_vcf, !str_detect(ALT, ","))
raw_mult <- filter(raw_vcf, str_detect(ALT, ","))

calc_singles <- raw_singles %>%
  separate_wider_delim(Dat, delim = ":", names = c(NA, "Depth", "Counts", NA, NA, NA, NA, NA)) %>%
  separate_wider_delim(Counts, delim = ",", names = c(NA, "Alt_Count")) %>%
  mutate(
    Depth = as.numeric(Depth),
    Alt_Count = as.numeric(Alt_Count),
    Alt_Freq = signif(Alt_Count / Depth, round(log10(Depth), 0))
  )

if (nrow(raw_mult) > 0) {
  split_calc_mult <- raw_mult %>%
    separate_wider_delim(Dat, delim = ":", names = c(NA, NA, "Counts", NA, NA, NA, NA, NA)) %>%
    mutate(Depth = DepthFromCounts(Counts)) %>%
    mutate(Counts = str_remove(Counts, "^[0-9]+,")) %>%
    separate_longer_delim(c(ALT, Counts), delim = ",") %>%
    mutate(
      Depth = as.numeric(Depth),
      Counts = as.numeric(Counts),
      Alt_Freq = signif(Counts / Depth, round(log10(Depth), 0))
    )
} else {
  split_calc_mult <- tibble(
    `#CHROM` = character(),
    POS = numeric(),
    REF = character(),
    ID = character(),
    ALT = character(),
    Depth = numeric(),
    Counts = numeric(),
    Alt_Freq = numeric()
  )
}

bind_rows(
  calc_singles,
  rename(split_calc_mult, Alt_Count = Counts)
) %>%
  # mutate(`#CHROM` = str_remove_all(`#CHROM`, "\\.[0-9]+")) %>%
  rowwise() %>%
  mutate(INFO = paste0("AF=", Alt_Freq, ";AD=", Alt_Count, ";DP=", Depth, collapse = "")) %>%
  ungroup() %>%
  select(-Depth, -Alt_Count, -Alt_Freq) %>%
  arrange(`#CHROM`, POS) %>%
  readr::write_tsv(str_replace(filename, "_variants.vcf", "_clean.vcf"))

mult_minority_consensus <- split_calc_mult %>%
  filter(Alt_Freq > 0.1, nchar(REF) == nchar(ALT)) %>% # Ignore alleles with less than 10% frequency, ignore indels (freebayes reports the depth weirdly)
  group_by(POS, REF) %>%
  filter(n() > 1, all(Alt_Freq < MinFreq), sum(Alt_Freq) >= (1 - MinFreq), all(Depth >= MinDepth)) %>% # Keep positions with multiple alts and with no consensus alt and the reference is also not consensus
  ungroup() %>%
  group_by(`#CHROM`, REF, ID, POS) %>%
  summarise(ALT = paste(ALT, collapse = ",")) %>%
  ungroup() %>%
  select(`#CHROM`, POS, ID, REF, ALT)

mult_majority_consensus <- split_calc_mult %>%
  mutate(UpperCI = qbeta(0.975, round(Depth * Alt_Freq), round(Depth * (1 - Alt_Freq)))) %>%
  filter((Alt_Freq >= MinFreq) | ((Alt_Freq > 0.5) & (UpperCI > MinFreq)), Depth >= MinDepth) %>%
  ungroup() %>%
  select(`#CHROM`, POS, ID, REF, ALT)

consensus_pre <- calc_singles %>%
  mutate(UpperCI = qbeta(0.975, round(Depth * Alt_Freq), round(Depth * (1 - Alt_Freq)))) %>%
  filter((Alt_Freq >= MinFreq) | ((Alt_Freq > 0.5) & (UpperCI > MinFreq)), Depth >= MinDepth) %>%
  select(`#CHROM`, POS, ID, REF, ALT) %>%
  bind_rows(mult_minority_consensus, mult_majority_consensus)

print(consensus_pre, n = Inf)

if (nrow(consensus_pre) >= 1) {
  consensus_pre %>%
    arrange(`#CHROM`, POS) %>%
    mutate(
      QUAL = NA,
      FILTER = NA,
      INFO = NA,
      FORMAT = "GT"
    ) %>%
    rowwise() %>%
    mutate(unknown = makeGT(ALT)) %>%
    ungroup() %>%
    readr::write_tsv(str_replace(filename, "_variants.vcf", "_clean_consensus.vcf"))
}
