# Read the frozen RNAfold results and measured GC content in data-raw/ and
# join them into one row per synthetic sequence.
#
# Each batch letter (a-h) covers a band of lengths. Sequence names encode the
# design: len_<nt>_GC_<target percent>_iter_<replicate>.

parse_sequence_name <- function(x) {
  tibble(
    length    = as.integer(str_extract(x, "(?<=len_)\\d+")),
    target_GC = as.integer(str_extract(x, "(?<=GC_)\\d+")),
    iteration = as.integer(str_extract(x, "(?<=iter_)\\d+"))
  )
}

read_rnafold_file <- function(path) {
  raw <- read_csv(path, col_types = cols(.default = "c"))
  bind_cols(parse_sequence_name(raw$sequence_details),
            MFE = as.numeric(raw$rnafold))
}

read_gc_file <- function(path) {
  raw <- read_tsv(path, col_names = c("sequence_details", "GC_content"),
                  col_types = "cd")
  bind_cols(parse_sequence_name(raw$sequence_details),
            GC_content = raw$GC_content)
}

load_mfe_data <- function(dir = "data-raw") {
  mfe_files <- list.files(dir, "^[a-z]_rnafold_results\\.csv$", full.names = TRUE)
  gc_files <- list.files(dir, "^[a-z]_actualGC\\.txt$", full.names = TRUE)
  if (length(mfe_files) == 0 || length(mfe_files) != length(gc_files)) {
    stop("Expected matching *_rnafold_results.csv and *_actualGC.txt files in ", dir)
  }

  mfe <- map_dfr(mfe_files, read_rnafold_file)
  gc <- map_dfr(gc_files, read_gc_file)

  data <- inner_join(mfe, gc, by = c("length", "target_GC", "iteration"),
                     relationship = "one-to-one") |>
    arrange(length, target_GC, iteration)

  if (nrow(data) != nrow(mfe) || anyNA(data)) {
    stop("RNAfold and GC files do not match one-to-one (",
         nrow(mfe), " MFE rows, ", nrow(data), " joined, ",
         sum(!complete.cases(data)), " incomplete).")
  }
  data
}
