# =============================================================================
# stacks_summary_to_mantel.R
#
# Parses a STACKS pairwise summary file, extracts the relevant "Means" block,
# builds a full symmetric distance matrix, and runs Mantel tests via ade4 and
# vegan.
#
# Two STACKS file types are supported (select via `dataset` below):
#   • populations.phistats_summary.tsv  ->  Fst' Means block (standardised Fst')
#   • populations.fst_summary.tsv       ->  Fst  Means block (Weir & Cockerham)
#
# Both share the same upper-triangle "Means" matrix layout, so the same
# parser handles either one.
#
# Handles:
#   • Multi-block STACKS summary file
#   • Upper-triangle-only matrix format
#   • Missing last-population data row (the SpringCreek problem)
#   • Windows line endings and UTF-8 BOM
#   • Prime vs apostrophe character in "Fst'" section header
#
# Requirements : ade4, vegan
# Optional     : geosphere (for Haversine geographic distances in Section 6)
# =============================================================================

library(ade4)
library(vegan)
rm(list=ls())

# =============================================================================
# 0. Choose which STACKS pairwise summary file to analyse
#
#   Both supported file types share the same upper-triangle "Means" matrix
#   layout, so the same parser is used for both. Set `dataset` to pick one:
#
#     "phistats" : populations.phistats_summary.tsv  ->  Fst' Means block
#                  (Hedrick's standardised Fst' — the original behaviour)
#     "fst"      : populations.fst_summary.tsv        ->  Fst  Means block
#                  (Weir & Cockerham Fst)
#
#   The phistats summary file contains several blocks (phi_st, Fst', D_est);
#   `pattern` selects the Fst' one. The fst summary file contains a single
#   Fst Means block. Both match the regex "Fst.* Means" (the ".*" also absorbs
#   the apostrophe/prime vs "Fst." spelling differences across encodings).
# =============================================================================
setwd("~/snails/L_virgata/population-genomics/virgata_m3M3n3_R80_maf025/")
dataset <- "fst"          # <-- "phistats" or "fst"

dataset_config <- list(
  phistats = list(
    filepath   = "populations.phistats_summary.tsv",
    pattern    = "Fst.* Means",
    stat_label = "Fst'",       # human-readable label for messages/plot titles
    out_prefix = "fst_prime",  # prefix for exported files
    use_prime  = TRUE          # whether to draw the prime symbol on plot axes
  ),
  fst = list(
    filepath   = "populations.fst_summary.tsv",
    pattern    = "Fst.* Means",
    stat_label = "Fst",
    out_prefix = "fst",
    use_prime  = FALSE
  )
)

if (!dataset %in% names(dataset_config))
  stop(sprintf("Unknown dataset '%s'. Choose one of: %s",
               dataset, paste(names(dataset_config), collapse = ", ")))

cfg        <- dataset_config[[dataset]]
stat_label <- cfg$stat_label
out_prefix <- cfg$out_prefix

# Plot-axis expression for the linearised ratio, with or without the prime
stat_ratio_expr <- if (cfg$use_prime) {
  expression(F[ST]*"'" / (1 - F[ST]*"'"))
} else {
  expression(F[ST] / (1 - F[ST]))
}

cat(sprintf("Analysing dataset '%s'  ->  file: %s  |  statistic: %s\n\n",
            dataset, cfg$filepath, stat_label))

# =============================================================================
# 1. Read the raw file
# =============================================================================

filepath  <- cfg$filepath   # set in Section 0 (adjust paths there if needed)

raw_lines <- readLines(filepath, encoding = "UTF-8", warn = FALSE)
raw_lines <- gsub("\r",    "", raw_lines)   # strip Windows CR
raw_lines <- gsub("^\uFEFF", "", raw_lines) # strip UTF-8 BOM if present

# Report every section header found
section_idx <- grep("^#", raw_lines)
if (length(section_idx) == 0) stop("No section headers found — check file format.")

cat("Sections found in file:\n")
cat(paste0("  [", seq_along(section_idx), "]  line ",
           formatC(section_idx, width = 3), ":  ",
           raw_lines[section_idx]),
    sep = "\n")
cat("\n")

# =============================================================================
# 2. Block extraction function
#
#   Arguments
#     lines   : character vector of raw file lines
#     pattern : PERL-compatible regex matching the target section header
#
#   Returns a data.frame with:
#     • rows    = all populations (column header is the authoritative list)
#     • columns = same populations in the same order
#     • values  = upper-triangle Fst values; lower triangle is NA before
#                 symmetrisation (handled in Section 4)
#
#   SpringCreek fix:
#     STACKS does not write a data row for the last population because it has
#     no upper-triangle values. The function detects any populations present
#     in the column header but absent from the parsed rows and injects them
#     as all-NA rows so the matrix is always n×n.
# =============================================================================

extract_block <- function(lines, pattern) {
  
  # -- Locate section header --------------------------------------------------
  header_idx <- grep(pattern, lines, perl = TRUE)
  
  if (length(header_idx) == 0)
    stop(paste0("No section header matched pattern: '", pattern, "'"))
  
  if (length(header_idx) > 1) {
    warning(paste0("Pattern '", pattern, "' matched ", length(header_idx),
                   " lines — using first match at line ", header_idx[1], "."))
    header_idx <- header_idx[1]
  }
  
  # -- Column-name row (line immediately after section header) ----------------
  col_row_idx <- header_idx + 1
  if (col_row_idx > length(lines))
    stop(paste0("Section header at line ", header_idx,
                " is the last line in the file — no data follows."))
  
  # -- Determine end of block -------------------------------------------------
  all_hdr  <- grep("^#", lines)
  next_hdr <- all_hdr[all_hdr > header_idx]
  end_idx  <- if (length(next_hdr) == 0) length(lines) else next_hdr[1] - 1
  
  # -- Parse population names from column header (AUTHORITATIVE) --------------
  col_fields <- strsplit(lines[col_row_idx], "\t")[[1]]
  pop_names  <- col_fields[-1]                          # drop leading empty cell
  pop_names  <- pop_names[nzchar(trimws(pop_names))]    # drop trailing empties
  n          <- length(pop_names)
  
  if (n == 0)
    stop(paste0("Could not parse any population names from the column header ",
                "at line ", col_row_idx, ". Check delimiter (must be tab)."))
  
  # -- Collect non-blank data rows --------------------------------------------
  data_lines <- lines[(col_row_idx + 1):end_idx]
  data_lines <- data_lines[nzchar(trimws(data_lines))]
  
  if (length(data_lines) == 0)
    stop(paste0("No data rows found in block matched by pattern: '",
                pattern, "'"))
  
  # -- Parse each data row ----------------------------------------------------
  # Each row:  col[1] = population name; col[2..] = numeric values
  # Upper-triangle rows are shorter than n values — pad to exactly n with NA
  parse_row <- function(line) {
    fields   <- strsplit(line, "\t")[[1]]
    row_name <- fields[1]
    vals     <- suppressWarnings(as.numeric(fields[-1]))
    length(vals) <- n    # pad short rows; truncates over-long rows (shouldn't happen)
    list(name = row_name, vals = vals)
  }
  
  parsed    <- lapply(data_lines, parse_row)
  row_names <- vapply(parsed, `[[`, character(1), "name")
  mat_data  <- do.call(rbind, lapply(parsed, `[[`, "vals"))
  
  df           <- as.data.frame(mat_data, stringsAsFactors = FALSE)
  colnames(df) <- pop_names
  rownames(df) <- row_names
  df[]         <- lapply(df, function(x) as.numeric(as.character(x)))
  
  # -- SpringCreek fix --------------------------------------------------------
  # Populations in the column header but absent from parsed row names had no
  # upper-triangle values and therefore no row in the file. Inject NA rows.
  missing_pops <- setdiff(pop_names, row_names)
  
  if (length(missing_pops) > 0) {
    cat(sprintf("NOTE [pattern = '%s']:\n", pattern))
    cat("  The following population(s) had no data row in the file.\n")
    cat("  This is expected for the last population in an upper-triangle\n")
    cat("  matrix. Injecting all-NA row(s) to maintain an", n, "x", n, "matrix:\n")
    cat(paste0("    • ", missing_pops), sep = "\n")
    cat("\n")
    
    na_rows <- as.data.frame(
      matrix(NA_real_,
             nrow     = length(missing_pops),
             ncol     = n,
             dimnames = list(missing_pops, pop_names))
    )
    df <- rbind(df, na_rows)
  }
  
  # -- Enforce canonical row and column order ---------------------------------
  # Reorders both dimensions to match the column-header population sequence
  df <- df[pop_names, pop_names]
  
  return(df)
}

# =============================================================================
# 3. Extract the relevant "Means" block (Fst' or Fst, per `dataset`)
#
#    The pattern is taken from the dataset config (Section 0). "Fst.* Means"
#    matches regardless of whether the file uses an apostrophe (U+0027), a
#    prime character (U+2032), or a plain "Fst." after "Fst", avoiding silent
#    failures due to encoding differences.
# =============================================================================

fst_raw   <- extract_block(raw_lines, cfg$pattern)
pop_names <- rownames(fst_raw)
n         <- length(pop_names)

cat("Populations in", stat_label, "block (n =", n, "):\n")
cat(paste0("  ", formatC(seq_len(n), width = 2), ".  ", pop_names), sep = "\n")
cat("\n")

# =============================================================================
# 4. Build full symmetric n×n matrix from upper triangle
#
#    Vectorised three-step approach:
#      (a) Replace NA with 0 — lower triangle and diagonal were absent/NA
#      (b) Add the matrix to its own transpose — mirrors upper → lower;
#          diagonal remains 0 + 0 = 0
#      (c) Restore dimnames (preserved through the operation but set explicitly)
# =============================================================================

fst_work              <- as.matrix(fst_raw)
fst_work[is.na(fst_work)] <- 0
fst_mat               <- fst_work + t(fst_work)
dimnames(fst_mat)     <- list(pop_names, pop_names)

# ---- Sanity checks ----------------------------------------------------------
stopifnot("Matrix is not square"           = nrow(fst_mat) == ncol(fst_mat))
stopifnot("Dimension does not equal n"     = nrow(fst_mat) == n)
stopifnot("Matrix is not symmetric"        = isSymmetric(fst_mat))
stopifnot("Diagonal contains non-zero"     = all(diag(fst_mat) == 0))

neg_cells <- which(fst_mat < 0 & row(fst_mat) != col(fst_mat), arr.ind = TRUE)
if (nrow(neg_cells) > 0) {
  cat("WARNING: Negative off-diagonal", stat_label, "values detected:\n")
  for (k in seq_len(nrow(neg_cells))) {
    r  <- neg_cells[k, 1]
    cl <- neg_cells[k, 2]
    cat(sprintf("  [%s, %s] = %.6f\n",
                pop_names[r], pop_names[cl], fst_mat[r, cl]))
  }
  cat("  These will produce NaN in the IBD linearisation step.",
      "Review your data.\n\n")
}

cat("Full symmetric", stat_label, "matrix:\n")
print(round(fst_mat, 6))
cat("\n")

# =============================================================================
# 5. Convert to dist object
#    This is the input currency for mantel.rtest() (ade4) and mantel() (vegan)
# =============================================================================

fst_dist <- as.dist(fst_mat)

cat(stat_label, "dist object summary:\n")
print(summary(fst_dist))
cat("\n")

# =============================================================================
# 6. Load YOUR geographic distance matrix
#
#    Uncomment ONE option. Population order is enforced against pop_names
#    automatically below regardless of which option you use.
#
#    Option A — already an R dist object saved as RDS:
#      geo_dist <- readRDS("geo_dist.rds")
#
# Option B — full symmetric matrix in a CSV, ORDER ASSUMED TO MATCH Fst' matrix
geo_mat_full <- as.matrix(read.csv("../L_virgata_pairwise-distances.csv",
                                   row.names   = 1,
                                   check.names = FALSE))

# Discard whatever names are in the CSV and impose the Fst' population names.
# This assumes row/column order in the CSV exactly matches the Fst' matrix.
dimnames(geo_mat_full) <- list(pop_names, pop_names)

geo_dist <- as.dist(geo_mat_full)
#
#    Option C — lon/lat table → Haversine distances (km):
#      library(geosphere)
#      coords        <- read.csv("coordinates.csv", row.names = 1)
#      coords        <- coords[pop_names, ]   # pre-order for safety
#      geo_mat       <- distm(as.matrix(coords[, c("lon", "lat")]),
#                             fun = distHaversine) / 1000
#      dimnames(geo_mat) <- list(pop_names, pop_names)
#      geo_dist      <- as.dist(geo_mat)
# =============================================================================


# ---- Enforce population order and validate ---------------------------------
geo_mat_full  <- as.matrix(geo_dist)
geo_labels    <- rownames(geo_mat_full)

missing_from_geo <- setdiff(pop_names, geo_labels)
if (length(missing_from_geo) > 0)
  stop(paste0(
    "The following populations are in the Fst' matrix but missing ",
    "from your geographic distance matrix:\n  ",
    paste(missing_from_geo, collapse = ", ")
  ))

extra_in_geo <- setdiff(geo_labels, pop_names)
if (length(extra_in_geo) > 0) {
  cat("NOTE: Geographic matrix contains populations not in Fst' matrix",
      "(will be dropped):\n  ", paste(extra_in_geo, collapse = ", "), "\n\n")
}

geo_mat_full <- geo_mat_full[pop_names, pop_names]
geo_dist     <- as.dist(geo_mat_full)

stopifnot(
  "Fst' and geographic dist objects are different sizes" =
    attr(fst_dist, "Size") == attr(geo_dist, "Size"),
  "Population labels do not match between Fst' and geographic matrices" =
    identical(attr(fst_dist, "Labels"), attr(geo_dist, "Labels"))
)

cat("Population order verified —", stat_label,
    "and geographic matrices match.\n\n")

# =============================================================================
# 7. Mantel test — ade4
# =============================================================================

set.seed(42)
mantel_ade4 <- mantel.rtest(fst_dist, geo_dist, nrepet = 9999)

cat("=== ade4::mantel.rtest()  | ", stat_label, "vs Geographic Distance ===\n")
print(mantel_ade4)
plot(mantel_ade4,
     main = paste0("Mantel Test (ade4)\n", stat_label,
                   " vs Geographic Distance"))

# =============================================================================
# 8. Mantel test — vegan  (also the engine used internally by adegenet)
# =============================================================================

set.seed(42)
mantel_pearson  <- vegan::mantel(fst_dist, geo_dist, 
                                 method = "pearson",  permutations = 9999)
set.seed(42)
mantel_spearman <- vegan::mantel(fst_dist, geo_dist,
                                 method = "spearman", permutations = 9999)

cat("\n=== vegan::mantel()  Pearson  | ", stat_label, "vs Geographic Distance ===\n")
print(mantel_pearson)

cat("\n=== vegan::mantel()  Spearman | ", stat_label, "vs Geographic Distance ===\n")
print(mantel_spearman)

# =============================================================================
# 9. Isolation by Distance — linearised Fst'  (Rousset 1997)
#    Recommended transformation for 2-D continuous populations:
#      Fst' / (1 - Fst')  vs  ln(geographic distance)
#
#    Diagonal: 0 / (1 - 0) = 0/1 = 0  — no special-casing needed
# =============================================================================

fst_lin_mat       <- fst_mat / (1 - fst_mat)
diag(fst_lin_mat) <- 0
fst_lin_dist      <- as.dist(fst_lin_mat)

non_finite <- !is.finite(as.vector(fst_lin_dist))
if (any(non_finite))
  warning(paste0(sum(non_finite), " non-finite value(s) in linearised Fst' dist. ",
                 "Check for Fst' values <= 0 or >= 1."))

# 1D model — use raw distance, not log distance
mantel_ibd <-vegan::mantel(fst_lin_dist, geo_dist, method = "pearson", permutations = 9999)
# (same matrices, but interpret the IBD plot x-axis as raw distance, not ln)

fst_lin_vec <- as.vector(fst_lin_dist)

# IBD plot for 1D comparison
geo_raw_vec <- as.vector(geo_dist)   # raw distance instead of log
plot(geo_raw_vec, fst_lin_vec,
     xlab = "Geographic Distance (km)",
     ylab = stat_ratio_expr,
     main = "IBD — 1D Model (Rousset 1997)\nraw distance",
     pch = 16, col = adjustcolor("steelblue", alpha.f = 0.7))
abline(lm(fst_lin_vec ~ geo_raw_vec), col = "firebrick", lwd = 2)
legend("topleft", bty = "n", cex = 0.9,
       legend = paste0("Mantel r = ", round(mantel_ibd$statistic, 4),
                       "\np = ",      mantel_ibd$signif,
                       "\nn perm = 9,999"))


# =============================================================================
# 10. Export
# =============================================================================

csv_out  <- paste0(out_prefix, "_symmetric.csv")
rds_out  <- paste0(out_prefix, "_dist.rds")

write.csv(as.data.frame(fst_mat),
          file      = csv_out,
          row.names = TRUE)

saveRDS(fst_dist,
        file = rds_out)

cat("\nExported:\n",
    " ", csv_out, " — full symmetric matrix\n",
    " ", rds_out, " — R dist object (reload with readRDS)\n")












###Run without Pigeon River
# =============================================================================
# Supplementary: remove one or more populations from both dist objects
# Run this AFTER Section 5 (fst_dist is built) and AFTER Section 6
# (geo_dist is built), but BEFORE Sections 7-9 (Mantel tests).
# =============================================================================

# Define populations to drop — add more names here if needed later
drop_pops <- "Pidgeon_River"

# Confirm the name exists before trying to drop it
missing_drop <- setdiff(drop_pops, pop_names)
if (length(missing_drop) > 0)
  stop(paste0("Population(s) not found in pop_names, check spelling:\n  ",
              paste(missing_drop, collapse = ", ")))

# ---- Drop from pop_names and fst_mat ---------------------------------------
keep_pops  <- setdiff(pop_names, drop_pops)
pop_names  <- keep_pops                          # update global pop_names
n          <- length(pop_names)                  # update global n

fst_mat    <- fst_mat[pop_names, pop_names]
fst_dist   <- as.dist(fst_mat)

# ---- Drop from geo matrix --------------------------------------------------
geo_mat_full <- as.matrix(geo_dist)
geo_mat_full <- geo_mat_full[pop_names, pop_names]
geo_dist     <- as.dist(geo_mat_full)

# ---- Drop from linearised Fst' if already built (Section 9) ---------------
if (exists("fst_lin_mat")) {
  fst_lin_mat  <- fst_lin_mat[pop_names, pop_names]
  fst_lin_dist <- as.dist(fst_lin_mat)
}

# ---- Verify ----------------------------------------------------------------
stopifnot(isSymmetric(fst_mat))
stopifnot(all(diag(fst_mat) == 0))
stopifnot(identical(attr(fst_dist,   "Labels"),
                    attr(geo_dist,   "Labels")))

cat("Dropped:", paste(drop_pops, collapse = ", "), "\n")
cat("Remaining populations (n =", n, "):\n")
cat(paste0("  ", formatC(seq_len(n), width = 2), ".  ", pop_names), sep = "\n")
cat("\n")

set.seed(42)
mantel_ade4 <- mantel.rtest(fst_dist, geo_dist, nrepet = 9999)

cat("=== ade4::mantel.rtest()  | ", stat_label, "vs Geographic Distance ===\n")
print(mantel_ade4)
plot(mantel_ade4,
     main = paste0("Mantel Test (ade4)\n", stat_label,
                   " vs Geographic Distance"))

# =============================================================================
# 8. Mantel test — vegan  (also the engine used internally by adegenet)
# =============================================================================

set.seed(42)
mantel_pearson  <- vegan::mantel(fst_dist, geo_dist, 
                                 method = "pearson",  permutations = 9999)
set.seed(42)
mantel_spearman <- vegan::mantel(fst_dist, geo_dist,
                                 method = "spearman", permutations = 9999)

cat("\n=== vegan::mantel()  Pearson  | ", stat_label, "vs Geographic Distance ===\n")
print(mantel_pearson)

cat("\n=== vegan::mantel()  Spearman | ", stat_label, "vs Geographic Distance ===\n")
print(mantel_spearman)

# =============================================================================
# 9. Isolation by Distance — linearised Fst'  (Rousset 1997)
#    Recommended transformation for 2-D continuous populations:
#      Fst' / (1 - Fst')  vs  ln(geographic distance)
#
#    Diagonal: 0 / (1 - 0) = 0/1 = 0  — no special-casing needed
# =============================================================================

fst_lin_mat       <- fst_mat / (1 - fst_mat)
diag(fst_lin_mat) <- 0
fst_lin_dist      <- as.dist(fst_lin_mat)

non_finite <- !is.finite(as.vector(fst_lin_dist))
if (any(non_finite))
  warning(paste0(sum(non_finite), " non-finite value(s) in linearised Fst' dist. ",
                 "Check for Fst' values <= 0 or >= 1."))

# 1D model — use raw distance, not log distance
mantel_ibd <-vegan::mantel(fst_lin_dist, geo_dist, method = "pearson", permutations = 9999)
# (same matrices, but interpret the IBD plot x-axis as raw distance, not ln)

fst_lin_vec <- as.vector(fst_lin_dist)

# IBD plot for 1D comparison
geo_raw_vec <- as.vector(geo_dist)   # raw distance instead of log
plot(geo_raw_vec, fst_lin_vec,
     xlab = "Geographic Distance (km)",
     ylab = stat_ratio_expr,
     main = "IBD — 1D Model (Rousset 1997)\nraw distance",
     pch = 16, col = adjustcolor("steelblue", alpha.f = 0.7))
abline(lm(fst_lin_vec ~ geo_raw_vec), col = "firebrick", lwd = 2)
legend("topleft", bty = "n", cex = 0.9,
       legend = paste0("Mantel r = ", round(mantel_ibd$statistic, 4),
                       "\np = ",      mantel_ibd$signif,
                       "\nn perm = 9,999"))

