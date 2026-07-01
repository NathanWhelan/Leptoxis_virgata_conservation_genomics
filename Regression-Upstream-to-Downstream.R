# ============================================================
#  CONFIGURATION — edit these for your dataset
# ============================================================
rm(list=ls())
csv_file      <- "distanceFromDownstream_geneticDiversity.csv"
predictor     <- "Distance"
site_col      <- "Site"
exclude_sites <- c("Pigeon_River")  # set to c() for none
# ============================================================

setwd("~/snails/L_virgata/population-genomics/virgata_m3M3n3_R80_maf025/")

data_full <- read.csv(csv_file)
data_full$Distance
metrics   <- setdiff(names(data_full)[sapply(data_full, is.numeric)], predictor)

datasets <- list(All_Sites = data_full)
if (length(exclude_sites) > 0)
  datasets[[paste0("Excl_", paste(exclude_sites, collapse = "_"))]] <-
  data_full[!data_full[[site_col]] %in% exclude_sites, ]

run_regressions <- function(data, label) {
  
  models <- setNames(lapply(metrics, function(y)
    lm(as.formula(paste(y, "~", predictor)), data = data)), metrics)
  
  # Print summaries and build stats table in one pass
  stats <- lapply(metrics, function(y) {
    m <- models[[y]]; s <- summary(m)
    p  <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    sp <- cor.test(data[[predictor]], data[[y]], method = "spearman", exact = FALSE)
    cat("\n---", y, "~", predictor, "---\n"); print(s)
    round(c(Intercept = unname(coef(m)[1]), Slope = unname(coef(m)[2]),
            R2 = s$r.squared, Adj_R2 = s$adj.r.squared,
            F_stat = unname(s$fstatistic[1]), P_value = p,
            Spearman_rho = unname(sp$estimate), Spearman_p = sp$p.value), 4)
  })
  write.table(
    data.frame(Metric = metrics, as.data.frame(do.call(rbind, stats))),
    file = paste0("summary_table_", label, ".txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  # Scatter plots
  png(paste0("regression_plots", label, ".png"), width = 1200, height = 1000)
  par(mfrow = c(ceiling(length(metrics) / 2), 2), mar = c(5, 5, 4, 2))
  for (y in metrics) {
    m <- models[[y]]; s <- summary(m)
    p <- pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    plot(data[[predictor]], data[[y]], pch = 19, col = "steelblue", cex = 1.3,
         xlab = predictor, ylab = y, main = paste0(y, " ~ ", predictor, "\n(", label, ")"))
    abline(m, col = "firebrick", lwd = 2)
    text(data[[predictor]], data[[y]], labels = data[[site_col]], cex = 0.55, pos = 3)
    legend("topright", bty = "n", cex = 0.9,
           legend = c(paste0("R² = ", round(s$r.squared, 3)),
                      ifelse(p < 0.001, "p < 0.001", paste0("p = ", round(p, 3)))))
  }
  dev.off()
  
  # Diagnostic plots
  png(paste0("diagnostics_", label, ".png"), width = 1600, height = 400 * length(metrics))
  par(mfrow = c(length(metrics), 4), mar = c(4, 4, 3, 2))
  for (y in metrics) plot(models[[y]], main = paste(y, "~", predictor, "|", label))
  dev.off()
  
  cat("Outputs saved with tag:", label, "\n")
}

for (label in names(datasets)) run_regressions(datasets[[label]], label)

