# Standalone scientific figures from validated CSV exports; base R only.
args <- commandArgs(trailingOnly=TRUE)
folder <- args[1]
runs <- read.csv(file.path(folder, "runs.csv"), stringsAsFactors=FALSE)
curves <- read.csv(file.path(folder, "evaluations.csv"), stringsAsFactors=FALSE)
finals <- read.csv(file.path(folder, "final_cuts.csv"), stringsAsFactors=FALSE)
methods <- sort(unique(runs$optimizer))
colors <- c("seeded_coordinate_search"="#0072B2", "COBYLA"="#D55E00")[methods]
short <- function(x) ifelse(x == "seeded_coordinate_search", "Coordinate search", x)
root <- subset(curves, n == root_n & target_size == floor(root_n/2) & attempt == 0 & phase == "training")
root_final <- subset(finals, n == root_n & target_size == floor(root_n/2) & attempt == 0)
paired_panel <- function(data, field, title, ylabel, baseline=NULL) {
  values <- data[[field]]
  good <- is.finite(values)
  yrange <- range(c(values[good], baseline, 0), na.rm=TRUE)
  if (diff(yrange) == 0) yrange <- yrange + c(0, 1)
  plot(NA, xlim=c(0.6, length(methods)+0.4), ylim=yrange * c(1, 1.12),
       xaxt="n", xlab="", ylab=ylabel, main=title)
  axis(1, at=seq_along(methods), labels=short(methods), cex.axis=0.85)
  for (seed in unique(data$seed)) {
    sub <- data[data$seed == seed & is.finite(data[[field]]), ]
    sub <- sub[order(match(sub$optimizer, methods)), ]
    if (nrow(sub) == 2) lines(match(sub$optimizer, methods), sub[[field]], col="#BBBBBB")
  }
  for (method in methods) {
    sub <- data[data$optimizer == method & is.finite(data[[field]]), ]
    points(rep(match(method, methods), nrow(sub)), sub[[field]],
           col=colors[method], pch=19, cex=1.1)
  }
  if (!is.null(baseline)) {
    abline(h=baseline, lty=3, col="#666666")
    legend("topright", "Uniform sampling (analytic)", lty=3, col="#666666", bty="n", cex=0.7)
  }
}
draw <- function() {
  par(mfrow=c(2,2), mar=c(4.2,4.2,3,1.2), oma=c(2.5,0,3,0), las=1)
  if (nrow(root)) {
    plot(NA, xlim=range(root$evaluation+1), ylim=range(root$best_training_so_far),
         xlab="Training circuit evaluations", ylab="Best scaled penalized mean energy",
         main="Root balanced-cut optimization")
    for (method in methods) {
      for (seed in unique(root$seed[root$optimizer == method])) {
        sub <- root[root$optimizer == method & root$seed == seed, ]
        lines(sub$evaluation+1, sub$best_training_so_far, col=adjustcolor(colors[method], 0.6), lwd=1.6)
      }
    }
    legend("topright", short(methods), col=colors, lwd=2, bty="n", cex=0.8)
  } else { plot.new(); title("No root training measurements") }
  root_final$optimal_percent <- 100 * root_final$optimal_shot_rate
  paired_panel(root_final, "optimal_percent", "Fresh final optimal-cut frequency", "Optimal shots (%)",
               if(nrow(root_final)) 100*root_final$uniform_optimal_rate[1] else NULL)
  paired_panel(runs, "rf_to_exact", "Tree agreement with exact NMcut", "Unrooted RF distance")
  paired_panel(runs, "submitted_jobs", "Total tree measurement cost", "Submitted circuit jobs")
  kind <- if ("dataset_kind" %in% names(runs) && runs$dataset_kind[1] == "protein_derived") "protein-derived" else "synthetic"
  mtext(paste0("QAOA optimizer pilot: ", runs$n[1], " ", kind, " taxa, ideal simulator"), outer=TRUE, side=3, line=1, font=2, cex=1.15)
  mtext(paste0("Recorded runs: ", nrow(runs), "; failures: ", sum(runs$status == "failed"),
               ". Paired initialization seeds joined where available; no significance claim."),
        outer=TRUE, side=1, line=1, cex=0.8)
}
png(file.path(folder,"optimizer_comparison.png"), width=1800, height=1300, res=180); draw(); dev.off()
pdf(file.path(folder,"optimizer_comparison.pdf"), width=10, height=7.2); draw(); dev.off()
cat("Saved optimizer_comparison.png and .pdf\n")
