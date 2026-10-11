# Called by qphylo_nj.py. Classical model distances, NJ, support and comparison.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Expected output folder, bootstrap count and seed")
folder <- normalizePath(args[1], mustWork=TRUE)
bs_count <- as.integer(args[2]); seed <- as.integer(args[3])
if (!is.finite(bs_count) || bs_count < 1 || !is.finite(seed)) stop("Invalid bootstrap settings")
for (pkg in c("ape", "phangorn")) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop(paste("Missing R package:", pkg, "- use qphylo Rscript"))
}
labels <- readLines(file.path(folder, "labels.txt"))
data <- phangorn::read.phyDat(file.path(folder, "alignment.fasta"), format="fasta", type="AA")
if (!identical(names(data), labels)) stop("R alignment labels/order changed")
distances <- function(x) {
  dm <- phangorn::dist.ml(x, model="WAG", k=5L, shape=0.5, exclude="pairwise")
  m <- as.matrix(dm)
  if (any(!is.finite(m)) || any(m < 0) || max(abs(m-t(m))) > 1e-10 || any(diag(m) != 0))
    stop("Invalid model distance matrix; no substitution or repair")
  dm
}
check_tree <- function(tree) {
  # ape retains Newick quote delimiters; our core quotes every tip identifier.
  tree$tip.label <- vapply(tree$tip.label, function(label) {
    if (label %in% labels) return(label)
    if (startsWith(label, "'") && endsWith(label, "'"))
      return(gsub("''", "'", substr(label, 2, nchar(label)-1), fixed=TRUE))
    label
  }, character(1), USE.NAMES=FALSE)
  if (anyDuplicated(tree$tip.label) || !setequal(tree$tip.label, labels)) stop("Tree taxa mismatch")
  ape::unroot(ape::collapse.singles(tree))
}
started <- proc.time()[["elapsed"]]
dm <- distances(data)
nj <- check_tree(ape::nj(dm))
write.csv(as.matrix(dm)[labels,labels], file.path(folder, "distances.csv"), quote=TRUE)
ape::write.tree(nj, file=file.path(folder, "nj_tree.newick"))
set.seed(seed)
bs <- phangorn::bootstrap.phyDat(data, FUN=function(x) check_tree(ape::nj(distances(x))),
                                bs=bs_count, jumble=TRUE, multicore=FALSE)
if (length(bs) != bs_count) stop("Incomplete bootstrap results")
ape::write.tree(bs, file=file.path(folder, "nj_bootstrap_trees.newick"))
supported <- phangorn::plotBS(nj, bs, type="none", method="FBP")
# phangorn 2.12.1 returns FBP proportions; store/display explicit percentages.
support_proportions <- as.numeric(supported$node.label)
if (any(!is.finite(support_proportions)) || any(support_proportions < 0 | support_proportions > 1))
  stop("Unexpected FBP support scale")
supported$node.label <- as.character(100 * support_proportions)
ape::write.tree(supported, file=file.path(folder, "nj_tree_with_support.newick"))
nj_elapsed <- proc.time()[["elapsed"]] - started
trees <- list(QAOA=check_tree(ape::read.tree(file.path(folder,"qaoa_tree.newick"))),
              Exact_NMcut=check_tree(ape::read.tree(file.path(folder,"exact_tree.newick"))),
              Model_NJ=nj,
              Provided_reference=check_tree(ape::read.tree(file.path(folder,"reference_tree.newick"))))
# [[ dispatch restores shared tip labels in compressed multiPhylo objects.
for (i in seq_along(bs)) check_tree(bs[[i]])
rows <- list()
rf_matrix <- matrix(0,4,4, dimnames=list(names(trees),names(trees)))
for (i in 1:3) for (j in (i+1):4) {
  rf <- as.numeric(phangorn::RF.dist(trees[[i]], trees[[j]], normalize=FALSE,
                                     check.labels=TRUE, rooted=FALSE))
  norm <- if (length(labels) <= 3) 0 else as.numeric(phangorn::RF.dist(
    trees[[i]], trees[[j]], normalize=TRUE, check.labels=TRUE, rooted=FALSE))
  rf_matrix[i,j] <- rf_matrix[j,i] <- rf
  rows[[length(rows)+1]] <- data.frame(method1=names(trees)[i], method2=names(trees)[j],
                                      rf=rf, normalized_rf=norm)
}
comparisons <- do.call(rbind, rows)
write.csv(comparisons, file.path(folder,"rf_comparison.csv"), row.names=FALSE)
write.csv(rf_matrix, file.path(folder,"rf_matrix.csv"))

# Export actual canonical nontrivial splits for independent inspection.
split_rows <- list()
for (method in names(trees)) {
  tree <- if (method == "Model_NJ") supported else trees[[method]]
  ntip <- length(tree$tip.label)
  descendants <- function(node) {
    if (node <= ntip) return(tree$tip.label[node])
    unlist(lapply(tree$edge[tree$edge[,1] == node,2], descendants), use.names=FALSE)
  }
  for (node in tree$edge[,2]) {
    if (node <= ntip) next
    a <- sort(descendants(node)); b <- sort(setdiff(labels,a))
    if (min(length(a),length(b)) < 2) next
    if (length(a)>length(b) || (length(a)==length(b) && paste(a,collapse="|")>paste(b,collapse="|"))) {
      tmp <- a; a <- b; b <- tmp
    }
    support <- if (method == "Model_NJ") as.character(tree$node.label[node-ntip]) else NA_character_
    split_rows[[length(split_rows)+1]] <- data.frame(method=method, side1=paste(a,collapse=","),
                                                     side2=paste(b,collapse=","), bootstrap_percent=support)
  }
}
write.csv(do.call(rbind,split_rows), file.path(folder,"splits.csv"), row.names=FALSE)
write.csv(data.frame(name=c("R","ape","phangorn"),
                     value=c(R.version.string,as.character(utils::packageVersion("ape")),
                             as.character(utils::packageVersion("phangorn")))),
          file.path(folder,"versions.csv"), row.names=FALSE)
write.csv(data.frame(n_taxa=length(labels), columns=sum(attr(data,"weight")),
                     bootstraps=length(bs), nj_and_bootstrap_seconds=nj_elapsed,
                     negative_nj_edges=sum(nj$edge.length < 0)),
          file.path(folder,"nj_diagnostics.csv"), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(folder,"sessionInfo.txt"))

png(file.path(folder,"tree_comparison.png"), width=1900, height=1450, res=180)
par(mfrow=c(2,2), mar=c(1.4,1.2,3.4,1.2), oma=c(3.6,0,3.2,0))
qaoa_settings <- read.csv(file.path(folder,"qaoa_settings.csv"))
titles <- c(paste0("QAOA NMcut (p=", qaoa_settings$p[1], ", seed ", qaoa_settings$seed[1], ")"), "Exact NMcut (same graph)",
            paste0("Model NJ: WAG + Gamma (", bs_count, " bootstraps)"), "Provided reference tree")
for (i in seq_along(trees)) {
  tree <- if (i == 3) supported else trees[[i]]
  # Common display anchor only; no evolutionary root or branch-length comparison.
  tree <- ape::ladderize(ape::root(tree, outgroup=tail(labels,1), resolve.root=TRUE,
                                  edgelabel=(i==3)))
  if (i == 3) tree$node.label[1] <- ""  # The display root has no supported incoming edge.
  ape::plot.phylo(tree, use.edge.length=FALSE, show.node.label=(i==3),
                  cex=0.9, edge.width=1.8, main=titles[i], label.offset=0.08)
}
mtext(paste0("Protein1822: ", length(labels), " taxa - three reconstruction methods and supplied reference"),
      outer=TRUE, side=3, line=1, font=2, cex=1.05)
mtext("Topology views: branch lengths hidden; last taxon anchors display only.", outer=TRUE, side=1, line=1, cex=0.82)
mtext("NJ numbers are site-bootstrap percentages. Reference provenance is not established here.", outer=TRUE, side=1, line=2.1, cex=0.78)
invisible(dev.off())
print(comparisons, row.names=FALSE)
cat("NJ baseline and topology comparison completed.\n")
