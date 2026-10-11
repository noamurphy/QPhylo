#!/usr/bin/env Rscript
# Original visualization from saved results; does not rerun experiments.
# Copyright (c) 2026 QPhylo contributors
# SPDX-License-Identifier: MIT
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) normalizePath(args[1]) else normalizePath(".")
source <- file.path(root, "results/protein1822_comparison")
destination <- file.path(root, "qollab_upload/assets")
dir.create(destination, recursive=TRUE, showWarnings=FALSE)
library(ape)
library(phangorn)
read_topology <- function(name) {
  tree <- read.tree(file.path(source, name))
  tree$tip.label <- gsub("^'|'$", "", tree$tip.label)
  unroot(collapse.singles(tree))
}
qaoa <- read_topology("qaoa_tree.newick")
exact <- read_topology("exact_tree.newick")
nj <- read_topology("nj_tree_with_support.newick")
stopifnot(RF.dist(qaoa, exact)==0, RF.dist(qaoa, nj)==2)
splits <- read.csv(file.path(source, "splits.csv"))
stopifnot(splits$bootstrap_percent[splits$method=="Model_NJ" &
                                  splits$side1=="taxon1,taxon8"]==94)
# Highlight an unrooted split, including its complementary representation.
edge_colors <- function(tree, target, accent) {
  tips_below <- function(node) {
    if (node <= length(tree$tip.label)) return(tree$tip.label[node])
    children <- tree$edge[tree$edge[,1]==node,2]
    unlist(lapply(children, tips_below))
  }
  vapply(tree$edge[,2], function(node) {
    side <- tips_below(node)
    if (setequal(side, target) || setequal(setdiff(tree$tip.label, side), target))
      accent else "#536779"
  }, character(1))
}
png(file.path(destination, "project_header.png"), width=1920, height=1080, res=160)
par(bg="#F6F8FB", family="sans", oma=c(3.4,1.0,5.2,1.0),
    mfrow=c(1,2), mar=c(4.0,2.5,3.0,4.0), fg="#182B40", col.main="#182B40")
trees <- list(qaoa, nj)
titles <- c("QAOA + exact graph cuts", "Model-based neighbor joining")
notes <- c("Identical topology  |  RF = 0", "One split differs from graph cuts  |  RF = 2")
targets <- list(c("taxon1","taxon6","taxon7"), c("taxon1","taxon8"))
accents <- c("#2169B0", "#B85B23")
for (i in 1:2) {
  tree <- ladderize(ape::root(trees[[i]], outgroup="taxon8", resolve.root=TRUE))
  tree$node.label <- NULL
  colors <- edge_colors(tree, targets[[i]], accents[i])
  plot.phylo(tree, use.edge.length=FALSE, show.node.label=FALSE,
             edge.color=colors, edge.width=ifelse(colors==accents[i],5,2.5),
             tip.color="#182B40", cex=1.3, label.offset=0.12,
             main=titles[i])
  mtext(notes[i], side=1, line=1.0, cex=1.02, col=accents[i])
  mtext(if (i==1) "Highlighted split: taxon1 / taxon6 / taxon7" else
        "Highlighted split: taxon1 / taxon8  |  94% bootstrap", side=1, line=2.5, cex=0.88)
}
mtext("QPhylo", side=3, outer=TRUE, line=2.9, cex=2.0, font=2)
mtext("Eight proteins. Two approaches. One different tree split.",
      side=3, outer=TRUE, line=1.0, cex=1.25)
mtext("Saved pilot results  |  Topology only; branch lengths hidden  |  Display anchor is not an evolutionary root",
      side=1, outer=TRUE, line=1.3, cex=0.85, col="#536779")
dev.off()
file.copy(file.path(source, "tree_comparison.png"),
          file.path(destination, "tree_comparison.png"), overwrite=TRUE)
file.copy(file.path(root, "results/eight_taxa_comparison/optimizer_comparison.png"),
          file.path(destination, "optimizer_comparison.png"), overwrite=TRUE)
cat("Prepared original header and unchanged results figures in", destination, "\n")
