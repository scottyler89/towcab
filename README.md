---
title: "The PercentMaxDiff user's guide"
author: "Scott Tyler"
date: "11/05/2020"
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: vignette
vignette: >
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{"The towcab user's guide"}
---



# Introduction




# Installation

`devtools::install_github('scottyler89/towcab')`


# Examples
## Running towcab analysis


```r
library("towcab")
library("Seurat")
library("SeuratData")
library("patchwork")

# install dataset
InstallData("ifnb")

# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

################
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

################
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)


###############
## import some stuffs for the TOWCAB analysis


## TOWCAB requires the original counts, so we'll need to regenerate those
# to illustrate how, here's a demo that this is actually the natural log of the
# counts per 10k normalized version of the data
colSums(as.matrix(exp(1)**immune.combined@assays$RNA[,1:10]-1))
# so let's undo that
exprs<-as.matrix(exp(1)**immune.combined@assays$RNA[,]-1)
## now we have to figure out what the "loading factors" were for each column.
## to do this we assume that the minimum non-zero count was a count of 1.
## Pretty reasonable I think. Never seen a UMI based dataset this wasn't true for.
nz_col_min<-function(in_vect){
  in_vect<-in_vect[which(in_vect>0)]
  return(min(in_vect))
}
nz_mins<-apply(exprs,2,nz_col_min)
for (i in seq(1,dim(exprs)[2])){
  exprs[,i]<-exprs[,i]/nz_mins[i]
}
## now we have nice & pretty integers again =)
print(exprs[1:10,1:10])
post_unnorm_nz_mins<-apply(exprs,2,nz_col_min)
post_unnorm_colsums<-colSums(exprs)
## except (on my machine at least...) changing it to integer mode 
## changes the damn numbers! Why oh why R... Shall I count the ways in which I hate thee
final_exprs<-matrix(data=exprs,nrow=dim(exprs)[1],ncol=dim(exprs)[2])
colnames(final_exprs)<-colnames(exprs)
rownames(final_exprs)<-rownames(exprs)
mode(final_exprs)<-"integer"
## ugh.. go and fix it
for (i in seq(1,dim(final_exprs)[2])){
  ## I checked, there aren't any non-integers, but in the coersion from float to int
  ## R still screws it up... my hair is falling out...
  final_exprs[,i]<-as.integer(round(exprs[,i],0))
}

## now we have all the info we need to run (TO)WCAB analyses!
## That's basically the raw expression matrix, the batch info, and the cluster info
batch_vect <- unlist(immune.combined@meta.data["orig.ident"])
clust_vect <- unlist(immune.combined@active.ident)



## if you want to do just the within cluster DEG analysis, you don't need the graph though:
wcab_res<-run_towcab_analysis(final_exprs,
                              batch_vect,
                              clust_vect,
                              topology_filter=FALSE,
                              species="hsapiens")## these are gprofiler codes

```

## Looking at the (to)wcab results

The resultant object is a list with several useful metrics that describe how similar batches are to each other based on their clustering results.

* _*example thing*_ : words and things



```

## now if we want to do the full TOWCAB analysis, we need the kNN as well!
# double check the name of the graphs
print(names(immune.combined@graphs))
## and make the adjacency matrix from the network of your choice
network<-igraph::graph_from_adjacency_matrix(immune.combined@graphs$integrated_nn)


towcab_res<-run_towcab_analysis(final_exprs,
                              batch_vect,
                              clust_vect,
                              topology_filter=TRUE,
                              species="hsapiens",
                              in_topology_graph=network)

```

## towcab vs wcab resutls

The resultant object is essentially the same - the only difference being that during the run, topological filtering was performed.


## comparing methods & getting a high-level overview

If you want to compare different methods, you can make a list with the method name corresponding to that methods (to)cab results object. 

```
## now we can get a high-level overview or compare different methods
towcab_res_lists<-list()
towcab_res_lists[["seurat_wcab"]]<-wcab_res
towcab_res_lists[["seurat_towcab"]]<-towcab_res

# Note that you could also feed in the results for different methods of batch correction!
## We also need to give each method's results a color for plotting
method_colors<-list()
method_colors[["seurat_wcab"]]="red"
method_colors[["seurat_towcab"]]="blue"

out_dir<-paste(getwd(),"towcab_results",sep="/")
collated_pathway_analyses<-analyze_all_pathway_results(towcab_res_lists,out_dir,method_colors,organism="hsapiens")

```

## interpreting the results

There could be a lot to digest here, so we'll start with the large scale take-away

```

```
