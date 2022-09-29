---
title: "towcab: Topologically Overlapped Within Cluster Across Batch DEG analysis"
author: "Scott Tyler"
date: "07/12/2022"
---



# Introduction

TOWCAB stands for Topologically Overlapped Within Cluster Across Batch. And that's excatly what it does!
It goes within each cluster, and finds regions of the within-cluster topology that look "entropically"
mixed together across pairwise batches. Then it performs differential expression on a downsampled version, using the LIBRA package's approach of pseudo-bulking.


## When towcab should be used: ##
* If you've done some form of topological alignment/batch-correction and want to perform an analysis that has the underlying assumption that topologically intermingled cells are biologically similar.
* Examples include pseudotime analysis on aligned topologies
* Trying to find biological heterogeneity within a topology to benchmark alignment/batch-correction algorithms

## When wcab should be used: ##
* If you've integrated batches & are interpreting clusters as homogenous populations. 
* Example: scRNAseq on a KO vs WT & you want to say that the cluster composition is very similar, so there is not much of a difference. This might not actually be the case, because you could have within cluster DEGs across your batches. WCAB will tell you that!
* If you've done an integration on samples that you *know* are biologically different, and expect that some aspect of cell-state was erased across batch, but you want to directly measure and quantify which aspects of cell identity are differentially expressed. This is what we show in the example below; towcab is also a useful technique here, just getting you more fine grained results (although less powered due to the toplogical filtering).

## When towcab and wcab don't make sense to use ##
* If you're looking at datasets that were generated across species or technologies. In these situations, doing differential expression to compare across these doesn't actually make sense, so you won't get much out of it.
* For example, because of the wildly different distributions that you get between full length sequencing and UMI counts, you would get uninteligable results trying to do DEG analysis across batches.


# Installation

`devtools::install_github('scottyler89/towcab')`


# Examples
## Running towcab analysis

Just as an example, we'll follow along Seurat's tutorial!
https://satijalab.org/seurat/articles/integration_introduction.html

Let's follow their suggestions to integrate a stimulated & unstimulated set of human PBMCs:

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

```

## Performing (to)wcab analyses

First - we need to have the original integer count matrix (which actually wasn't distributed in their data). But we can reconstruct it with a little bit of digging into the data!

```r
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

```

Phew - now that we actually have the original count data, we can dig into the (to)wcab analyses! First we'll do just a within cluster across batch (wcab) analysis.


```r
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

## A side-note that I've noticed is that R doesn't actually correctly handly float to integer conversion correctly! Sometimes adding 1 or subtracting 1 seemingly arbitrarily. This leads to errors since towcab assumes it's reading the integer count matrix. So right after running, you see a warning saying that it's doing the coersion to integers safely using a custom little loop I wrote. I honestly have no idea why R is so finiky with this. It seems like the sort of thing that might actually be injecting lots of weird noise that people might not realize, even with something as simple as calculating the colSums - you get the wrong answer... R is strange...
```

We can also do the same analysis, but specifically looking at the toplologically overlapping areas of the dataset! Note that since these two datasets get extremely well mixed using the above integration methods, the results should come out pretty similar.


```r
## now if we want to do the full TOWCAB analysis, we need the kNN as well!
# double check the name of the graphs
print(names(immune.combined@graphs))

## and make the adjacency matrix from the network
network<-igraph::graph_from_adjacency_matrix(immune.combined@graphs$integrated_nn)

towcab_res<-run_towcab_analysis(final_exprs,
                              batch_vect,
                              clust_vect,
                              topology_filter=TRUE,
                              species="hsapiens",
                              in_topology_graph=network)

```

## Looking at the (to)wcab results

Now that we've got those results - what's in these lists?

* _*sig_table*_ : A boolean table indicating whether a was differentially expressed across batches in the noted cluster.
* _*p_table*_ : The p-values for those comparisons
* _*avg_logFC*_ : The average log fold change for those comparisons
* _*percent_sig*_ : Percent of clusters in which the given gene was differentiallye expressed
* _*num_degs*_ : For each cluster, how many genes were differentially expressed?
* _*pathway_results*_ : This is a list containing 4 objects below:
** _pathway_results$clusters_ : A list of dataframes for each cluster containing a custom annotated gProfiler results table that also annotates if it was up or down in each batch.
** _pathway_results$num_paths_sig_ : For each cluster, how many pathways were significant
** _pathway_results$geom_mean_sig_ : Of the significant pathways, what was the geometric mean of significance?
** _pathway_results$global_degs_ : Gprofiler pathway analysis on the genes that were significantly different in >50% of clusters.
* _*full_stats_list*_ : A list of objects containing all of the raw pathway analyses and inptus including custom backgrounds filtered for only the expressed genes for each cluster.

## towcab vs wcab resutls

As I mentioned before, since these were integrated via Seurat, the topologies get aligned in a way that yields a very intermingled topology. This means that the results between the wcab and towcab analyses should be fairly similar. Let's get a high level overview of what pathways were differentially expressed via the wcab and towcab approaches.

_With the code below, you can compare them, but if you want to compare different methods (liger, harmony, downsampling, etc), you can make a list with the method name corresponding to that methods (to)cab results object._

```r
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

The first thing we can look at is the collated results of pathway significance:
```r
head(collated_pathway_analyses$deg_pathway_table[,c("correction_method","cluster_id","direction","p_value","term_name")])
#       correction_method cluster_id            direction      p_value
# 1100       seurat_wcab  cluster_1 lower_in_IMMUNE_CTRL 4.662104e-46
# 2100       seurat_wcab  cluster_1 lower_in_IMMUNE_CTRL 4.662104e-46
# 3100       seurat_wcab  cluster_1 lower_in_IMMUNE_CTRL 1.592834e-45
# 317        seurat_wcab  cluster_1 lower_in_IMMUNE_CTRL 1.821440e-42
# 425        seurat_wcab  cluster_1 lower_in_IMMUNE_CTRL 2.336789e-42
# 912        seurat_wcab  cluster_1        all_sig_diffs 1.283823e-40
#                                                                      term_name
# 1100                                      response to external biotic stimulus
# 2100                                                response to other organism
# 3100                                               response to biotic stimulus
# 317                                                              Immune System
# 425  biological process involved in interspecies interaction between organisms
# 912                                                 response to other organism
```

These results are great! It shows us that in cluster_1, genes related to the immune response to pathogens is lower in the control compared to the stimulated condition. This is what you would expect when using a class-2 algorithm to integrate datasets. There's meaninful biological signal that is heterogeneous in expression along the topology. This means that two adjacent cells within the toplogy can actually be in quite different biological states making interpretation of the data something that must be done with lots of care and caution.

## A rapid interpretation for comparing two integration/batch-correction methods
There is also a meta-analysis of pathway results is performed by comparing how similar the intersection genes are between all pathways in your input that were significantly different. This gives a graph network for which nodes are the pathways & edges connect differentially expressed pathways whose gene-lists were quite similar. This graph is then clustered based on Louvain modularity. The file titled "<out_dir>/results/pathway_degs/WInClust_AcrossBatch_DEG_pathway_meta_analysis_02.png" and should look something like this:

![pathway networks](https://github.com/scottyler89/towcab/blob/main/towcab_results/results/pathway_DEGs/WInClust_AcrossBatch_DEG_pathway_meta_analysis_02.png?raw=true)

The colors correspond to the pathway clusters (not the clustering of the cells!). This is just to quickly get a high-level view of which pathways are significantly different in each of the methods (in this case wcab and towcab both from the same Seurat run).

Now how do we interpret them? Well we can look back at that returned object & see the number of pathways and percentage of pathways that within the pathway clusters that are significantly different:

```r
print(head(collated_pathway_analyses$path_sig_res))
#  correction_method pathway_cluster num_sig_pathways percent_sig_pathways
#1       seurat_wcab               1              300            0.9316770
#2     seurat_towcab               1              215            0.6677019
#3       seurat_wcab               2               97            0.8738739
#4     seurat_towcab               2               93            0.8378378
#5       seurat_wcab               3              103            0.7803030
#6     seurat_towcab               3               94            0.7121212
```

If we want to check out these results visually, there are a few plots in the output directory as well:

![quantified deg num and degree](https://github.com/scottyler89/towcab/blob/main/towcab_results/results/combined_DEG_analysis_1.png?raw=true)
![path cluster sig](https://github.com/scottyler89/towcab/blob/main/towcab_results/results/pathway_DEGs/NonTechnical_WInClust_AcrossBatch_DEG_analysis_01.png?raw=true)
![stacked percent of clusters sig](https://github.com/scottyler89/towcab/blob/main/towcab_results/results/pathway_DEGs/NonTechnical_WInClust_AcrossBatch_DEG_analysis_02.png?raw=true)
![stacked number of paths](https://github.com/scottyler89/towcab/blob/main/towcab_results/results/pathway_DEGs/NonTechnical_WInClust_AcrossBatch_DEG_analysis_03.png?raw=true)


Now what about the biological interpretation? Check out the "towcab_results/results/pathway_DEGs/DEG_pathway_descriptions.txt" file. This is a print out of each pathway_cluster & which pathways are contained in it sorted alphabetically. Personally, I've found that to be the easiest way to rapidly get to interpration.

In this case we see:
* *1*: Immune response, interferons, antigen processing/presentation, etc
* *2*: Proteasome, cell stress, translation (perhaps a bit of replication going on here)
* *3*: Things related to glandular cells, lymphocytes, and germinal centers
* *4*: MHC-II antigen processing, lysosomes, & granules 
* *5*: IRF binding sites in the promoter, and some immune pathways & response to stimuli

## Identification of locations with errant integrations

We can also look at specific genes of interest that should not be topologically overlapped using the identify_any_mixed_clusters function. This particular run below is more for illustrative purposes. See our manuscript to see how we used this function to identify the errant integration of urothelial- and B-cells.

```r
## first make a list of cell-type or -state markers that should not be mixed together

## since this is a silly example with known biological differences, we'll look for the mixing of inflammatory cytokines in T and B cells, recognizing that the cytokines very well may belong in these populations. It's only for illustrative purposes. But this is the function that led us to identify the errant mixture of urothelial cells and B-cells in our manuscript.

marker_list<-list()
marker_list[["T-cell"]]<-c("CD4","CD3E","CD8A")
marker_list[["B-cell"]]<-c("CD79A","CD79B","CD19","VPREB3")
marker_list[["Macrophage"]]<-c("CD14")
marker_list[["inflam"]]<-c("IFNG", "IL1B", "CXCL10", "IL8")


mixed_results<-identify_any_mixed_clusters(towcab_res$p_table, marker_list, alpha_cutoff = 1e-3)## default is super-conservative at 1e-6

## Here we se that the there are 3 clusters that show heterogeneous expression across batches for expression of CD14, IL1B, and CXCL10
# 
# Macrophage inflam
# [1] "cluster_0" "cluster_8" "cluster_3"
#   cluster_0
#       CD14 IL1B CXCL10 IL8 
#   cluster_8
#       CD14 IL1B CXCL10 IL8 
#   cluster_3
#       CD14 IL1B CXCL10 IL8 
# inflam Macrophage
# [1] "cluster_0" "cluster_8" "cluster_3"
#   cluster_0
#       IL1B CXCL10 IL8 CD14 
#   cluster_8
#       IL1B CXCL10 IL8 CD14 
#   cluster_3
#       IL1B CXCL10 IL8 CD14 



```



