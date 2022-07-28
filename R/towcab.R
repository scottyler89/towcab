###################################################################################
###################################################################################
###################################################################################

'%!in%' <- function(x,y)!('%in%'(x,y))


get_ds_transcript_vect<-function(in_col_vect1, target_depth){
  trans_vect<-get_sim_transcripts(in_col_vect1, target_depth)
  new_vect<-rep(0,length(in_col_vect1))
  for (i in trans_vect)  {
    new_vect[i]<-new_vect[i]+1
  }
  return(new_vect)
}


get_sim_transcripts<-function(in_col_vect, target_depth){
  print(sum(in_col_vect))
  in_col_vect<-as.integer(in_col_vect)
  message("target_depth",target_depth)
  print(sum(in_col_vect))
  nz_rows<-which(in_col_vect>0)
  transcript_vect<-c()
  for (i in 1:length(nz_rows)){
    temp_id<-nz_rows[i]
    transcript_vect<-c(transcript_vect,rep(temp_id,in_col_vect[nz_rows[i]]))
  }
  #message(length(transcript_vect)," & ",sum(in_col_vect))
  #print(table(transcript_vect))
  if (length(transcript_vect)<target_depth){
    print(paste(c(length(transcript_vect),"transcripts total for depth of ",target_depth), collapse=" "))
  }
  return(sample(transcript_vect,target_depth))
}


#' Downsample columns of a matrix
#'
#' This function takes in a matrix of counts & will downsample 
#' each column such that their colSums are equivalent. Note that you can't just
#' randomly subtract 1 from cells until they're equivalent. To recapitulate the
#' actual biases from Poisson sampling, we actually create the "bag of marbles"
#' and sample each until we get the the target depth. You can either pass in 
#' the target depth, or if not, it will be automatically set to the minimum of 
#' the colSums, or if you are doing an integration & want to downsample 
#' different datasets, you should first FILTER OUT the columns whose colSums 
#' are less than the cutoff, then set this value to the cutoff.
#'
#' @param in_mat An input matrix
#' @param target_depth This can either be nothing or a specific target to 
#' @param min_count The minimum number of observed genes for a cell to be included
#' @param max_sum The maximum number of counts per cell for a cell to be included
#' @param max_count Maximum count of observed genes for cells to be included
#' @param rm_less Remove cells with less than the total count target_depth
#' @param parallel Run in parallel? Default=T
#' @param cores If running in parallel, how many cores to use, default=NULL=all available
#' @param quiet Runs a bit quieter, default=F
#' which we will normalize to. Note that this should be the minimum of colsum 
#' across all datasets.
#' @return A matrix of the normalized matrix
#' @examples
#'    in_mat<-matrix(rnbinom(1000*20000,.1,c(0.5,0.5,0.5)),ncol=1000,nrow=20000)
#'    out_mat<-downsample_mat(in_mat)
#' @importFrom parallel detectCores makeCluster parApply stopCluster clusterExport
#' @name local_downsample_mat
local_downsample_mat<-function(in_mat, 
                               target_depth=NULL, 
                               min_count=0, 
                               max_sum=Inf, 
                               max_count=Inf, 
                               rm_less=TRUE, 
                               do_parallel=TRUE, 
                               cores=NULL, 
                               quiet=FALSE){
  #message("do_parallel: ",do_parallel)
  all_cell_sums<-colSums(in_mat)
  #message("min colSums:",min(all_cell_sums))
  in_mat_copy<-in_mat
  in_mat_copy[in_mat_copy>0]<-1
  all_cell_counts<-colSums(in_mat_copy)
  rm(in_mat_copy)
  gc()
  if (is.null(target_depth)){
    target_depth<-min(colSums(in_mat))
  }
  keep_idxs <- which(all_cell_counts >= min_count & 
                     all_cell_counts <= max_count & 
                     all_cell_sums >= target_depth & 
                     all_cell_sums <= max_sum )
  in_mat <- in_mat[,keep_idxs]
  print(min(colSums(in_mat)))
  if (do_parallel){
    if (is.null(cores)){
        cores<- detectCores()
    }
    cl <- makeCluster(cores)
    clusterExport(cl, c("get_ds_transcript_vect", "get_sim_transcripts"),
              envir=environment())
    if (!quiet){    
      message("        running with ",cores," cores")
    }
    out_mat <- parApply(cl, X=in_mat, MARGIN=2, FUN=get_ds_transcript_vect, target_depth=target_depth)
    stopCluster(cl)
  } else {
    out_mat <- apply(in_mat, 2, function(x) get_ds_transcript_vect(x, target_depth))
  }
  rownames(out_mat)<-rownames(in_mat)
  colnames(out_mat)<-colnames(in_mat)
  return(out_mat)
}
#downsample_mat(matrix(rep(5,50),ncol=5,nrow=10))

##############################################
############################################################################
###################################################################################
###################################################################################

offset_reps<-function(meta, num_pseudo_replicates=3){
    # print(head(meta))
    batch_vect<-as.character(meta[,"label"])
    unique_batch_vect<-sort(unique(batch_vect))
    # message("unique_batch_vect:")
    # print(unique_batch_vect)
    num_batches<-length(unique_batch_vect)
    num_reps<-num_pseudo_replicates#length(unique(as.character(meta$replicate)))
    meta[,"replicate"]<-as.numeric(meta[,"replicate"])
    all_possible<-seq(num_pseudo_replicates)
    batch_rep_lookup_table<-list()
    if (num_batches==1){
        return(meta)
    } else {
        batch_rep_lookup_table[[unique_batch_vect[1]]]<-as.character(seq(num_pseudo_replicates))
        for (temp_idx in seq(2,num_batches)){
            temp_possible<-c(seq(num_pseudo_replicates))+((temp_idx-1)*num_reps)
            batch_rep_lookup_table[[unique_batch_vect[temp_idx]]]<-as.character(temp_possible)
            all_possible<-c(all_possible,temp_possible)
            temp_batch<-unique_batch_vect[temp_idx]
            temp_batch_idxs<-which(batch_vect==temp_batch)
            meta[temp_batch_idxs,"replicate"]<-meta[temp_batch_idxs,"replicate"]+(temp_idx-1)*num_reps
        }
    }
    ## remove the inappropriate batches from the levels because this
    ## screws up making the contrast matrix
    meta[,"label"]<-factor(meta[,"label"],
                           levels=unique_batch_vect)
    meta[,"replicate"]<-factor(meta[,"replicate"],
                               levels=all_possible)
    # print(table(meta[,"replicate"],meta[,"label"]))
    # print(names(meta))
    meta[,"cell_type"]<-as.factor(meta[,"cell_type"])
    out_list<-list()
    out_list[["meta"]]<-meta
    out_list[["batch_rep_lookup_table"]]<-batch_rep_lookup_table
    return(out_list)
}



#' @importFrom igraph neighbors as_edgelist graph_from_data_frame
get_per_cell_percentage<-function(local_toplology_graph_in, meta_in, do_shuffle=FALSE){
    ### If we're shuffling to get the background, then do that
    if (do_shuffle){
        temp_edge_df<-igraph::as_edgelist(local_toplology_graph_in)
        temp_edge_df[,1]<-sample(temp_edge_df[,1])
        temp_edge_df[,2]<-sample(temp_edge_df[,2])
        local_toplology_graph<-igraph::graph_from_data_frame(temp_edge_df, directed=FALSE)
    } else {
        local_toplology_graph<-local_toplology_graph_in
    }
    ### now calculate the percentage 
    # print(head(meta_in))
    all_local_verts<-as.character(meta_in$cell_id)#
    # print(head(all_local_verts))
    all_local_verts2<-as.character(unlist(igraph::V(local_toplology_graph)$name))
    all_local_verts<-sort(unique(c(all_local_verts, all_local_verts2)))
    percent_host_batch<-c()
    for (vert in all_local_verts){
        if (vert %!in% all_local_verts2){
            temp_pcent_host_batch<-0
        } else {
            temp_neigh<-unlist(igraph::neighbors(local_toplology_graph, vert, mode="all"))
            total_neigh<-length(temp_neigh)
            host_batch<-meta_in[vert,"label"]
            batch_table<-table(meta_in[temp_neigh,"label"])
            temp_pcent_host_batch<-batch_table[host_batch]/total_neigh
        }
        if (is.nan(temp_pcent_host_batch)){
            temp_pcent_host_batch<-0
        }
        if (is.na(temp_pcent_host_batch)){
            temp_pcent_host_batch<-0
        }
        percent_host_batch<-c(percent_host_batch, temp_pcent_host_batch)
    }
    out_df<-data.frame(node=as.character(all_local_verts),
                       percent_host_batch=percent_host_batch)
    rownames(out_df)<-as.character(all_local_verts)
    # if (do_shuffle==FALSE){
    #     message("head percent host batch df:")
    #     print(head(out_df))
    # }
    return(out_df)
}


par_get_per_cell_percentage<-function(dummy,local_toplology_graph_in, meta_in, do_shuffle){
    return(get_per_cell_percentage(local_toplology_graph_in=local_toplology_graph_in, meta_in=meta_in, do_shuffle=do_shuffle))
}


#' run_towcab_analysis
#' @description \code{run_towcab_analysis} Performs the full towcab analysis
#' @param local_temp_mat input matrix for only a single cluster
#' @param local_meta the meta data frame for only the cells contained in the pertinent cluster. A couple important notes: needs columns called "cell_id" and "label" which are consistent cell_ids across the three inputs, matching the node names in the network, colnames in the matrix, and rownames in the meta. "Label" is the batch ID.
#' @param local_topology_graph the kNN graph used for clustering (may either be subset of this cluster or the full network)
#' @param num_boots The number of bootstrap shuffled versions of the kNN used to estimate each cell's null distribution for what that cell's (default=100)
#' @param abs_z_cut The cutoff of the absolute value of the Z away from the mean for a cell to be excluded as non-overlapping. For example, a cell with a cross-batch connection fraction z-score > 1.5 or < -1.5 away from the bootstrap shuffled edges from the null background (default=1.5)
#' @param clip_val The lower clip value to include a sample; the upper clip is calculated as 1-clip_val. For example, if a sample doesn't have at least 10% connections across batch, exclude it because it could be a spurios thing that doesn't represent broad patterns. Same at 90%. (default = 0.10)
#' @return a filtered version of the input matrix and meta containing only the topologically overlapped subset
#' @importFrom igraph induced_subgraph
#' @importFrom parallel makeCluster detectCores clusterExport parLapply stopCluster
#' @name do_topology_filter
#' @export
do_topology_filter<-function(local_temp_mat, local_meta, local_topology_graph, 
                             num_boots=100, abs_z_cut=1.5, clip_val=0.1){
    local_toplology_graph<-induced_subgraph(local_topology_graph, colnames(local_temp_mat))
    # print(head(local_meta))
    # print(colnames(local_temp_mat))
    real_percent_host<-get_per_cell_percentage(local_toplology_graph_in=local_toplology_graph,
                                               meta_in=local_meta,
                                               do_shuffle=FALSE)
    #print("getting boots")
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("%!in%","get_per_cell_percentage","as_edgelist","neighbors","V","graph_from_data_frame"), envir=environment())
    boots<-parLapply(cl, 
                    X=seq(num_boots),
                    fun=par_get_per_cell_percentage, 
                    local_toplology_graph_in=local_toplology_graph,
                    meta_in=local_meta,
                    do_shuffle=TRUE)
    stopCluster(cl)
    null_mat<-get_null_mat(boots)
    include_vect<-rep(TRUE,dim(null_mat)[1])
    mean_vect<-c()
    sd_vect<-c()
    # message("null_mat dims:")
    # print(dim(null_mat))
    # print(dim(real_percent_host))
    for (i in seq(dim(null_mat)[1])){
        temp_mean<-mean(na.omit(null_mat[i,]))
        mean_vect<-c(mean_vect,temp_mean)
        temp_sd<-sd(na.omit(null_mat[i,]))
        sd_vect<-c(sd_vect, temp_sd)
        temp_pct<-real_percent_host[i,"percent_host_batch"]
        if (is.na(temp_pct) | is.nan(temp_pct)){
            include_vect[i]<-FALSE
        } else {
            if ( temp_pct < (temp_mean-abs_z_cut*temp_sd) ){
                include_vect[i]<-FALSE
            }
            if (temp_pct>(temp_mean+abs_z_cut*temp_sd)) {
                include_vect[i]<-FALSE
            }
            ## add some clipping
            if ( temp_pct < clip_val ) {
                include_vect[i]<-FALSE
            }
            if (temp_pct > (1-clip_val)) {
                include_vect[i]<-FALSE
            }
        }
    }
    # print(real_percent_host)
    # print(dim(real_percent_host))
    # print(length(mean_vect))
    # print(length(sd_vect))
    # print(length(include_vect))
    keep_cells<-as.character(real_percent_host[which(include_vect),"node"])
    real_percent_host<-as.data.frame(cbind(real_percent_host,mean_vect,sd_vect,include_vect))
    message("keeping ",length(keep_cells)," after topology filtering")
    #print(head(real_percent_host,20))
    out_list<-list()
    out_list[["temp_mat"]]<-local_temp_mat[,keep_cells]
    out_list[["meta"]]<-local_meta[keep_cells,]
    return(out_list)
}




#################################################################################################


check_cross_tabs<-function(cross_tab, num_pseudo_replicates=3, min_cells_per_pseudo_rep=3){
    all_good<-TRUE
    #print(cross_tab)
    for (batch_idx in seq(dim(cross_tab)[2])){
        cur_rep_idxs<-seq(num_pseudo_replicates)+((batch_idx-1)*num_pseudo_replicates)
        batch<-colnames(cross_tab)[batch_idx]
        # message("cross_tab[cur_rep_idxs,batch_idx]")
        # print(cross_tab[cur_rep_idxs,batch_idx])
        # message("min(cross_tab[cur_rep_idxs,batch_idx]) < min_cells_per_pseudo_rep")
        # message(min(cross_tab[cur_rep_idxs,batch_idx]) < min_cells_per_pseudo_rep)
        if (min(cross_tab[cur_rep_idxs,batch_idx]) < min_cells_per_pseudo_rep) {
            all_good<-FALSE
        }
    }
    #message("all_good:",all_good)
    return(all_good)
}


try_to_update_meta<-function(meta, cross_tab, num_pseudo_replicates=3, min_cells_per_pseudo_rep=3){
    ## wow - sometimes the cross_tab rows come out in non-sorted order...
    # print("before resort")
    # print(cross_tab)
    cross_tab<-cross_tab[order(as.integer(rownames(cross_tab))),]
    # print("after resort")
    # print(cross_tab)
    for (batch_idx in seq(dim(cross_tab)[2])){
        cur_rep_idxs<-seq(num_pseudo_replicates)+((batch_idx-1)*num_pseudo_replicates)
        batch<-colnames(cross_tab)[batch_idx]
        if (min(cross_tab[cur_rep_idxs,batch_idx]) < min_cells_per_pseudo_rep){
            ## find the rep that has the highest
            lowest_row_idx<-c(cur_rep_idxs[which(cross_tab[cur_rep_idxs,batch_idx]==min(cross_tab[cur_rep_idxs,batch_idx]))])[1]
            lowest_rep_id<-rownames(cross_tab)[lowest_row_idx]
            highest_row_idx<-c(cur_rep_idxs[which(cross_tab[cur_rep_idxs,batch_idx]==max(cross_tab[cur_rep_idxs,batch_idx]))])[1]
            highest_rep_id<-rownames(cross_tab)[highest_row_idx]
            ## pick out a cell at random from this batch & highest label category
            message("highest_rep_idx:",highest_rep_id)
            message("highest_rep_id:",highest_rep_id)
            message("batch:",batch)
            right_batch_bool<-as.integer(as.character(meta[,"label"])==as.character(batch))
            print(right_batch_bool)
            right_rep_bool<-as.integer(as.character(meta[,"replicate"])==as.character(highest_rep_id))
            print(right_rep_bool)
            all_pertinent_cell_idxs <- which((right_rep_bool+right_batch_bool)==2)
            #all_pertinent_cell_idxs <- which(as.character(meta[,"label"])==as.character(batch) & as.character(meta[,"replicate"])==as.character(highest_rep_id))
            message("all_pertinent_cell_idxs")
            print(all_pertinent_cell_idxs)
            if (length(all_pertinent_cell_idxs)>0){
                random_cell_idx<-sample(c(all_pertinent_cell_idxs),1)
                print(meta[random_cell_idx,])
                meta[random_cell_idx,"replicate"]<-lowest_rep_id
                print(meta[random_cell_idx,])
            } else {
                print(cross_tab)
            }
        }
    }
    return(meta)
}


make_sure_enough_cells_in_each_rep<-function(meta, num_pseudo_replicates=3, min_cells_per_pseudo_rep=3, offset_lookup=NULL){
    # cross_tab<-table(meta$replicate,meta$label)
    # print("CT1")
    # print(cross_tab)
    # print(summary(meta))
    all_rep_ids<-c()
    for (temp_batch_id in sort(as.character(unique(meta$label)))){
        all_rep_ids<-c(all_rep_ids,offset_lookup[[temp_batch_id]])
    }
    meta$replicate<-factor(meta$replicate,levels=all_rep_ids)#seq(length(unique(meta$label))*num_pseudo_replicates))#sort(as.integer(as.numeric(as.character(unique(meta$replicate))))))
    cross_tab<-table(meta$replicate,meta$label)
    # print("CT2")
    # print(cross_tab)
    # print(summary(meta))
    counter<-0
    #print("first check:")
    res<-check_cross_tabs(cross_tab,min_cells_per_pseudo_rep=min_cells_per_pseudo_rep)
    #message("first_check result:",res)
    while (check_cross_tabs(cross_tab,min_cells_per_pseudo_rep=min_cells_per_pseudo_rep)==FALSE){
        counter<-counter+1
        #message("updating meta...")
        meta<-try_to_update_meta(meta, cross_tab, 
                                 num_pseudo_replicates=num_pseudo_replicates,
                                 min_cells_per_pseudo_rep=min_cells_per_pseudo_rep)
        cross_tab<-table(meta$replicate,meta$label)
        if (counter>20){
            print(cross_tab)
            stop("something's wrong....")
        }
    }
    return(meta)
}


########################################################################################################
## pathway analysis


get_pcnt_deg_df<-function(deg_table){
    #print(colSums(deg_table))
    percent_of_times_significant <- rowSums(deg_table)/dim(deg_table)[2]
    percent_significant_categories<-unique(percent_of_times_significant)
    percent_of_times_sig_df<-cbind(rownames(deg_table), percent_of_times_significant)
    colnames(percent_of_times_sig_df)<-c("gene","percent_deg")
    #print(mean(percent_of_times_sig_df[which(percent_of_times_significant>0)]))
    return(percent_of_times_sig_df)
}

#############


#' @importFrom gprofiler2 gost
get_single_path<-function(temp_degs, temp_bg, species){
    message("inputs:",length(temp_degs))
    print(head(temp_degs))
    message("background:",length(temp_bg))
    print(head(temp_bg))
    if (length(temp_degs)!=0){
        message("    ",species)
        gprof_res<-gost(query = temp_degs, organism = species, custom_bg=temp_bg, significant = TRUE, evcodes=TRUE)#evcodes returns the intersection of genes
        if ( !is.null(gprof_res$result) ){
            gprof_res$result<-gprof_res$result[order(gprof_res$result$p_value),]
            #print(colnames(gprof_res$result))
            print(head(gprof_res$result[,c("p_value","term_name")]))
            temp_paths <-gprof_res$result
            num_paths_sig<-dim(gprof_res$result)[1]
            geom_mean_sig<-exp(mean(log(-log10(unlist(gprof_res$result$p_value)))))
        } else {
            num_paths_sig<-0
            geom_mean_sig<-0
            temp_paths<-NULL
        }
    } else {
        num_paths_sig<-0
        geom_mean_sig<-0
        temp_paths <-NULL
    }
    return(list(num_paths_sig,
                geom_mean_sig,
                temp_paths))
}



merge_up_down<-function(up_res, down_res, up_or_down_res, factor_levels){
    num_paths_sig <- up_res[[1]] + down_res[[1]] + up_or_down_res[[1]]
    geom_mean_sig <- mean(c(up_res[[2]], down_res[[2]], up_or_down_res[[2]]))
    first_lab<-paste("higher_in_",factor_levels[1],sep="")
    second_lab<-paste("lower_in_",factor_levels[1],sep="")
    third_lab<-"all_sig_diffs"
    ## 
    are_null_vect<-c(is.null(up_res[[3]]), is.null(down_res[[3]]), is.null(up_or_down_res[[3]]))
    num_null<-sum(as.integer(are_null_vect))
    # message("num_null:",num_null)
    if (num_null==3){
        message("  no significant pathways")
        ## this means that all are null (no pathways enriched)
        return(list(num_paths_sig,
                    geom_mean_sig,
                    NULL))
    }
    ## need to handle the null result for each one
    non_null_idxs<-which(are_null_vect==FALSE)
    labs <- c(first_lab, second_lab, third_lab)
    ## pull the header from it; we know that there must be at least one entry here
    ## otherwise we would have returned earlier
    all_out_pathway_dfs<-list()
    all_out_pathway_dfs[[1]]<-as.data.frame(up_res[[3]])
    all_out_pathway_dfs[[2]]<-as.data.frame(down_res[[3]])
    all_out_pathway_dfs[[3]]<-as.data.frame(up_or_down_res[[3]])
    temp_idx<-as.integer(non_null_idxs[1])
    # message("index being used: ",temp_idx)
    out_df<-all_out_pathway_dfs[[ temp_idx ]] 
    direction<-rep(labs[non_null_idxs[1]],dim(out_df)[1])
    header<-c("direction",out_df)
    if (length(non_null_idxs)>1){
        ## this should only be entered if there's more than significant of the three
        for (non_null_idx in seq(2,length(non_null_idxs))){
            temp_df<-all_out_pathway_dfs[[non_null_idx]]
            direction<-c(direction,rep(labs[non_null_idxs[non_null_idx]],dim(temp_df)[1]) )
            out_df<-rbind(out_df,temp_df)
        }
    }
    out_df<-as.data.frame(cbind(direction, out_df))
    out_df<-out_df[order(out_df[,"p_value"]),]
    return(list(num_paths_sig,
                geom_mean_sig,
                out_df))
}



#' @importFrom gprofiler2 gost
do_whole_dataset_meta_analysis<-function(deg_table,
                                         avg_logFC_table,
                                         factor_levels,
                                         percent_of_times_sig_df,
                                         species='mmusculus',
                                         dict_with_bg=NULL,
                                         global_bg=NULL){
    ## geg_table: genes x cluster, populated by T/F on whether it's a DEG across batches w/in that cluster
    print(head(deg_table))
    all_pathway_list<-list()
    num_paths_sig_vect<-c()
    geom_mean_sig_vect<-c()
    all_pathway_list[["clusters"]]<-list()
    for (clust in seq(1,dim(deg_table)[2])){
        clust_id<-paste("cluster",clust,sep="_")
        message("working on clust: ",clust)
        temp_degs<-rownames(deg_table)[which(deg_table[,clust]==TRUE)]
        up_genes <- rownames(deg_table)[which(avg_logFC_table[,clust] > 0 )]
        down_genes <- rownames(deg_table)[which(avg_logFC_table[,clust] < 0 )]
        ## double check them for significance
        up_genes <- intersect(up_genes,temp_degs)
        down_genes <- intersect(down_genes,temp_degs)
        #############
        ## get the res
        up_res <- get_single_path(up_genes, dict_with_bg[[clust_id]][["genes_expressed_in_clust"]], species=species)
        down_res <- get_single_path(down_genes, dict_with_bg[[clust_id]][["genes_expressed_in_clust"]], species=species)
        up_or_down_res <- get_single_path(temp_degs, dict_with_bg[[clust_id]][["genes_expressed_in_clust"]], species=species)
        temp_paths_res<-merge_up_down(up_res, down_res, up_or_down_res, factor_levels)
        num_paths_sig<-temp_paths_res[[1]]
        geom_mean_sig<-temp_paths_res[[2]]
        temp_paths<-temp_paths_res[[3]]
        ## log them locally
        message("        num paths sig:",num_paths_sig)
        num_paths_sig_vect<-c(num_paths_sig_vect,num_paths_sig)
        geom_mean_sig_vect<-c(geom_mean_sig_vect,geom_mean_sig)
        all_pathway_list[["clusters"]][[clust_id]]<-temp_paths
    }
    all_pathway_list[["num_paths_sig"]]<-num_paths_sig_vect
    all_pathway_list[["geom_mean_sig"]]<-geom_mean_sig_vect
    ## 
    high_degs<-percent_of_times_sig_df[which(percent_of_times_sig_df[,"percent_deg"]>0.5),"gene"]
    if (length(high_degs)>10){
        gprof_res<-gost(query = high_degs, organism = species, custom_bg=global_bg, significant = TRUE, evcodes=TRUE)#evcodes returns the intersection of genes
        gprof_res$result<-gprof_res$result[order(gprof_res$result$p_value),]
        all_pathway_list[["global_degs"]]<-gprof_res$result
    } else {
        all_pathway_list[["global_degs"]]<-NULL
    }
    return(all_pathway_list)
}



get_null_mat<-function(boots){
    temp_cells<-as.character(boots[[1]][,"node"])
    null_mat<-matrix(rep(0,length(temp_cells)*length(boots)),
                        nrow=length(temp_cells),ncol=length(boots))
    rownames(null_mat)<-temp_cells
    for (i in seq(length(boots))){
        null_mat[,i]<-boots[[i]][,"percent_host_batch"]
    }
    return(null_mat)
}


########################################################################################################

#' run_towcab_analysis
#' @description \code{run_towcab_analysis} Performs the full towcab analysis
#' @param exprs an expression matrix with cells in columns, genes in rows
#' @param batch_vect a vector that can be interpreted as a factor corresponding to batches (same order as columns in exprs)
#' @param clust_vect an integer vector (indexed either from 0 or 1) with cluster ID numbers.
#' @param alpha The BH corrected significance cutoff (default=0.01)
#' @param n_threads number of threads to use (default=10)
#' @param min_express_for_inclusion_percentile the percentile cutoff for determining the downsampling threshold (default is 0.10 (meaning 90% of cells will be included))
#' @param num_pseudo_replicates Number of pseudoreplicates for Libra DE analysis (default=3)
#' @param min_cells_per_pseudo_rep Minimum number of cells to include in each pseudoreplicate (default=3)
#' @param species The species code from gprofiler for pathway analysis (default="mmusculus", could also be something like hsapiens, etc)
#' @param topology_filter Boolean of whether to perform filtering for only topologically overlapping areas.
#' @param in_topology_graph An igraph network - this should be the kNN or similar object that went into the clustering algorithm.
#' @param method_name The name of the method (empty string "" by default). Not required except for printing out some things while the function is running.
#' @return a list of all results from the towcab analysis
#' @importFrom Libra run_de
#' @importFrom stats quantile
#' @importFrom igraph V induced_subgraph
#' @importFrom matrixStats rowMaxs
#' @importFrom Matrix colSums rowSums
#' @name run_towcab_analysis
#' @export
run_towcab_analysis<-function(exprs, 
                              batch_vect,
                              clust_vect,
                              alpha=0.01,
                              n_threads=10,
                              min_express_for_inclusion_percentile=0.10,
                              num_pseudo_replicates=3,
                              min_cells_per_pseudo_rep=3,
                              species='mmusculus',
                              topology_filter=TRUE,
                              in_topology_graph=NULL,
                              method_name=""){
    batch_vect<-factor(batch_vect, levels=unique(as.character(batch_vect)))
    clust_vect<-factor(clust_vect, levels=unique(as.character(clust_vect)))
    results_list<-list()
    counts_per_cell<-Matrix::colSums(exprs)
    min_express_for_inclusion<-as.integer(unlist(quantile(counts_per_cell,c(min_express_for_inclusion_percentile))))
    keep<-which(counts_per_cell>min_express_for_inclusion)
    message("keeping ",length(keep)," out of ",length(counts_per_cell)," cells for DEG analysis")
    exprs<-exprs[,keep]
    keep_cells<-colnames(exprs)
    if (topology_filter){
        graph_nodes<-as.character(unlist(V(in_topology_graph)$name))
        percent_vertex_in_keep_cells<-sum(graph_nodes %in% keep_cells)/length(keep_cells)
        message(percent_vertex_in_keep_cells*100,"% of keep_cells were in the input graph")
        in_topology_graph<-induced_subgraph(in_topology_graph, keep_cells)
        if (percent_vertex_in_keep_cells==0){
            for (blah in seq(10)){
                message()
            }
            message("WARNING: NONE OF THE CELLS WERE FOUND IN THE GRAPH")
            print("head(graph_nodes)")
            print(head(graph_nodes))
            print("head(keep_cells)")
            print(head(keep_cells))
            for (blah in seq(10)){
                message()
            }
        }
    }
    batch_vect<-batch_vect[keep]
    clust_vect<-clust_vect[keep]
    min_clust_vect<-min(as.integer(clust_vect))
    if (min_clust_vect==0){
    	clust_vect<-clust_vect+min_clust_vect+1## plus 1 because R does 1 indexing instead of zero...
    }
    keep_genes<-rownames(exprs)[which(Matrix::rowSums(exprs)>10)]
    exprs<-exprs[keep_genes,]
    #################################
    #################################
    sig_table<-matrix(rep(FALSE,length(unique(clust_vect))*dim(exprs)[1]),
                      ncol=length(unique(clust_vect)),
                      nrow=dim(exprs)[1] )
    p_table<-matrix(rep(1,length(unique(clust_vect))*dim(exprs)[1]),
                      ncol=length(unique(clust_vect)),
                      nrow=dim(exprs)[1] )
    avg_logFC_table<-matrix(rep(0,length(unique(clust_vect))*dim(exprs)[1]),
                      ncol=length(unique(clust_vect)),
                      nrow=dim(exprs)[1] )
    rownames(sig_table)<-rownames(exprs)
    rownames(p_table)<-rownames(exprs)
    rownames(avg_logFC_table)<-rownames(exprs)
    unique_clust_names<-sort(unique(clust_vect))
    colnames(sig_table)<-paste("cluster",unique_clust_names,sep="_")
    colnames(p_table)<-paste("cluster",unique_clust_names,sep="_")
    colnames(avg_logFC_table)<-paste("cluster",unique_clust_names,sep="_")
    ####
    results_list[["global_bg"]]<-rownames(exprs)
    message("num genes in global bg:",length(results_list[["global_bg"]]))
    factor_levels<-sort(unique(batch_vect))
    results_list[["batches_in_logFC_order"]]<-factor_levels
    ## go through each cluster & subset exprs for it
    message("all clusters:")
    print(unique_clust_names)
    if (length(unique_clust_names)>0){
        for (clust in unique_clust_names){#c(unique_clust_names[1])){
            clust_id<-paste("cluster",clust,sep="_")
            message("  working on clust: ",clust_id," ... ",method_name)
            results_list[[clust_id]]<-list()
            ## now downsample the whole cluster
            current_clust_idxs<-which(clust_vect==clust)
            message("     num cells in clust:   ",length(current_clust_idxs))
            message("     current batches:      ",unique(batch_vect[current_clust_idxs]))
            results_list[[clust_id]][["nBatches_per_cluster"]]<-length(unique(batch_vect[current_clust_idxs]))
            message("     nBatches_per_cluster: ",results_list[[clust_id]][["nBatches_per_cluster"]])
            results_list[[clust_id]][["cluster_size"]]<-length(current_clust_idxs)
            message("     cluster size:        ",results_list[[clust_id]][["cluster_size"]])
            ## just setting these as empty, but will fill them in later if the pass the cutoffs
            results_list[[clust_id]][["deg_adj_list"]]<-NULL
            results_list[[clust_id]][["genes_expressed_in_clust"]]<-c()
            results_list[[clust_id]][["deg_coexpression_adj_list"]]<-data.frame(gene_1=c(NA),gene_2=c(NA),spearman_rho=c(NA))
            clust_cell_names<-as.character(colnames(exprs)[current_clust_idxs])
            if (length(unique(batch_vect[current_clust_idxs]))>1){
                ## randomly assign replicates
                replicate<-sample(seq(num_pseudo_replicates),length(current_clust_idxs),replace=TRUE)
                ## looks weird, but Libra does DE across "label" so that's where we'll put the batch ID
                meta<- data.frame(cell_id=clust_cell_names, 
                                  replicate=replicate, 
                                  cell_type=rep(clust_id,length(current_clust_idxs)),
                                  label= as.character(batch_vect)[current_clust_idxs])
                meta[,"label"]<-factor(meta[,"label"],levels=unique(unlist(as.character(meta[,"label"]))))
                # message("meta head:")
                #print(head(meta))
                #print("trying to set rownames for meta")
                cell_table<-table(clust_cell_names)
                if (max(cell_table)>1){
                    print(cell_table[which(cell_table>1)])
                } else {
                    print("no duplicates")
                }
                #print(cbind(make.names(clust_cell_names),clust_cell_names))
                rownames(meta)<-clust_cell_names
                #print("finished setting rownames for meta")
                cross_tab<-table(meta$cell_type,meta$label)
                print(cross_tab)
                ## make sure we have enough cells to do the DE, otherwise it'll throw an error
                second_highest_count<-rev(sort(cross_tab))[2]
                #####################
                ## only include batches that have >9 cells
                min_cells_per_batch <- num_pseudo_replicates*min_cells_per_pseudo_rep
                # message("min_cells_per_batch:",min_cells_per_batch)
                batches_with_enough <- colnames(cross_tab)[which(Matrix::colSums(cross_tab)>(min_cells_per_batch))]
                batches_without_enough <- colnames(cross_tab)[colnames(cross_tab) %!in% batches_with_enough]
                # message("batches with enough:")
                # print(batches_with_enough)
                # message("batches without enough:")
                # print(batches_without_enough)
                meta <- meta[which(meta$label %in% batches_with_enough),]
                keep_cells <- as.character(meta$cell_id)
                ## make sure the cells are lined up
                keep_cells<-intersect(keep_cells, colnames(exprs))
                meta <- meta[keep_cells,]
                #####################
                ## control factor order so that 
                ## a positive log FC means that it's higher in the first one
                ## 
                meta$label<-factor(meta$label, levels=factor_levels)
                message("     second highest count:",second_highest_count)
                if (second_highest_count > min_cells_per_batch){
                    ## make sure the matrix is dense & subset
                    ## now we'll actually run the DEG analysis
                    blah_mat<-as.matrix(exprs[,keep_cells])
                    # message("minimum colSums:")
                    # print(min(colSums(blah_mat)))
                    message("min_express_for_inclusion: ",min_express_for_inclusion)
                    temp_mat <- as.matrix(local_downsample_mat(blah_mat,
                                                         target=as.integer(unlist(min_express_for_inclusion)),
                                                         quiet=TRUE,
                                                         do_parallel=TRUE))
                    results_list[[clust_id]][["genes_expressed_in_clust"]]<-rownames(temp_mat)[which(Matrix::rowSums(temp_mat)>10)]
                    message("      background for clust :",length(results_list[[clust_id]][["genes_expressed_in_clust"]]))
                    offset_list<-offset_reps(meta)
                    meta<-offset_list[["meta"]]
                    offset_lookup<-offset_list[["batch_rep_lookup_table"]]
                    # print(head(meta))
                    # print(names(meta))
                    # print(summary(meta))
                    ## final sanity check
                    keep_cells<-intersect(as.character(meta$cell_id),colnames(temp_mat))
                    meta<-meta[keep_cells,]
                    temp_mat<-temp_mat[,keep_cells]
                    meta<-make_sure_enough_cells_in_each_rep(meta, num_pseudo_replicates=num_pseudo_replicates, min_cells_per_pseudo_rep=min_cells_per_pseudo_rep, offset_lookup=offset_lookup)
                    print(table(meta$replicate,meta$label))
                    DE <- NULL
                    sig_gene_idxs <- NULL
                    sig_genes <- NULL
                    if (length(unique(meta$label))==2){
                        ######################################
                        ######################################
                        ######################################
                        enough_cells<-TRUE
                        if (topology_filter){
                            temp_mat_meta<-do_topology_filter(temp_mat, meta, in_topology_graph)
                            temp_mat<-temp_mat_meta[["temp_mat"]]
                            meta<-temp_mat_meta[["meta"]]
                            if (min(table(meta$label))>min_cells_per_batch){
                                meta<-make_sure_enough_cells_in_each_rep(meta, num_pseudo_replicates=num_pseudo_replicates, min_cells_per_pseudo_rep=min_cells_per_pseudo_rep, offset_lookup=offset_lookup)
                                keep_cells<-as.character(meta$cell_id)
                                enough_cells<-TRUE
                            } else{
                                enough_cells<-FALSE
                            }
                            #print(meta)
                            #print(keep_cells)
                        }
                        ######################################
                        ######################################
                        ######################################
                        if (enough_cells){
                            print(table(meta$replicate,meta$label))
                            DE <- run_de(temp_mat, 
                                        meta = meta, 
                                        replicate_col = "replicate",
                                        cell_type_col = "cell_type",
                                        de_family = 'pseudobulk', 
                                        de_type = 'LRT',
                                        n_threads = n_threads)
                            sig_gene_idxs<-which(DE$p_val_adj<alpha)
                            sig_genes<-DE$gene[sig_gene_idxs]
                            message("    ",method_name," -- ",length(sig_gene_idxs)," significant genes; ","    mean logFC: ",round(mean(DE$avg_logFC[sig_gene_idxs]),2)," +/- ",round(sd(DE$avg_logFC[sig_gene_idxs])))
                            # print(DE[1:5,"p_val_adj"])
                            DE<-as.data.frame(DE)
                            #print(head(DE))
                            #print(dim(DE))
                            #print("head(colnames(sig_table))")
                            #print(head(colnames(sig_table)))
                            #print("head(rownames(sig_table))")
                            #print(head(rownames(sig_table)))
                            for (temp_idx in seq(dim(DE)[1])){
                                temp_gene<-DE$gene[temp_idx]
                                ## wrote this just in case the genes come out in different orders
                                ## NOTE: They were in different orders...
                                # message("temp_gene: ",temp_gene)
                                # message("    ",as.logical(DE[temp_idx,"p_val_adj"]<alpha)," : ",DE[temp_idx,"p_val_adj"])
                                sig_table[temp_gene,clust_id]<-as.logical(DE[temp_idx,"p_val_adj"]<alpha)
                                p_table[temp_gene,clust_id]<-DE[temp_idx,"p_val_adj"]
                                avg_logFC_table[temp_gene,clust_id]<-DE[temp_idx,"avg_logFC"]
                            }
                        } else {
                            message ("      not quite enough mixing to do DE after filtering on topology")
                        }
                    } else {
                        ###############################################################
                        ###############################################################
                        ## run_de appears to only work when there are only 2 categories
                        unique_batches<-sort(unique(as.character(meta$label)))
                        num_combos<-as.integer((length(unique_batches)**2 - length(unique_batches))/2)
                        temp_sig_table<-matrix(rep(0,dim(temp_mat)[1]*num_combos),ncol=num_combos)
                        temp_p_table<-matrix(rep(1,dim(temp_mat)[1]*num_combos), ncol=num_combos)
                        rownames(temp_sig_table)<-rownames(temp_mat)
                        rownames(temp_p_table)<-rownames(temp_mat)
                        counter<-0
                        for (i in seq(length(unique_batches))){
                            temp_b1<-unique_batches[i]
                            for (j in seq(i,length(unique_batches))){
                                temp_b2<-unique_batches[j]
                                if (i!=j){
                                    counter<-counter+1
                                    temp_run_id<-paste(temp_b1,temp_b2,sep="_vs_")
                                    # print("old meta")
                                    # print(head(meta))
                                    cell_subset_id_idxs<-which(as.character(meta$label)==temp_b1 | as.character(meta$label)==temp_b2)
                                    # print("cell_subset_id_idxs")
                                    # print(cell_subset_id_idxs)
                                    temp_meta2<-meta[cell_subset_id_idxs,]
                                    temp_meta2[,"label"]<-factor(temp_meta2[,"label"],
                                                                 levels=sort(c(temp_b1,temp_b2)) )
                                    temp_meta2[,"replicate"]<-factor(temp_meta2[,"replicate"],
                                                                levels = unique(as.character(temp_meta2[,"replicate"])) )
                                    # print("temp_meta2")
                                    # print(temp_meta2)
                                    #print(table(temp_meta2$replicate,temp_meta2$label))
                                    temp_mat2<-temp_mat[,as.character(temp_meta2$cell_id)]
                                    # print(dim(temp_mat2))
                                    ######################################
                                    ######################################
                                    ######################################
                                    enough_cells<-TRUE
                                    if (topology_filter){
                                        temp_mat_meta2<-do_topology_filter(temp_mat2, temp_meta2, in_topology_graph)
                                        temp_mat2<-temp_mat_meta2[["temp_mat"]]
                                        temp_meta2<-temp_mat_meta2[["meta"]]
                                        if (min(table(temp_meta2$label))>min_cells_per_batch){
                                            temp_meta2<-make_sure_enough_cells_in_each_rep(temp_meta2, num_pseudo_replicates=num_pseudo_replicates, min_cells_per_pseudo_rep=min_cells_per_pseudo_rep, offset_lookup=offset_lookup)
                                            #keep_cells<-as.character(meta$cell_id)
                                            enough_cells<-TRUE
                                        } else{
                                            enough_cells<-FALSE
                                        }
                                    }
                                    ######################################
                                    ######################################
                                    ######################################
                                    if (enough_cells){
                                        DE <- run_de(temp_mat2, 
                                                        meta = temp_meta2, 
                                                        replicate_col = "replicate",
                                                        cell_type_col = "cell_type",
                                                        de_family = 'pseudobulk', 
                                                        de_type = 'LRT',
                                                        min_cells=1,
                                                        min_reps=1,
                                                        n_threads = n_threads)
                                        # print(head(DE))
                                        sig_gene_idxs<-which(DE$p_val_adj<alpha)
                                        sig_genes<-DE$gene[sig_gene_idxs]
                                        message("    ",length(sig_gene_idxs)," significant genes; ","    mean logFC: ",round(mean(DE$avg_logFC[sig_gene_idxs]),2)," +/- ",round(sd(DE$avg_logFC[sig_gene_idxs])))
                                        # print(DE[1:5,"p_val_adj"])
                                        DE<-as.data.frame(DE)
                                        for (temp_idx in seq(dim(DE)[1])){
                                            #message(as.logical(DE[temp_idx,"p_val_adj"]<alpha),DE[temp_idx,"p_val_adj"])
                                            temp_gene<-DE$gene[temp_idx]
                                            ## wrote this just in case the genes come out in different orders
                                            temp_sig_table[temp_gene,counter]<-as.integer(DE[temp_idx,"p_val_adj"]<alpha)
                                            temp_p_table[temp_gene,counter]<-DE[temp_idx,"p_val_adj"]
                                            #temp_avg_logFC_table[temp_gene,counter]<-DE[temp_idx,"avg_logFC"]
                                        }
                                    } else {
                                        message ("      not quite enough mixing to do DE after filtering on topology")
                                    }
                                }
                            }
                        }

                        sig_table[,clust_id]<-as.logical(rowMaxs(temp_sig_table))
                        p_table[,clust_id]<-rowMins(temp_p_table)
                        avg_logFC_table[,clust_id] <-rep(NA,dim(temp_mat)[1])
                        ###############################################################
                        ###############################################################
                    }
                } else {
                    message ("      not quite enough mixing to do DE")
                }
                #############################
                #############################
            } else {
                message("     Either only one batch, or not quite enough mixing to do DE")
            }
            #print(results_list)
        }
    } else {
        stop("nothing found?")
    }
    #print(results_list)
    percent_sig<-get_pcnt_deg_df(sig_table)
    pathway_results<-do_whole_dataset_meta_analysis(sig_table,
                                                    avg_logFC_table,
                                                    factor_levels,
                                                    percent_sig,
                                                    species=species,
                                                    dict_with_bg=results_list, 
                                                    global_bg=results_list[["global_bg"]])
    out_list<-list(sig_table=sig_table,
                p_table=p_table,
                avg_logFC=avg_logFC_table,
                percent_sig=percent_sig,
                num_degs=colSums(sig_table),
                pathway_results=pathway_results,
                full_stats_list=results_list)
    return(out_list)
}






########################################################################


analyze_all_pathway_list_results<-function(top_deg_lists, input_lists, top_out_dir, method_colors, organism, technical_file = "data/technical_studies/all_technical_pathways.tsv"){
    all_results_list<-list()
    for (run in names(top_deg_lists)){
        temp_out_dir<-paste(c(top_out_dir, run), collapse="/")#, "results","pathway_DEGs/"
        dir.create(temp_out_dir,showWarnings=FALSE)
        all_results_list[[run]]<-analyze_all_pathway_results(top_deg_lists[[run]], temp_out_dir, method_colors, organism=organism[[run]][["SPECIES"]], technical_file = technical_file)
    }
    return(all_results_list)
}


#' analyze_all_pathway_results
#' @description \code{analyze_all_pathway_results} Performs the full towcab method comparison and high level overview
#' @param sig_res a list whose names correspond to the names of the methods, coding the towcab results in that entry
#' @param out_dir the output directory to write all the results into
#' @param method_colors a list with the method names and the colors to code them by
#' @param organism takes the gprofiler codes ('hsapiens, mmusculus, etc')
#' @param technical_file a file containing the names of pathways that should be filtered out & not \
#'                       counted against a method. These should be defined as technical in nature, \
#'                       and need to have evidence through technical replicates that these are technical. \
#'                       There is a default file distributed with the package that has a draft generated \
#'                       from technical replicates of mouse brain and heart with Chromium V3.
#' @return writes outputs in the out_dir
#' @importFrom Libra run_de
#' @importFrom stats quantile
#' @importFrom igraph V
#' @importFrom ggplot2 aes geom_boxplot scale_fill_manual
#' @importFrom cowplot plot_grid
#' @name analyze_all_pathway_results
#' @export
analyze_all_pathway_results<-function(deg_lists, out_dir, method_colors, organism="mmusculus", technical_file=NULL){
    if (is.null(technical_file)){
        technical_file<- paste0(path.package("towcab"), "/all_technical_pathways.tsv")
        ## technical_file<-paste(dirname(sys.frame(1)$ofile),"/all_technical_pathways.tsv",sep="")
    }
    dir.create(out_dir,showWarnings=FALSE)
    res_dir<-paste(out_dir,"results",sep='/')
    dir.create(res_dir,showWarnings=FALSE)
    out_res_file<-paste(out_dir,"results/pathway_DEGs/all_deg_pathway_results.Rds",sep="/")
    if (!file.exists(out_res_file)){
        temp_path_analysis_results<-list()
        plot_winClust_acrossBatch_DEG_results(deg_lists, out_dir, method_colors)
        ############
        all_technical_pathways<-as.character(unlist(read.csv(technical_file)))
        deg_pathway_table<-write_combined_significance_table(deg_lists,
                                                             paste(out_dir,"results/",sep="/"),
                                                             technical_file=technical_file)
        ############
        non_technical_table<-deg_pathway_table[which(deg_pathway_table[,"term_id"] %!in% all_technical_pathways),]
        num_non_technical <- length(unique(non_technical_table[,"term_name"]))
        #all_intestine_pathways<-sort(unique(intestine_different_state_pathway_results[,"term_id"]))
        non_technical_pathways<-analyze_DEG_pathways(non_technical_table, out_dir, organism=organism)
        path_sig_res<-quantify_non_technical_pathways(deg_pathway_table,
                                        non_technical_pathways,
                                        out_dir,
                                        names(deg_lists),
                                        technical_file=technical_file)
        out_file<-paste(out_dir,"results/pathway_DEGs/all_significant_pathways_collapsed.tsv",sep="/")
        write.table(path_sig_res, file=out_file, quote=FALSE, row.names=FALSE, sep="\t")
        temp_path_analysis_results[["deg_pathway_table"]]<-deg_pathway_table
        temp_path_analysis_results[["path_sig_res"]]<-path_sig_res
        ############
        saveRDS(temp_path_analysis_results, file=out_res_file)
    } else {
        temp_path_analysis_results<-readRDS(out_res_file)
    }
    return(temp_path_analysis_results)
}

#########




plot_winClust_acrossBatch_DEG_results_lists<-function(deg_results,out_dir){
    for (run in names(deg_results)){
        temp_sub_out_dir<-paste(out_dir,run,sep="/")
        plot_winClust_acrossBatch_DEG_results(deg_results[[run]],temp_sub_out_dir)
        all_deg_res<-deg_results[[run]]
        save(all_deg_res,file=paste(temp_sub_out_dir,"results/all_deg_res.Rds",sep="/"))
    }
}



get_num_sig_df<-function(sig_res){
    meth_vect<-c()
    num_paths_sig<-c()
    num_degs<-c()
    geom_mean_sig<-c()
    for (method in names(sig_res)){
        #print(head(sig_res[[method]]$pathway_result$num_paths_sig))
        num_paths_sig <-c(num_paths_sig,sig_res[[method]]$pathway_results$num_paths_sig)
        num_degs<-c(num_degs,sig_res[[method]]$num_degs)
        meth_vect<-c(meth_vect,rep(method,length(sig_res[[method]]$num_degs)))
        geom_mean_sig<-c(geom_mean_sig, sig_res[[method]]$pathway_results$geom_mean_sig)
    }
    return(data.frame(method=meth_vect,
                      num_paths_sig=num_paths_sig,
                      num_degs=num_degs,
                      geometric_mean_of_significance=geom_mean_sig))
}




#' @importFrom data.table fwrite
write_combined_significance_table<-function(sig_res, out_dir, technical_file="data/technical_studies/all_technical_pathways.tsv"){
    all_technical_pathways<-as.character(unlist(read.csv(technical_file)))
    final_out_table<-NULL
    for (method in names(sig_res)){
        for (temp_cluster in names(sig_res[[method]]$pathway_results$clusters)){
            temp_sub_table<-sig_res[[method]]$pathway_results$clusters[[temp_cluster]]
            correction_method<-rep(method, dim(temp_sub_table)[1])
            cluster_id<-rep(temp_cluster, dim(temp_sub_table)[1])
            temp_sub_table<-as.data.frame(cbind(correction_method,cluster_id, temp_sub_table))
            if (is.null(final_out_table)){
                final_out_table<-temp_sub_table
            } else {
                final_out_table<-rbind(final_out_table, temp_sub_table)
            }
        }
    }
    in_technical_paths<-rep(TRUE,dim(final_out_table)[1])
    in_technical_paths[which(final_out_table[,"term_id"] %!in% all_technical_pathways)]<-FALSE
    final_out_table<-as.data.frame(cbind(final_out_table,in_technical_paths))
    print("head(final_out_table)")
    print(head(final_out_table))
    out_file<-paste(out_dir,"all_significant_pathways_raw.tsv",sep="/")
    #write.table(final_out_table, file=out_file, quote=FALSE, row.names=FALSE, sep="\t")
    data.table::fwrite(final_out_table, file=out_file, quote=FALSE, row.names=FALSE, sep="\t")
    return(final_out_table)
}
########



########################################################################














#####################################################################################


#' @importFrom stringr str_split
get_single_comparison<-function(line_1, line_2){
    gene_list_1<-unlist(str_split(unlist(line_1["intersection"]),","))
    gene_list_2<-unlist(str_split(unlist(line_2["intersection"]),","))
    return(length(intersect(gene_list_1,gene_list_2))/min(c(length(gene_list_1),length(gene_list_2))))
}


get_par_comparisons<-function(line_num, comparison_idx_df, in_df){
    line_1<-in_df[comparison_idx_df[line_num,1],]
    line_2<-in_df[comparison_idx_df[line_num,2],]
    term_1<-in_df[comparison_idx_df[line_num,1],"term_id"]
    term_2<-in_df[comparison_idx_df[line_num,2],"term_id"]
    gene_list_1<-unlist(str_split(unlist(line_1["intersection"]),","))
    gene_list_2<-unlist(str_split(unlist(line_2["intersection"]),","))
    overlap<-length(intersect(gene_list_1,gene_list_2))/min(c(length(gene_list_1),length(gene_list_2)))
    return(data.frame(path_1=c(as.character(term_1)),
                      path_2=c(as.character(term_2)),
                      overlap=overlap)
                    )
}


make_comparison_idx_df<-function(num_verts){
    idx_1<-c()
    idx_2<-c()
    for (i in seq(num_verts-1)){
        idx_1<-c(idx_1,rep(i,num_verts-i))
        idx_2<-c(idx_2,seq(i+1,num_verts))
    }
    return(data.frame(idx_1,idx_2))
}


do_single_collate<-function(X, temp_subset){
    return(do.call(rbind.data.frame, temp_subset[X]))
}


#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster detectCores clusterExport stopCluster
par_collate<-function(adj_list_list){
    num_rows<-length(names(adj_list_list))
    num_cores<-detectCores()-1
    cl <- makeCluster(num_cores)
    rows_per_core<-as.integer(num_rows/num_cores)
    ####################
    ## get the rows allocated
    row_lists<-list()
    for (i in seq(num_cores-1)){
        temp_start_idx<-((i-1)*rows_per_core)+1
        row_lists[[i]]<-seq(temp_start_idx, temp_start_idx+rows_per_core)
    }
    ## append the last
    temp_start_idx<-max(row_lists[[i]])+1
    row_lists[[i+1]]<-seq(temp_start_idx,num_rows)
    ####################
    adj_collated_lists<-pblapply(X=row_lists,
                    FUN=do_single_collate, 
                    cl=cl,
                    temp_subset=adj_list_list)
    stopCluster(cl)
    print(head(adj_collated_lists[[1]]))
    adj_list_df<-do.call(rbind.data.frame, adj_collated_lists)
    return(adj_list_df)
}


#' @importFrom pbapply pblapply
#' @importFrom parallel makeCluster detectCores clusterExport stopCluster
get_pathway_network_parallel<-function(in_df){
    #in_df<-in_df[sample(seq(dim(in_df)[1]),1500),]
    comparison_idx_df<-make_comparison_idx_df(dim(in_df)[1])
    cl <- makeCluster(detectCores())
    clusterExport(cl, c("str_split"), envir=environment())
    adj_list_list<-pblapply(X=seq(dim(comparison_idx_df)[1]),
                    FUN=get_par_comparisons, 
                    cl=cl,
                    comparison_idx_df=comparison_idx_df,
                    in_df=in_df)
    stopCluster(cl)
    print("collating adj list in parallel")
    #adj_list_df<-par_collate(adj_list_list)
    adj_list_df<-do.call(rbind.data.frame, adj_list_list)
    print(head(adj_list_df))
    return(adj_list_df)
}


get_pathway_network<-function(in_df){
    message("getting adj list:")
    path_1<-c()
    path_2<-c()
    overlap<-c()
    for (i in seq(1,dim(in_df)[1])){
        if (i %% 50==0){
            message("    ",i/dim(in_df)[1])
        }
        line_1<-in_df[i,]
        term_1<-as.character(unlist(line_1["term_id"]))
        message("        ",term_1)
        for (j in seq(i,dim(in_df)[1])){
            if (i!=j){
                line_2<-in_df[j,]
                term_2<-as.character(unlist(line_2["term_id"]))
                temp_overlap<-get_single_comparison(line_1, line_2)
                if (temp_overlap>0){
                    path_1<-c(path_1,term_1)
                    path_2<-c(path_2,term_2)
                    overlap<-c(overlap,temp_overlap)
                }
            }
        }
    }
    return(data.frame(path_1, path_2, overlap))
}


pre_process_pathway_df<-function(in_df){
    all_paths<-sort(unique(in_df[,"term_id"]))
    term_id<-c()
    term_name<-c()
    intersection<-c()
    for (temp_path in all_paths){
        temp_genes<-c()
        all_intersections_list<-in_df[which(in_df[,"term_id"]==temp_path),"intersection"]
        term_name<-c(term_name,in_df[which(in_df[,"term_id"]==temp_path)[1],"term_name"])
        for (temp_list in all_intersections_list){
            temp_genes <- c(temp_genes, unlist(str_split(temp_list, ",")))
        }
        term_id <- c(term_id, temp_path)
        #print(temp_path)
        # print(temp_genes)
        temp_gene_str <- paste(sort(unique(temp_genes)),collapse=",")
        #print(temp_gene_str)
        intersection<-c(intersection,temp_gene_str)
    }
    return(data.frame(term_id, term_name, intersection))
}


#' @importFrom grDevices png
#' @importFrom igraph graph_from_data_frame cluster_louvain get.edgelist vcount
#' @importFrom qgraph qgraph.layout.fruchtermanreingold
analyze_DEG_pathways<-function(non_technical_table, out_dir, organism="mmusculus"){
    dir.create(paste(out_dir,"results/pathway_DEGs/",sep="/"),showWarnings=FALSE)
    png(paste(out_dir,"results/pathway_DEGs/WInClust_AcrossBatch_DEG_pathway_meta_analysis_%02d.png",sep="/"),
        width=5,height=4,units="in",res=450)
    vertex_table<-pre_process_pathway_df(non_technical_table)
    non_technical_adj<-get_pathway_network_parallel(vertex_table)
    #non_technical_adj<-get_pathway_network(vertex_table[sample(seq(dim(vertex_table)[1]),500),])
    hist(non_technical_adj$overlap)
    filtered_non_technical_adj<-non_technical_adj[which(non_technical_adj$overlap>0.5),]
    colnames(filtered_non_technical_adj)<-c("term_1","term_2","weight")
    g <- graph_from_data_frame(filtered_non_technical_adj, directed=FALSE, vertices=vertex_table)
    clust_results <- cluster_louvain(g, weights = NA)
    clust_result_df <- as.data.frame(cbind(clust_results$names,clust_results$membership), stringsAsFactors=F)
    colnames(clust_result_df) <- c("term_id","cluster")
    clust_result_df <- merge(clust_result_df, vertex_table, by="term_id")
    nverts<-dim(vertex_table)[1]
    # mins<-rep(0,nverts)
    # maxs<-rep(10,nverts)
    #temp_layout <- layout.fruchterman.reingold(g, minx=mins,maxx=maxs,miny=mins,maxy=maxs, niter=2000)
    #################
    nclust<-length(unique(clust_results$membership))
    color_pal <- get_set_colors()[seq(nclust)]#rainbow(nclust, alpha=0.6)
    #temp_layout <- layout.fruchterman.reingold(g, niter=100)
    e <- get.edgelist(g,names=FALSE)
    temp_layout <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g))
    plot(g, 
         layout=temp_layout, 
         vertex.size=3, 
         vertex.label=NA, 
         vertex.color=color_pal[clust_results$membership],
         vertex.frame.width=0,#rep(100,nverts),
         edge.color=rgb(0,0,0,alpha=.05))
    plot(g, 
         layout=temp_layout, 
         vertex.size=3, 
         vertex.label=clust_results$membership, 
         #vertex.label.size=100, 
         vertex.color=color_pal[clust_results$membership],
         vertex.frame.width=rep(100,nverts),
         edge.color=rgb(0,0,0,alpha=.05))
    plot(rep(0,nclust),
         seq(nclust),
         col="black",
         bg=color_pal, 
         pch=21, 
         cex=2,
         ylim=c(-1,nclust+1))
    #################
    clust_result_df<-clust_result_df[order(clust_result_df$cluster),]
    reorganized_cluster_result_df<-NULL
    message("reorganizing results based on clusters")
    clust_sizes<-table(as.character(clust_result_df$cluster))
    # print(clust_sizes)
    # print(max(clust_sizes))
    # big_clusts<-sort(colnames(clust_sizes)[which(clust_sizes>1)])
    # print(big_clusts)
    # clust_result_df<-clust_result_df[,]
    for (temp_clust in sort(unique(clust_result_df$cluster)) ){
        for(i in seq(10)){
            message("")
        }
        message(temp_clust)
        for(i in seq(3)){
            message("")
        }
        ## get the unified genes and annotate them
        #temp_gene_str<-paste(unlist(clust_result_df[which(clust_result_df$cluster==temp_clust),"intersection"]),collapse=",")
        #temp_gene_str<-sort(unique(unlist(str_split(temp_gene_str,","))))
        #conversion_res<-gconvert(temp_gene_str, organism=organism)
        #print(colnames(conversion_res))
        # temp_gene_table<-conversion_res[,c("target", "name","description")]
        # colnames(temp_gene_table)<-c("gene","name","description")
        # pathway_cluster<-rep(temp_clust,dim(temp_gene_table)[1])
        # temp_gene_table<-as.data.frame(cbind(pathway_cluster,temp_gene_table))
        #print(temp_gene_table)
        ## get the subset dataframe for the reorganization
        subset_df<-clust_result_df[which(clust_result_df$cluster==temp_clust),]
        subset_df<-subset_df[order(subset_df$term_name),]
        if (is.null(reorganized_cluster_result_df)){
            reorganized_cluster_result_df<-subset_df
        } else {
            reorganized_cluster_result_df<-as.data.frame(rbind(reorganized_cluster_result_df,subset_df))
        }
        #print(as.character(clust_result_df[which(clust_result_df$cluster==temp_clust),"term_name"]))
        # print(temp_gene_str)
    }
    write.table(apply(non_technical_table,2,as.character), file=paste(out_dir,"results/pathway_DEGs/DEG_pathway_non_technical_full_annotations.tsv",sep="/"), quote=FALSE, sep="\t", row.names=F)
    write.table(apply(filtered_non_technical_adj,2,as.character), paste(out_dir,"results/pathway_DEGs/DEG_pathway_non_technical_adj_list.tsv",sep="/"), quote=FALSE, sep="\t")
    write.table(apply(reorganized_cluster_result_df,2,as.character), paste(out_dir,"results/pathway_DEGs/DEG_pathway_clusters.tsv",sep="/"), quote=FALSE, sep="\t")
    sink(file=paste(out_dir,"results/DEG_pathway_descriptions.txt",sep="/"))
    for (path_clust in sort(unique(reorganized_cluster_result_df[,c("cluster")]))){
        for (i in seq(5)){
            print("")
        }
        print(path_clust)
        temp_subset<-reorganized_cluster_result_df[which(reorganized_cluster_result_df[,c("cluster")]==path_clust),c("term_id","term_name")]
        for (i in seq(dim(temp_subset)[1])){
            print(paste0(temp_subset[i,c("term_id")],":",temp_subset[i,c("term_name")]))
        }
    }
    sink()
    dev.off()
    return(reorganized_cluster_result_df)
}

#' @importFrom grDevices png
#' @importFrom ggplot2 ggplot aes geom_bar theme labs element_text
plot_conservation_of_tm_tr<-function(unified_summary, out_dir){
    png(paste(out_dir,"/Tm_Tr_conservation_%02d.png",sep="/"),
        width=5,height=4,units="in",res=450)
    tm_overlap<-unlist(unified_summary[unified_summary$effect_class=="Tm","relative_overlap"])
    tr_overlap<-unlist(unified_summary[unified_summary$effect_class=="Tr","relative_overlap"])
    t_res<-t.test(tm_overlap,mu=tr_overlap)
    tm_v_tr_title<-paste0("T=",round(t_res$statistic,2),"\n","P=",format(t_res$p.value,scientific=T, digits=3))
    p<-ggplot(unified_summary,
           aes(x=term, y=relative_overlap, fill=effect_class) )
    p<-p+geom_bar(stat="identity") 
    p<-p+theme(text=element_text(size=14), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    p<-p+labs(title=tm_v_tr_title)
    print(p)
    dev.off()
}


#' @importFrom RColorBrewer brewer.pal
get_set_colors<-function(){
    #col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    col_vector <- c(brewer.pal(n = 9, name = c("Set1")),brewer.pal(n = 8, name = c("Dark2")))
    return(col_vector)
}

###################################################

get_path_to_path_clust_lookup<-function(non_technical_pathways){
    lookup<-list()
    for (i in seq(dim(non_technical_pathways)[1])){
        path=non_technical_pathways[i,"term_id"]
        lookup[[path]]=as.character(non_technical_pathways[i,"cluster"])
    }
    return(lookup)
}


get_min_p_path_table<-function(non_tech_subset, pathway_lookup, all_methods){
    correction_method<-c()
    pathway<-c()
    pathway_cluster<-c()
    p_val<-c()
    for (temp_meth in all_methods){
        subset_df<-non_tech_subset[which(non_tech_subset[,"correction_method"]==temp_meth),]
        message(temp_meth)
        print(dim(subset_df))
        for (temp_path in unique(subset_df[,"term_id"])){
            method_temp_path_subset_df<-subset_df[which(subset_df[,"term_id"]==temp_path),]
            temp_method_path_min_p<-min(method_temp_path_subset_df[,"p_value"])
            correction_method<-c(correction_method, temp_meth)
            pathway<-c(pathway, temp_path)
            pathway_cluster<-c(pathway_cluster,pathway_lookup[[temp_path]])
            p_val<-c(p_val, temp_method_path_min_p)
        }
    }
    pathway_cluster<-factor(pathway_cluster, levels=sort(unique(pathway_cluster)))
    correction_method<-factor(correction_method, levels=all_methods)
    out_df <- data.frame(correction_method,
                      pathway,
                      pathway_cluster,
                      p_val)
    return(out_df)
}


collapse_min_p<-function(min_p_path_table, all_methods){
    min_p_path_table
    correction_method<-c()
    pathway_cluster<-c()
    num_sig_pathways<-c()
    percent_sig_pathways<-c()
    all_clusters<-sort(unique(min_p_path_table$pathway_cluster))
    for (temp_clust in all_clusters){
        clust_subset<-min_p_path_table[which(min_p_path_table$pathway_cluster==temp_clust),]
        num_paths_per_clust<-length(unique(clust_subset[,"pathway"]))
        for (temp_method in all_methods){
            clust_meth_subset<-clust_subset[which(clust_subset$correction_method==temp_method),]
            num_paths_per_meth_clust<-dim(clust_meth_subset)[1]
            correction_method<-c(correction_method, temp_method)
            pathway_cluster<-c(pathway_cluster, as.character(temp_clust))
            num_sig_pathways<-c(num_sig_pathways,num_paths_per_meth_clust)
            percent_sig_pathways<-c(percent_sig_pathways,num_paths_per_meth_clust/num_paths_per_clust)
        }
    }
    correction_method<-factor(correction_method,levels=all_methods)
    pathway_cluster<-factor(pathway_cluster,levels=all_clusters)
    return(data.frame(correction_method,
                      pathway_cluster,
                      num_sig_pathways,
                      percent_sig_pathways))
}


#' @importFrom grDevices png
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot aes geom_bar theme labs element_text scale_fill_manual geom_boxplot
quantify_non_technical_pathways<-function(deg_pathway_table, non_technical_pathways, out_dir, all_methods, technical_file="data/technical_studies/all_technical_pathways.tsv"){
    dir.create(out_dir, showWarnings=F)
    png(paste(out_dir,"/results/pathway_DEGs/NonTechnical_WInClust_AcrossBatch_DEG_analysis_%02d.png",sep="/"),width=11,height=4,units="in",res=450)
    all_technical_paths<-as.character(unlist(read.csv(technical_file, stringsAsFactors=F)))
    all_non_technical_path_ids<-as.character(unlist(non_technical_pathways[,"term_id"]))
    nclust<-length(unique(non_technical_pathways$cluster))
    print(nclust)
    color_pal <- get_set_colors()[seq(nclust)]
    print(color_pal)
    non_tech_subset<-deg_pathway_table[deg_pathway_table[,"term_id"] %in% all_non_technical_path_ids,]
    print("get_path_to_path_clust_lookup")
    pathway_to_clust_lookup<-get_path_to_path_clust_lookup(non_technical_pathways)
    print("get_min_p_path_table")
    min_p_path_table<-get_min_p_path_table(non_tech_subset, pathway_to_clust_lookup, all_methods)
    print(ggplot(min_p_path_table, aes(x=correction_method, y=-log10(p_val),fill=pathway_cluster)) + geom_boxplot()) + scale_fill_manual(values=color_pal)
    print("collapse_min_p")
    path_sig_df<-collapse_min_p(min_p_path_table, all_methods)
    print(head(path_sig_df))
    p_percent<-ggplot(path_sig_df, aes(x=correction_method, y=percent_sig_pathways,fill=pathway_cluster)) + geom_bar(stat="identity")
    p_percent<-p_percent + theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) + scale_fill_manual(values=color_pal)
    print(p_percent)
    p_num<-ggplot(path_sig_df, aes(x=correction_method, y=num_sig_pathways,fill=pathway_cluster)) + geom_bar(stat="identity")
    p_num<-p_num+theme(text=element_text(size=14), axis.text.x = element_text(angle = 20, vjust = 1, hjust=1)) + scale_fill_manual(values=color_pal)
    print(p_num)
    plot_grid(p_num,p_percent,nrow=2)
    dev.off()
    return(path_sig_df)
}

















########################################################################################


















##############################################################################################
## plot_win_clust_acrossBatch_DEG_results and support
get_num_sig_df<-function(sig_res){
    meth_vect<-c()
    num_paths_sig<-c()
    num_degs<-c()
    geom_mean_sig<-c()
    for (method in names(sig_res)){
        #print(head(sig_res[[method]]$pathway_result$num_paths_sig))
        num_paths_sig <-c(num_paths_sig,sig_res[[method]]$pathway_results$num_paths_sig)
        num_degs<-c(num_degs,sig_res[[method]]$num_degs)
        meth_vect<-c(meth_vect,rep(method,length(sig_res[[method]]$num_degs)))
        geom_mean_sig<-c(geom_mean_sig, sig_res[[method]]$pathway_results$geom_mean_sig)
    }
    return(data.frame(method=meth_vect,
                      num_paths_sig=num_paths_sig,
                      num_degs=num_degs,
                      geometric_mean_of_significance=geom_mean_sig))
}


#' run_towcab_analysis
#' @description \code{run_towcab_analysis} Performs the full towcab analysis
#' @param sig_res an expression matrix with cells in columns, genes in rows
#' @param out_dir a vector that can be interpreted as a factor corresponding to batches (same order as columns in exprs)
#' @param method_colors TODO
#' @return a list of all results from the towcab analysis
#' @examples
#'    x <- rep( c( rep(1, 50), rep(2, 50) ), 2)
#'    y <- c(rep(1, 100), rep(2, 100))
#'    cont_table <- get_cont_table(x, y)
#' @importFrom Libra run_de
#' @importFrom stats quantile
#' @importFrom igraph V
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_manual
#' @importFrom cowplot plot_grid
#' @importFrom grDevices png
#' @name plot_winClust_acrossBatch_DEG_results
#' 
plot_winClust_acrossBatch_DEG_results<-function(sig_res,out_dir, method_colors=NULL){
    out_dir<-paste(out_dir,"results",sep="/")
    png(paste(out_dir,"WInClust_AcrossBatch_DEG_analysis_%02d.png",sep="/"),width=5,height=4,units="in",res=450)
    num_sig_df<-get_num_sig_df(sig_res)
    print(sort(unique(num_sig_df$method)))
    if (!is.null(method_colors)){
    	num_sig_df$method<-factor(num_sig_df$method, levels=names(method_colors))
    }
    print(sort(unique(num_sig_df$method)))
    print(num_sig_df)
    p1<-ggplot(num_sig_df, aes(x=method, y=log(num_paths_sig+1), fill=method)) + geom_boxplot() + scale_fill_manual(values=c(method_colors))
    p3<-ggplot(num_sig_df, aes(x=method, y=log(num_degs+1), fill=method)) + geom_boxplot() + scale_fill_manual(values=c(method_colors))
    p4<-ggplot(num_sig_df, aes(x=method, y=geometric_mean_of_significance, fill=method)) + geom_boxplot() + scale_fill_manual(values=c(method_colors))
    if (!is.null(method_colors)){
    	p1<- p1 + scale_fill_manual(values=c(method_colors))
    	p3<- p3 + scale_fill_manual(values=c(method_colors))
    	p4<- p4 + scale_fill_manual(values=c(method_colors))
    }
    print(p1)
    print(p3)
    print(p4)
    dev.off()
    outfile <- paste(out_dir, "combined_DEG_analysis_%01d.png", sep = '/')
    print(outfile)
    png(outfile,width=15, height=4, res=600, units="in")
    print(plot_grid(p3,p1,p4, nrow=1))
    dev.off()
}
###########################################################################################







































