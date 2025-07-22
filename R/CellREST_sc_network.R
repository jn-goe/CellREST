#' knn_tree
#'
#' @param obj Seurat object.
#' @param tree_index Index of tree in list of trees `tree_list`.
#' @param tree_list List of trees or `multiPhylo` object.
#' @param pre_comp_path_dists_mat Pre-computed distances output by `plot_tree_UMAP()` (optional).
#' @param tree Phylogenetic tree of class `phylo`. If `NULL`, tree in `tree_list[[tree_index]]` is used.
#' @param k Number of neighbors for the nearest neighbor graph construction. If `NULL` the maximal neighbor distance in `k_dist` is considered.
#' @param k_dist Maximal distance for the nearest neighbor graph construction, if `k = NULL`.
#' @param branch_counts If `TRUE` branch counts are used for the evaluation of the maximal distance `k_dist`, if `FALSE`, path lengths are used.
#' @param return_graph If `TRUE` the `igraph` object is returned, if `FALSE` a list of edges is returned.
#' @description
#' Creates a nn-tree-graph from a phylogenetic tree using either `k` nearest neighbors or all neighbors closer than `k_dist`. The phylogenetic tree can be either input directly via `tree`, or a via specifiying a list of trees `tree_list` and
#' `tree_index`. Optionally, pre-computed patristic distances `pre_comp_path_dists_mat` output by `plot_tree_UMAP()` can be input to save computational time.
#'
#' @export
knn_tree <- function(obj,
                     tree_index = 1,
                     tree_list = tree_list,
                     pre_comp_path_dists_mat = NULL,
                     tree = NULL,
                     k = NULL,
                     k_dist = 1,
                     branch_counts = T,
                     return_graph = F) {

  if(is.null(tree)) {
    tree <- tree_list[[tree_index]]

    if(!is.null(pre_comp_path_dists_mat)) {
      distmat_branchlength <- matrix(pre_comp_path_dists_mat[,tree_index], ncol = dim(obj)[2], nrow = dim(obj)[2])
      colnames(distmat_branchlength) <- rownames(distmat_branchlength) <- colnames(obj)
    } else {
      distmat_branchlength <- castor::get_all_pairwise_distances(tree, only_clades = colnames(obj))
      colnames(distmat_branchlength) <- rownames(distmat_branchlength) <- colnames(obj)
    }
  } else {
    distmat_branchlength <- castor::get_all_pairwise_distances(tree, only_clades = colnames(obj))
    colnames(distmat_branchlength) <- rownames(distmat_branchlength) <- colnames(obj)
  }

  if(branch_counts) {
    distmat_branch <- castor::get_all_pairwise_distances(tree, only_clades = colnames(obj), as_edge_counts = TRUE)
    colnames(distmat_branch) <- rownames(distmat_branch) <- colnames(obj)

  } else {
    distmat_branch <- distmat_branchlength
  }

  if(!is.null(k)) {
    knn <- t(apply(distmat_branch, 1, function(x) colnames(distmat_branch)[order(x)[1:k]]))
    edgelist <- data.frame("x" = colnames(distmat_branch),
                           "y" = c(knn))
  } else {
    nn <- t(apply(distmat_branch, 1, function(x) colnames(distmat_branch)[which(x <= k_dist)]))
    names(nn) <- colnames(distmat_branch)
    if(is.list(nn)) {
      edgelist <- data.frame("x" = rep(names(nn), lapply(nn, length)),
                             "y" = unlist(nn))
    } else {
      edgelist <- data.frame("x" = colnames(distmat_branch),
                             "y" = nn[1,])
    }
  }
  rownames(edgelist) <- NULL

  # remove one vertex loops
  edgelist <- edgelist[which(edgelist[,1] != edgelist[,2]),]

  if(dim(edgelist)[1]>0) {
    # remove duplicate edges
    edgelist[,1:2] <- t(apply(edgelist[,1:2], 1,function(x) sort(x[])))
    edgelist <- edgelist[!duplicated(edgelist),]

    edgelist <- data.frame("node_1" = edgelist[,1],
                           "node_2" = edgelist[,2],
                           "length" = distmat_branchlength[as.matrix(edgelist[,1:2])],
                           "branch_count" = distmat_branch[as.matrix(edgelist[,1:2])])
  }

  if(return_graph) {
    graph <- igraph::graph_from_edgelist(as.matrix(edgelist[,1:2]), directed = F)
    if(dim(edgelist)[1]>0) {
      igraph::E(graph)$length <- edgelist[,3]
      igraph::E(graph)$branch_count <- edgelist[,4]
      igraph::E(graph)$weight <- rep(1, length(igraph::E(graph)$length))
    }
    return(graph)
  } else {
    return(edgelist)
  }
}

#' make_knn_networks
#'
#' @param obj Seurat object.
#' @param tree_list List of trees or `multiPhylo` object.
#' @param k Number of neighbors for the nearest neighbor graph construction. If `NULL` the maximal neighbor distance in `k_dist` is considered.
#' @param k_dist Maximal distance for the nearest neighbor graph construction, if `k = NULL`.
#' @param use_branch_counts If `TRUE` branch counts are used for the evaluation of the maximal distance `k_dist`, if `FALSE`, path lengths are used.
#' @param pre_comp_path_dists_mat Pre-computed distances output by `plot_tree_UMAP()` (optional).
#' @param edge_length_representations Summary statistic to compute representative edge length across nn-tree-graphs.
# #' @param parallel boolean, parallel computation of knn-tree-graphs
# #' @param cores number of cores for parallel computation
# #' @importFrom foreach "%dopar%"
#' @description
#' Constructs network given a list of trees by transforming trees into nnn-tree-graphs and uniting them.
#' @returns Seurat object in which the network is stored as `igraph` object in `obj@misc$cellREST_igraph_obj`.
#'
#' @export
make_knn_networks <- function(obj,
                              tree_list,
                              k = NULL,
                              k_dist = 4,
                              use_branch_counts = T,
                              pre_comp_path_dists_mat = NULL,
                              edge_length_representations = c("min","mean","median")[1]#,
                              #parallel = F,
                              #cores = parallel::detectCores()-2
                              ) {

  # if(parallel == T) {
  #   cl <<- parallel::makeCluster(cores)
  #   doParallel::registerDoParallel(cl)
  #
  #   edgelist_long <- foreach::foreach(tree_index=1:length(tree_list), .packages = c("scPhylo"), .combine = bind_rows) %dopar% {
  #     res <- knn_tree(tree_index,
  #                     obj = obj,
  #                     k = k,
  #                     k_dist = k_dist,
  #                     tree_list=tree_list,
  #                     pre_comp_path_dists_mat = pre_comp_path_dists_mat,
  #                     branch_counts = use_branch_counts)
  #
  #   }
  #
  #   on.exit(parallel::stopCluster(cl))
  # } else {
    edgelist_long <- dplyr::bind_rows(lapply(1:length(tree_list),
                                      knn_tree,
                                      obj = obj,
                                      k = k,
                                      k_dist = k_dist,
                                      tree_list=tree_list,
                                      pre_comp_path_dists_mat = pre_comp_path_dists_mat,
                                      branch_counts = use_branch_counts))
  # }

  if(dim(edgelist_long)[1]>0) {
    weight_df <- (edgelist_long %>% dplyr::group_by(node_1,node_2) %>% dplyr::summarise(weight = dplyr::n(), .groups = 'drop'))
    network <- igraph::graph_from_edgelist(as.matrix(weight_df[,c("node_1","node_2")]), directed = F)
    E(network)$weight <- weight_df$weight

    for(e_rep in edge_length_representations) {
      e_rep_df <- (edgelist_long %>% dplyr::group_by(node_1,node_2) %>% dplyr::summarise(repr_length = eval(call(e_rep, length)),
                                                                                         repr_branch_count = eval(call(e_rep, branch_count)),
                                                                                       .groups = 'drop'))
      network <- set_edge_attr(network, paste0(e_rep,"_length"),
                               value = e_rep_df$repr_length)
      network <- set_edge_attr(network, paste0(e_rep,"_branch_count"),
                               value = e_rep_df$repr_branch_count)
    }

  } else {
    network <- igraph::graph_from_edgelist(as.matrix(edgelist[,1:2]), directed = F)
    E(network)$weight <- 0
  }

  unconnected_cells <- colnames(obj)[which(!(colnames(obj) %in% V(network)$name))]
  network <- igraph::add_vertices(network, length(unconnected_cells), attr = list("name" = unconnected_cells))

  # node degree
  node_degree <- igraph::degree(network)
  obj$node_degree <- 0
  obj$node_degree[names(node_degree)] <- node_degree

  obj@misc$cellREST_igraph_obj <- network
  adj_mat <- as_adjacency_matrix(network, attr = "weight")
  adj_mat <- adj_mat[colnames(obj),colnames(obj)]
  obj@graphs$cellREST <- Seurat::as.Graph(adj_mat)

  return(obj)
}

#' embed_network
#'
#' @param obj Seurat object.
#' @param network `igraph` object.
#' @param umap_key Name of dimensionality reduction to use in the Seurat object.
#' @param edge_length_attr Edge attribute of `network` to use for computation of the shortest path distance between cells. If `constant` constant edge lengths are used.
#' @param min_dist_umap Passed to `uwot::umap()`.
#' @param n_neighbors_umap Passed to `uwot::umap()`.
#' @param spread_umap Passed to `uwot::umap()`.
#' @param seed_umap Passed to `uwot::umap()`.
#' @description
#' Embeds the nodes of a network based on edge lengths using UMAP.
#' @returns Seurat object in which the resulting UMAP is stored in `obj@reductions[[umap_key]]` and the
#' the shortest path distance matrix is stored in `obj@misc$cellREST_distmat`.
#'
#' @export
embed_network <- function(obj,
                          network,
                          umap_key = "network_umap",
                          edge_length_attr = c("min_length","min_branch_count")[1],
                          min_dist_umap = .3,
                          n_neighbors_umap = 30,
                          spread_umap = 1,
                          seed_umap = 42) {

  if(dim(as_edgelist(network))[1] == 0 ) {
    warning("The input network has no edges. All cell-to-cell distance are Inf.\nCell-to-cell distances will be set to 1 for UMAP embedding.")
    distmat_graph <- igraph::distances(network)
    distmat_graph_orig <- distmat_graph
    distmat_graph[which(distmat_graph == Inf, arr.ind = T)] <- 1
  } else {
    if(!edge_length_attr == "constant") {
      if(!(edge_length_attr %in% names(edge_attr(network)))) {
        warning(paste0("\"",edge_length_attr, "\" is not in edge attributes. \"", names(edge_attr(network))[1], "\" is chosen instead."))
        edge_length_attr <- names(edge_attr(network))[1]
      }

      if(edge_length_attr == "weight") {
        warning("Inverse of ",edge_length_attr," edge attribute values will be used distances.")
        weight_vec <- max(edge_attr(network)[[edge_length_attr]])/edge_attr(network)[[edge_length_attr]]
      } else {
        weight_vec <- edge_attr(network)[[edge_length_attr]]
      }
    } else {
      weight_vec <- rep(1, igraph::ecount(network))
    }

    distmat_graph <- igraph::distances(network, weight = weight_vec)
    distmat_graph <- distmat_graph[colnames(obj),colnames(obj)]
    distmat_graph_orig <- distmat_graph

    #-- if unconnected, add high distance ensuring triangle inequality --#
    if(any(c(distmat_graph) == Inf)) {
      warning(paste0("Network has unconnected components, hence some cell-to-cell distances are Inf.\nInfinite cell-to-cell distances will be set to ",
                     round(3*sort(unique(c(distmat_graph)), decreasing = T)[2],2)," for UMAP embedding."))
      distmat_graph[which(distmat_graph == Inf, arr.ind = T)] <- round(3*sort(unique(c(distmat_graph)), decreasing = T)[2],2)
    }

  }
  distmat_graph <- as.dist(distmat_graph)

  umap_new <- uwot::umap(distmat_graph,
                         n_neighbors = n_neighbors_umap,
                         spread = spread_umap,
                         min_dist = min_dist_umap,
                         seed = seed_umap)
  suppressWarnings(obj[[umap_key]] <- CreateDimReducObject(embeddings = umap_new, key = umap_key))

  if(edge_length_attr == "min_length") {
    obj@misc$cellREST_distmat <- as.dist(distmat_graph_orig)
  } else {
    obj@misc$distmat <- as.dist(distmat_graph_orig)
  }
  return(obj)
}

#' plot_network
#'
#' @param obj Seurat object.
#' @param network `igraph` object.
#' @param umap_key Name of dimensionality reduction that is used for either plotting the network or for initializing a graph layout.
#' @param layout_by_weight If `TRUE`, a force-directed graph layout is computed using network edge weights.
#' @param layout_model Either `layout_with_kk`, or `layout_with_fr` and is passed to `igraph` if `layout_by_weight = TRUE`.
#' @param node_size Size of nodes.
#' @param edge_width_factor Scaling factor of edge widths.
#' @param col_by Name of metadata column to color cells by.
#' @param col_pal Color palette whose names match the unique values of `col_by`.
#' @param edge.color Vector of colors for edges.
#' @param seed_graph_layout Seed for computation of force-directed graph layout if `layout_by_weight = TRUE`.
#' @param show_plots boolean
#' @param save_plots boolean
#' @param pre Prefix used for name of created PNG.
#' @param dir_plot Directory in which PNG file is saved.
#' @param im_width Plot width.
#' @param im_height Plot height.
#' @param im_res Plot resolution.
#' @description
#' Plots `network` and saves plot as PNG.
#'
#' @export
plot_network <- function(obj,
                         network,
                         umap_key = "network_umap",
                         igraph_layout = FALSE,
                         igraph_layout_type = c("layout_with_kk",
                                                "layout_with_fr")[1],
                         layout_weight_attr = c("weight",
                                                "min_branch_count",
                                                "min_length")[1],
                         edge_width_attr = c("weight",
                                             "min_branch_count",
                                             "min_length")[1],
                         node_size = 6,
                         edge_width_factor = .2,
                         seed_graph_layout = 1,
                         col_by = "cluster_names",
                         col_pal = NULL,
                         edge.color = NULL,
                         show_plots = T,
                         save_plots = F,
                         pre = "",
                         dir_plot = getwd(),
                         im_width = 1,
                         im_height = 1,
                         im_res = 900) {

  colorvertex <- rep("grey",length(igraph::V(network)))
  names(colorvertex) <- names(igraph::V(network))

  if(!is.null(col_by)) {
    if(is.numeric(obj@meta.data[,col_by])) {
      col_vec <- numerical_to_viridis(obj@meta.data[,col_by])
      names(col_vec) <- colnames(obj)
    } else {
      if(!is.null(col_pal)) {
        col_vec <- col_pal[as.character(obj@meta.data[,col_by])]
        names(col_vec) <- colnames(obj)
      } else {

        col_pal <- scales::hue_pal()(length(unique(as.character(obj@meta.data[,col_by]))))
        names(col_pal) <- unique(as.character(obj@meta.data[,col_by]))

        col_vec <- col_pal[as.character(obj@meta.data[,col_by])]
        names(col_vec) <- colnames(obj)
      }
    }
    colorvertex[names(col_vec)] <- col_vec
  } else {
    colorvertex <- scales::alpha(colour = "grey",alpha = .3)
  }

  umap_df <- obj@reductions[[umap_key]]@cell.embeddings
  umap_df <- umap_df[V(network)$name,]

  if(igraph_layout) {
    set.seed(seed_graph_layout)

    if(igraph_layout_type == "layout_with_kk") {
      if(layout_weight_attr == "weight") {
        warning("Since `layout_with_kk` is chosen, inverse ",layout_weight_attr," edge attribute values will be used as weights.")
        weight_vec <- max(edge_attr(network)[[layout_weight_attr]])/
          edge_attr(network)[[layout_weight_attr]]
      } else {
        weight_vec <- edge_attr(network)[[layout_weight_attr]]/
          max(edge_attr(network)[[layout_weight_attr]])
      }
      graph_layout <- layout_with_kk(network,
                                     coords = umap_df,
                                     weights = weight_vec)
    } else if(igraph_layout_type == "layout_with_fr") {
      if(layout_weight_attr %in% c("min_length","min_branch_count")) {
        warning("Since `layout_with_kk` is chosen, inverse ",layout_weight_attr," edge attribute values will be used as weights.")
        weight_vec <- max(edge_attr(network)[[layout_weight_attr]])/
          edge_attr(network)[[layout_weight_attr]]
      } else {
        weight_vec <- edge_attr(network)[[layout_weight_attr]]/
          max(edge_attr(network)[[layout_weight_attr]])
      }
      graph_layout <- layout_with_fr(network,
                                     coords = umap_df,
                                     weights = weight_vec)
    }

    layout_model <- paste0(igraph_layout_type, "_on_", umap_key)
  } else {
    layout_model <- umap_key
    graph_layout <- umap_df
  }

  if(is.null(edge_width_attr)) {
    width_vec <- edge_width_factor
  } else if(edge_width_attr %in% c("min_length","min_branch_count")) {
    warning("Inverse of ",edge_width_attr," edge attribute values will be used for edge width.")
    width_vec <- edge_width_factor/edge_attr(network)[[edge_width_attr]]
  } else {
    width_vec <- edge_width_factor*edge_attr(network)[[edge_width_attr]]/max(edge_attr(network)[[edge_width_attr]])
  }

  if(is.null(edge.color)) edge.color <- "grey"

  if(save_plots) {
    if(nchar(pre)>0) paste0(pre,"_")
    png(paste0(dir_plot,"/",pre,"network_on_",layout_model,".png"),
        width = im_width, height = im_height,
        units = "in", bg = "transparent", res = im_res)
    par(mar=c(0,0,0,0))
    plot(network,
         vertex.size = node_size,
         vertex.label = NA,
         vertex.color = colorvertex,
         vertex.frame.width = 0,
         vertex.frame.color = "transparent",
         edge.color = edge.color,
         edge.width = width_vec,
         layout = graph_layout,
         rescale = T)
    dev.off()
  }

  if(show_plots) {
    plot(network,
         vertex.size = node_size,
         vertex.label = NA,
         vertex.color = colorvertex,
         vertex.frame.width = 0,
         vertex.frame.color = "transparent",
         edge.color = edge.color,
         edge.width = width_vec,
         layout = graph_layout,
         rescale = T)
  }
}

#' cluster_heatmap
#'
#' @param obj Seurat object.
#' @param distmat Distance matrix.
#' @param group_by Name of metadata column to group cells by.
#' @param col_by Name of metadata column to color cells by.
#' @param col_pal Color palette whose names match the unique values of `col_by`.
#' @param sum_stats Summary statistics for heatmap entries if `group_by` groups multiple cells into one category.
#' @param hierarchical_clustering If `TRUE` columns and rows of the heatmap are hierarchically clustered.
#' @param fixed_order Vector of characters indicating a fixed order of heatmap columns and rows. If `NULL` the heatmap columns and rows are not ordered.
#' @param cl_method Clustering method for the hierarchical clustering of heatmap columns and rows (see `pheatmap::pheatmap()` for details).
#' @param show_plots boolean
#' @param save_plots boolean
#' @param pre Prefix used for name of created PNG.
#' @param dir_plot Directory in which PNG file is saved.
#' @param im_width Plot width.
#' @param im_height Plot height.
#' @param ... Additional parameters passed to `pheatmap::pheatmap()`.
#' @description
#' Wrapper of `pheatmap::pheatmap()` visualizing `distmat` values as a heatmap and saving it as PNG.
#'
#' @export
cluster_heatmap <- function(obj,
                            distmat,
                            group_by = "cluster_names",
                            col_by = NULL,
                            col_pal = NULL,
                            sum_stats = c("min","mean")[1],
                            hierarchical_clustering = FALSE,
                            fixed_order = NULL,
                            cl_method = c("average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid")[1],
                            show_plots = T,
                            save_plots = F,
                            return_matrix = F,
                            pre = "",
                            dir_plot = getwd(),
                            save_height = 4,
                            save_width = 4.6,
                            ... ){

  if(any(is.infinite(c(distmat)))) {
    warning("Network has unconnected components, hence some cell-to-cell distances are Inf.
                    Infinite cell-to-cell distances will be set to NA")
    distmat[which(is.infinite(distmat), arr.ind = T)] <- NA
  }

  if(!(any(obj@meta.data[, group_by] != colnames(obj)))) {
    if(is.null(fixed_order)) {
      clusters <- as.character(obj@meta.data[, col_by])
    } else {
      clusters <- as.character(obj@meta.data[fixed_order, col_by])
    }

    distmat_percluster <- distmat
    breaks <- NULL
    sum_stats <- ""

  } else {
    if(is.null(col_by)) {
      col_by <- group_by
    }
    if(col_by != group_by) {
      warning("Heatmap columns and columns will be coloured by ", group_by)
      col_by <- group_by
    }
    if(is.factor(obj@meta.data[, group_by])) {
      clusters <- as.character(levels(obj@meta.data[, group_by]))
    } else {
      clusters <- as.character(unique(obj@meta.data[, group_by]))
    }

    distmat_percluster <- matrix(0, nrow = length(clusters), ncol = length(clusters))
    colnames(distmat_percluster) <- rownames(distmat_percluster) <- clusters

    for(cl_i in 1:(length(clusters)-1)) {
      for(cl_j in cl_i:length(clusters)) {
        cl_i_cells <- colnames(obj)[which(as.character(obj@meta.data[, group_by]) == clusters[cl_i])]
        cl_j_cells <- colnames(obj)[which(as.character(obj@meta.data[, group_by]) == clusters[cl_j])]

        distmat_percluster[cl_i, cl_j] <- distmat_percluster[cl_j,cl_i] <- eval(call(sum_stats, as.matrix(distmat)[cl_i_cells,cl_j_cells]))
      }
    }

    if(length(unique(c(distmat_percluster))) == 1) {
      breaks <- c(unique(c(distmat_percluster)), unique(c(distmat_percluster)) + 1)
    } else {
      breaks <- NULL
    }
  }

  if(!is.null(fixed_order)) {
    distmat_percluster <- as.dist(as.matrix(distmat_percluster)[fixed_order, fixed_order])
  }

  distmat_percluster[which(is.na(distmat_percluster), arr.ind = T)] <- NA

  if(nchar(sum_stats) > 0) sum_stats <- paste0(sum_stats,"_")

  if(is.null(col_by)) {
    if(save_plots) {
      if(nchar(pre)>0) paste0(pre,"_")
      png(paste0(dir_plot,"/",pre,"",sum_stats,"heatmap.png"),
          width = save_width, height = save_height,
          units = "in", bg = "transparent", res = 150)
      par(mar=c(0,0,0,0)+.1)
      pheatmap::pheatmap(distmat_percluster,
                         scale = "none",
                         color = pals::coolwarm(50),
                         annotation_names_row = F,
                         annotation_names_col = F,
                         cluster_cols = hierarchical_clustering,
                         cluster_rows = hierarchical_clustering,
                         breaks = breaks,
                         na_col = "grey",
                         clustering_method = cl_method,
                         show_rownames = F, show_colnames = F,
                         fontsize = 20,
                         ...)
      dev.off()
    }
    if(show_plots) {
      print(pheatmap::pheatmap(distmat_percluster,
                               scale = "none",
                               color = pals::coolwarm(50),
                               annotation_names_row = F,
                               annotation_names_col = F,
                               cluster_cols = hierarchical_clustering,
                               cluster_rows = hierarchical_clustering,
                               breaks = breaks,
                               na_col = "grey",
                               clustering_method = cl_method,
                               show_rownames = F, show_colnames = F,
                               ...))
    }
  } else {
    anno_df <- data.frame("clusters" = clusters)
    rownames(anno_df) <- colnames(as.matrix(distmat_percluster))
    colnames(anno_df) <- col_by

    if(is.null(col_pal)) {
      if(!is.numeric(obj@meta.data[, col_by])) {
        col_pal <- scales::hue_pal()(length(clusters))
        names(col_pal) <- clusters
      } else {
        col_pal <- numerical_to_viridis(obj@meta.data[, col_by])
        names(col_pal) <- as.character(obj@meta.data[, col_by])
      }
    }

    anno_col <- list("clusters" = col_pal)
    names(anno_col) <- col_by
    if(save_plots) {
      if(nchar(pre)>0) paste0(pre,"_")
      png(paste0(dir_plot,"/",pre,"",sum_stats,"heatmap.png"),
          width = save_width, height = save_height,
          units = "in", bg = "transparent", res = 150)
      par(mar=c(0,0,0,0)+.1)
      pheatmap::pheatmap(distmat_percluster,
                         scale = "none",
                         color = pals::coolwarm(50),
                         annotation_names_row = F,
                         annotation_names_col = F,
                         annotation_row = anno_df,
                         annotation_col = anno_df,
                         annotation_colors = anno_col,
                         annotation_legend = F,
                         cluster_cols = hierarchical_clustering,
                         cluster_rows = hierarchical_clustering,
                         breaks = breaks,
                         na_col = "grey",
                         clustering_method = cl_method,
                         show_rownames = F, show_colnames = F,
                         fontsize = 20,
                         ...)
      dev.off()
    }
    if(show_plots) {
      print(pheatmap::pheatmap(distmat_percluster,
                               scale = "none",
                               color = pals::coolwarm(50),
                               annotation_names_row = F,
                               annotation_names_col = F,
                               annotation_row = anno_df,
                               annotation_col = anno_df,
                               annotation_colors = anno_col,
                               annotation_legend = F,
                               cluster_cols = hierarchical_clustering,
                               cluster_rows = hierarchical_clustering,
                               breaks = breaks,
                               na_col = "grey",
                               clustering_method = cl_method,
                               show_rownames = F, show_colnames = F,
                               ...))
    }
  }
  if(return_matrix) return(distmat_percluster)
}


#' edge_weight_length_plot
#'
#' @param obj Seurat object.
#' @param network `igraph` object.
#' @param group_by Name of metadata column to group cells by.
#' @param col_by Name of metadata column to color cells by.
#' @param col_pal Color palette whose names match the unique values of `col_by`.
#' @param x_attr_name Edge attribute of `network` plotted on x-axis.
#' @param y_attr_name Edge attribute of `network` plotted on y-axis.
#' @param size_attr_name Edge attribute of `network` determining size of points.
#' @param pt_maxsize Maximal point size.
#' @param pt_stroke Passed to `ggplot2::geom_point()`.
#' @param pt_alpha Passed to `ggplot2::geom_point()`.
#' @param show_plots boolean
#' @param save_plots boolean
#' @param pre Prefix used for name of created plot file.
#' @param dir_plot Directory in which plot file is saved.
#' @param im_width Plot width.
#' @param im_height Plot height.
#' @param format File format.
#' @param ... Additional parameters passed to `ggplot2::geom_point()`.
#' @description
#' Plots edge lengths vs edge weights (or other attribute combinations) of a network.
#'
#' @export
edge_weight_length_plot <- function(obj,
                                    network,
                                    group_by = "cluster_names",
                                    col_by = "cluster_names",
                                    col_pal = NULL,
                                    x_attr_name = "min_length",
                                    y_attr_name = "weight",
                                    size_attr_name = "min_branch_count",
                                    pt_maxsize = 2,
                                    pt_stroke = 0,
                                    pt_alpha = 0.3,
                                    show_plots = T,
                                    save_plots = F,
                                    pre = "",
                                    dir_plot = getwd(),
                                    format = c("png","pdf")[1],
                                    im_width = 3,
                                    im_height = 2.5) {

  edgelist <- as_edgelist(network)
  edgelist_cl <- data.frame("out_cl" = as.character(obj@meta.data[edgelist[,1],group_by]),
                            "in_cl" = as.character(obj@meta.data[edgelist[,2],group_by]),
                            "x" = edge_attr(network)[[x_attr_name]],
                            "y" = edge_attr(network)[[y_attr_name]],
                            "size" = edge_attr(network)[[size_attr_name]])
  edgelist_cl <- edgelist_cl[,c("out_cl","in_cl","x","y","size")]

  if(is.null(col_pal)) {
    col_pal <- scales::hue_pal()(length(unique(as.character(obj@meta.data[,group_by]))))
    names(col_pal) <- clusters
  }

  for(cl in unique(as.character(obj@meta.data[,group_by]))) {

    data <- edgelist_cl[edgelist_cl$out_cl == cl |
                          edgelist_cl$in_cl == cl,]
    data <- data[!(data$out_cl == cl &
                     data$in_cl == cl), ]

    data[which(data$in_cl == cl),1:2] <- data[which(data$in_cl == cl),2:1]

    size_vec <- data$size
    names(size_vec) <- data$size

    p <- ggplot(data,
                aes(x=x, y=y, col=in_cl, size=size)) +
      geom_point(stroke = pt_stroke, alpha = pt_alpha) +
      scale_size(range=c(pt_maxsize,1),
                 breaks = 1:max(data$size), name = "branch count") +
      scale_color_manual(values = col_pal, name = "", ) +
      ggtitle(paste0("edges from/to ",cl)) +
      xlab("edge length") + ylab("edge weight") +
      theme_minimal() +
      guides(col = guide_legend(override.aes = list(alpha = 1)))

    if(show_plots) print(p)

    if(save_plots) {
      p <- p +
        theme_minimal_plots() +
        guides(size=guide_legend(keywidth=0.1,
                                 keyheight=0.15,
                                 default.unit="inch",
                                 order = 2),
               color=guide_legend(keywidth=0.1,
                                  keyheight=0.15,
                                  default.unit="inch",
                                  override.aes = list(size = 1.5,
                                                      alpha = 1),
                                  order = 1))

      if(nchar(pre)>0) paste0(pre,"_")
      try(ggsave(plot = p, filename = paste0(dir_plot,"/", pre,cl,"_weight_length_plot.",format),
                 width = im_width, height = im_height, dpi = 900), silent = T)
    }
  }
}

#' prune_network
#'
#' @param network `igraph` object.
#' @param quantile Theoretical quantile of the fitted normal distribution on logarithmized edge lengths.
#' @param method If `IQR` standard deviation is estimated based on the Interquantile Range, i.e. IQR/1.349. If `left` standard deviation is estimated using only values smaller than the median.
#' @param show_plots boolean
#' @param save_plots boolean
#' @param pre Prefix used for name of created plot files.
#' @param dir_plot Directory in which plot files is saved.
#' @param format File format.
#' @param pt.size Point size of edge length vs edge weight plot.
#' @param pt.alpha Point alpha of edge length vs edge weight plot.
#' @param im_width_hist Width of edge length histogram image.
#' @param im_height_hist Height of edge length histogram image.
#' @param im_width_dot Width of edge length vs edge weight plot image.
#' @param im_height_dot Height of edge length vs edge weight plot image.
#' @param im_res Plot resolution.
#' @description
#' Prunes edges of a network if their logarithmized lengths exceed the theoretical quantile of a fitted normal distribution
#' and visualizes edge length distribution as histogram and edge length vs edge weight plots.
#' @returns Pruned network as `igraph` object.
#'
#' @export
prune_network <- function(network,
                          quantile = 0.9999,
                          method = c("IQR", "left")[1],
                          show_plots = T,
                          save_plots = F,
                          pre = "",
                          dir_plot = getwd(),
                          format = c("png","pdf")[1],
                          pt.size = 1,
                          pt.alpha = .5,
                          im_width_hist = 1.5,
                          im_height_hist = 1.5,
                          im_width_dot = 1.5,
                          im_height_dot = 1.5,
                          im_res = 900) {

  # log-transform edge lengths of sc-network
  df_plot <- data.frame("log_min_length" = log(E(network)$min_length),
                        "weight" = E(network)$weight)
  ylog <- df_plot$log_min_length[!df_plot$log_min_length == -Inf]

  # median as mean estimator
  norm_median <- median(ylog)

  if(method == "IQR") {
    # estimate standard deviation using IQR
    norm_sd <- IQR(ylog)/1.349
  } else if(method == "left") {
    # estimate standard deviation on values smaller median
    left_vec <- ylog[ylog <= norm_median]
    norm_sd <- (sum((left_vec-norm_median)^2)/length(left_vec))^0.5
  }

  # use theoretical quantile as threshold
  q_lnorm <- qnorm(quantile, norm_median, norm_sd)
  edge.col_red <- rep("grey", dim(as_edgelist(network))[1])
  edge.col_red[which(log(E(network)$min_length) > q_lnorm)] <- "darkred"
  df_plot$edge.col_red <- edge.col_red

  percentage <- (1-as.numeric(table(edge.col_red)["grey"]/length(edge.col_red)))*100
  message(paste0(round(percentage,2) ,"% of edges were pruned based on the theoretical ",quantile*100,"% quantile of the fitted normal distribution."))

  network <- delete_edges(network, which(log(E(network)$min_length) > q_lnorm ))
  obj@misc$cellREST_igraph_obj <- network

  xvec <- seq(min(ylog),max(ylog),0.001)
  hist <- ggplot(data = df_plot,
                 aes(x = log_min_length,
                     y = after_stat(density))) +
    geom_histogram(col = "white", fill = "#93b2b8ff", lwd = .1) +
    geom_line(data = data.frame("x" = xvec,
                                "y" = dnorm(xvec, norm_median, norm_sd)),
              aes(x = x, y = y), col = "#3e4461ff", lwd = .3) +
    geom_vline(xintercept = q_lnorm, col = "darkred", lwd = .3, lty = 2) +
    xlab("log(edge length)") +
    ylab("ratio of cells") +
    xlim(min(ylog),max(ylog,q_lnorm)) +
    theme_minimal_plots(legend.position = "none")

  p <- ggplot(df_plot, aes(x = log_min_length, y = weight, col = edge.col_red)) +
    xlab("log(edge length)") +
    ylab("edge weight") +
    geom_point(size = pt.size, stroke = 0, alpha = pt.alpha) +
    scale_color_manual(name = "", values = c("darkred" = "darkred",
                                             "grey" = "#93b2b8ff")) +
    geom_vline(xintercept = q_lnorm, col = "darkred", lwd = .3, lty = 2) +
    xlim(min(ylog),max(ylog,q_lnorm)) +
    theme_minimal_plots() + NoLegend()

  if(save_plots) {
    if(nchar(pre)>0) paste0(pre,"_")
    try(ggsave(plot = hist,
               filename = paste0(dir_plot,"/", pre, "histogram_cellREST_dists.",format),
               width = im_width_hist, height = im_height_hist, dpi = im_res), silent = T)
    try(ggsave(plot = p, filename = paste0(dir_plot,"/", pre, "length_vs_weight_cellREST_dists.",format),
               width = im_width_dot, height = im_height_dot, dpi = im_res), silent = T)
  }

  if(show_plots) print(cowplot::plot_grid(hist, p, ncol = 1))

  return(network)
}
