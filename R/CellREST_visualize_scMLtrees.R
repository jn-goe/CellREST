#' tree_list_pngs
#'
#' @param obj Seurat object.
#' @param tree_list List of trees or `multiPhylo` object.
#' @param col_by Name of metadata column to color cells by.
#' @param col_pal Color palette whose names match the unique values of `col_by`.
#' @param plot_OTU Color trees using `ggtree::groupOTU()`.
#' @param verbose boolean
#' @param branch_width Factor of branch widths.
#' @param plot_tip_as_point boolean
#' @param point_size Tip point size.
#' @param check_for_existing_PNGs boolean
#' @param pre Prefix used for name of created plot file.
#' @param dir_plot Directory in which plot file is saved.
#' @param im_width Plot width.
#' @param im_height Plot height.
#' @param im_res Plot resolution.
#' @description
#' Creates and stores tree plots as PNGs.
#'
#' @importFrom ggtree "%<+%"
# #' @param parallel boolean, parallel computation of knn-tree-graphs
# #' @param cores number of cores for parallel computation
# #' @importFrom foreach "%dopar%"
#'
#' @export
tree_list_pngs <- function(obj,
                           tree_list,
                           col_by = "cluster_names",
                           col_pal = NULL,
                           plot_OTU = F,
                           verbose = T,
                           branch_width = 1,
                           plot_tip_as_point = T,
                           point_size = 0.5,
                           check_for_existing_PNGs = T,
                           pre = "",
                           dir_plots = getwd(),
                           im_width = 3,
                           im_height = 2.5,
                           im_res = 900
                           #parallel = F,
                           #cores = parallel::detectCores()-2
                           ) {


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

  col_vec <- c(col_vec,"grey")
  names(col_vec) <- c(as.character(obj@meta.data[,col_by]),NA)

  ##### Create Tree Plot PNG for each Tree #####
  png_list <- rep(NA, length(tree_list))

  if(check_for_existing_PNGs) {
    if(nchar(pre)>0) pre <- paste0(pre,"_")
    png_exist <- which(unlist(lapply(1:length(tree_list), function(x) file.exists(paste0(dir_plots,"/",pre,"treeplot_",x,".png")))))
    for(tree_index in png_exist) {
      png_list[tree_index] <- paste0(dir_plots, "/",pre,"treeplot_",tree_index,".png")
    }
    if(!(sum(is.na(png_list)) == 0)) {
      message("Missing PNGs will be generated.")
    }
  }

  png_treeindex <- function(tree_index) {
    tree <- tree_list[[tree_index]]
    tree$node.label <- paste0("Node",1:tree$Nnode)

    inner_edges <- which(tree$edge[,2] > length(tree$tip.label))
    outer_edges <- which(tree$edge[,2] <= length(tree$tip.label))

    if(plot_OTU) {
      tree_plot <- suppressMessages(groupOTU(.data = ggtree(tree, layout="equal_angle") ,
                                             .node = split(x = colnames(obj),
                                                           f = as.character(obj@meta.data[,col_by])),
                                             'cluster') +
                                      aes(color=cluster, size=branch_width) +
                                      (if(plot_tip_as_point) geom_tippoint(size = point_size)) +
                                      scale_color_manual(values = col_vec) +
                                      scale_size_continuous(range = c(0, branch_width)) +
                                      theme_minimal_plots() +
                                      theme(rect=element_rect(fill="transparent"),
                                            panel.grid = element_blank(),
                                            axis.ticks = element_line(colour = "transparent"),
                                            axis.text = element_blank(),
                                            panel.border = element_rect(colour = "transparent"),
                                            panel.background= element_blank()) +
                                      NoLegend())
      tree_plot <- tree_plot  +
        xlim(min(tree_plot$data$x)-as.numeric(dist(range(tree_plot$data$x))/5),
             max(tree_plot$data$x)+as.numeric(dist(range(tree_plot$data$x))/5)) +
        ylim(min(tree_plot$data$y)-as.numeric(dist(range(tree_plot$data$y))/5),
             max(tree_plot$data$y)+as.numeric(dist(range(tree_plot$data$y))/5))
    } else {

      if(is.numeric(obj@meta.data[tree$tip.label,col_by])) {
        d <- data.frame(node=1:(tree$Nnode+length(tree$tip.label)),
                        color = c(obj@meta.data[tree$tip.label,col_by],
                                  rep(NA, tree$Nnode)),
                        size = rep(branch_width, length(tree$tip.label)+tree$Nnode))
      } else {
        d <- data.frame(node=1:(tree$Nnode+length(tree$tip.label)),
                        color = c(as.character(obj@meta.data[tree$tip.label,col_by]),
                                  rep(NA, tree$Nnode)),
                        size = rep(branch_width, length(tree$tip.label)+tree$Nnode))
      }

      tree_plot <- suppressMessages((ggtree(tree, layout="equal_angle") %<+%
                                       d + aes(color=color, size=size)) +
                                      (if(plot_tip_as_point) geom_tippoint(size = point_size)) +
                                      (if(is.numeric(obj@meta.data[tree$tip.label,col_by])) scale_color_viridis_c()) +
                                      scale_size_continuous(range = c(0, branch_width)) +
                                      (if(!is.numeric(obj@meta.data[tree$tip.label,col_by])) scale_color_manual(values = col_pal)) +
                                      theme_minimal_plots() +
                                      theme(rect=element_rect(fill="transparent"),
                                            panel.grid = element_blank(),
                                            axis.ticks = element_line(colour = "transparent"),
                                            axis.text = element_blank(),
                                            panel.border = element_rect(colour = "transparent"),
                                            panel.background= element_blank()) +
                                      NoLegend())
      tree_plot <- tree_plot +
        xlim(2*min(tree_plot$data$x), 2*max(tree_plot$data$x)) +
        ylim(2*min(tree_plot$data$y), 2*max(tree_plot$data$y))

    }

    if(nchar(pre)>0) paste0(pre,"_")
    try(ggsave(plot = tree_plot,
               filename = paste0(dir_plots, "/",pre,"treeplot_",tree_index,".png"),
               width = im_width, height = im_height, dpi = im_res), silent = T)

    if(verbose) {
      message(paste0("Creating PNG for tree ",tree_index," out of ",length(tree_list)," trees."))
    }
  }

  # if(parallel == T) {
  #   cl <<- parallel::makeCluster(cores)
  #   doParallel::registerDoParallel(cl)
  #
  #   foreach::foreach(tree_index=which(is.na(png_list)), .packages = c("scPhylo")) %dopar% {
  #     png_treeindex(tree_index)
  #   }
  #   on.exit(parallel::stopCluster(cl))
  # } else {
    lapply(which(is.na(png_list)), png_treeindex)
  # }

  png_list <- c()
  for(tree_index in 1:length(tree_list)) {
    if(nchar(pre)>0) paste0(pre,"_")
    png_list[tree_index] <- paste0(dir_plots, "/",pre,"treeplot_",tree_index,".png")
  }

  return(png_list)
}

#' plot_tree_UMAP
#'
#' @param obj Seurat object.
#' @param tree_list List of trees or `multiPhylo` object.
#' @param distance_measures Distance measure that is used to evaluate similarity between trees. Can be a measure compatible with `treespace::treespace()` or `patristic_correlation`.
#' @param save_patristic_dists If `TRUE` returns patristic distances per tree.
#' @param png_list List of paths to tree PNGs.
#' @param lh_numbers Vector of likelihoods that will specify a color code if `png_list = NULL`.
#' @param size_image Size of the PNGs in the tree-embedding visualization.
#' @param min_dist_umap Passed to `uwot::umap()`.
#' @param n_neighbors_umap Passed to `uwot::umap()`.
#' @param spread_umap Passed to `uwot::umap()`.
#' @param seed_umap Passed to `uwot::umap()`.
#' @param show_plots boolean
#' @param save_plots boolean
#' @param pre Prefix used for name of created plot file.
#' @param dir_plot Directory in which plot file is saved.
#' @param im_width Plot width.
#' @param im_height Plot height.
#' @param im_res Plot resolution.
#' @param format File format.
#' @param title Plot title.
# #' @param parallel boolean, parallel computation of knn-tree-graphs
# #' @param cores number of cores for parallel computation
# #' @importFrom foreach "%dopar%"
#' @description
#' Plots and stores tree-embeddings using tree-to-tree similarities based on `distance_measures`.
#'
#' @export
plot_tree_UMAP <- function(obj,
                           tree_list,
                           distance_measures = c("patristic_correlation","RF","wRF","patristic","nNodes")[1],
                           save_patristic_dists = FALSE,
                           png_list = NULL,
                           lh_numbers = NULL,
                           size_image = .1,
                           umap_on_PCs = F,
                           n_PCs = 50,
                           min_dist_umap = .3,
                           n_neighbors_umap = 30,
                           spread_umap = 1,
                           seed_umap = 42,
                           show_plots = F,
                           save_plots = F,
                           pre = "",
                           dir_plot = getwd(),
                           im_width = 4,
                           im_height = 4,
                           im_res = 900,
                           format = c("png", "pdf")[1],
                           title = NULL#,
                           #parallel = F,
                           #cores = parallel::detectCores()-2
                           ) {

  bool_over5000 <- FALSE
  res <- list()
  cell_order <- colnames(obj)

  for(dm in distance_measures) {
    if(dm %in% c("RF","wRF","patristic","nNodes")) {
      distres <- suppressPackageStartupMessages(try(treespace::treespace(tree_list, method = dm, nf = n_PCs),silent=T))
      if(class(distres) == "try-error") {
        message("Some trees have a distance of zero, treespace cannot be called.")
        distance_obj <- matrix(0)
      } else {
        distance_obj <- as.dist(distres$D)
      }
      path_dists_mat <- paste0("No path distances were calculated for ",dm)
    } else if(dm == "patristic_correlation") {

      if(dim(obj)[2] >= 5000) {
        bool_over5000 <- TRUE
        warning("The data contains more than 5000 cells, hence the patristic distances will not be saved.")

        indx <- list()
        indx_df <- data.frame()
        counter <- 0
        for(i in 1:(length(tree_list)-1)) {
          for(j in (i+1):length(tree_list)) {
            counter <- counter + 1
            indx[[counter]] <- c(i,j)
            indx_df <- rbind(indx_df,c(i,j))
          }
        }

        # if(parallel == T) {
        #
        #   cl <<- parallel::makeCluster(cores)
        #   doParallel::registerDoParallel(cl)
        #
        #   path_dists_mat <- foreach::foreach(ind_tree_i=indx, .packages = c("scPhylo"), .inorder = T, .combine = c) %dopar% {
        #     res <- cor(c(castor::get_all_pairwise_distances(tree_list[[ind_tree_i[1]]], only_clades = cell_order)),
        #                c(castor::get_all_pairwise_distances(tree_list[[ind_tree_i[2]]], only_clades = cell_order)))
        #   }
        #
        #   on.exit(parallel::stopCluster(cl))
        # } else {
          path_dists_mat <- c()

          for(ind_tree_i in indx) {
            path_dists_mat <- c(path_dists_mat,
                                cor(c(castor::get_all_pairwise_distances(tree_list[[ind_tree_i[1]]], only_clades = cell_order)),
                                    c(castor::get_all_pairwise_distances(tree_list[[ind_tree_i[2]]], only_clades = cell_order))))
          }

        # }

        distres <- matrix(NA, ncol = length(tree_list), nrow = length(tree_list))
        distres[as.matrix(indx_df)] <- distres[as.matrix(indx_df)[,2:1]] <- path_dists_mat
        diag(distres) <- 1

        distres <- (1 - abs(distres))^0.5 # to ensure triangle inequality
        distance_obj <- as.dist(distres)

      } else {
        # if(parallel == T) {
        #   cl <<- parallel::makeCluster(cores)
        #   doParallel::registerDoParallel(cl)
        #
        #   path_dists_mat <- foreach::foreach(ind_tree_i=1:length(tree_list), .packages = c("scPhylo"), .inorder = T, .combine = cbind) %dopar% {
        #     res <- c(castor::get_all_pairwise_distances(tree_list[[ind_tree_i]], only_clades = cell_order))
        #     res
        #   }
        #
        #   on.exit(parallel::stopCluster(cl))
        # } else {
          path_dists_mat <- matrix(NA, ncol = length(tree_list), nrow = length(cell_order)*length(cell_order))

          for(ind_tree_i in 1:length(tree_list)) {
            path_dists_mat[,ind_tree_i] <- c(castor::get_all_pairwise_distances(tree_list[[ind_tree_i]], only_clades = cell_order))
          }

        # }

        distres <- (1 - abs(cor(path_dists_mat)))^0.5 # to ensure triangle inequality
        distance_obj <- as.dist(distres)
      }
    }
    rm(distres)

    if(sum(unique(c(as.matrix(distance_obj)))) > 0) {
      prcomp_dist <- as.data.frame(cmdscale(distance_obj, k = n_PCs))
      if(dim(prcomp_dist)[2] != n_PCs) {
        warning(paste0("Using only ",dim(prcomp_dist)[2]," Principal Components, since remaining eigenvalues are <= 0."))
        n_PCs <- dim(prcomp_dist)[2]
      }
      colnames(prcomp_dist) <- paste0("PCo_",1:n_PCs)

      if(umap_on_PCs) {
        distance_obj <- dist(prcomp_dist)
      }

      umap_dist <- as.data.frame(uwot::umap(distance_obj,
                                            seed = seed_umap,
                                            min_dist = min_dist_umap,
                                            n_neighbors = min(n_neighbors_umap, length(tree_list)),
                                            spread = spread_umap))
      colnames(umap_dist) <- c("UMAP_1", "UMAP_2")

      if(!is.null(png_list)) {
        plot_df <- cbind(prcomp_dist, umap_dist, "image" = png_list)

        pcoa_plot <- ggplot(plot_df, aes(x = PCo_1, y = PCo_2)) +
          ggimage::geom_image(aes(image=image), size = size_image) +
          xlim(min(plot_df$PCo_1)-as.numeric(dist(range(plot_df$PCo_1))/5),
               max(plot_df$PCo_1)+as.numeric(dist(range(plot_df$PCo_1))/5)) +
          ylim(min(plot_df$PCo_2)-as.numeric(dist(range(plot_df$PCo_2))/5),
               max(plot_df$PCo_2)+as.numeric(dist(range(plot_df$PCo_2))/5)) +
          theme_minimal_plots() + ggtitle(paste0(dm, " - PCoA plot"))

        umap_plot <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
          ggimage::geom_image(aes(image=image), size = size_image) +
          xlim(min(plot_df$UMAP_1)-as.numeric(dist(range(plot_df$UMAP_1))/5),
               max(plot_df$UMAP_1)+as.numeric(dist(range(plot_df$UMAP_1))/5)) +
          ylim(min(plot_df$UMAP_2)-as.numeric(dist(range(plot_df$UMAP_2))/5),
               max(plot_df$UMAP_2)+as.numeric(dist(range(plot_df$UMAP_2))/5)) +
          theme_minimal_plots() + ggtitle(paste0(dm, " - UMAP plot"))
      } else if (!is.null(lh_numbers)){
        plot_df <- cbind(prcomp_dist, umap_dist)
        plot_df$lh_numbers <- as.numeric(lh_numbers)

        pcoa_plot <- ggplot(plot_df, aes(x = PCo_1, y = PCo_2, col = lh_numbers)) +
          geom_point(size = size_image) +
          scale_colour_gradient2(low = "darkred", mid = "white", high = "darkblue", midpoint = mean(lh_numbers), name = "likelihood") +
          xlim(min(plot_df$PCo_1)-as.numeric(dist(range(plot_df$PCo_1))/5),
               max(plot_df$PCo_1)+as.numeric(dist(range(plot_df$PCo_1))/5)) +
          ylim(min(plot_df$PCo_2)-as.numeric(dist(range(plot_df$PCo_2))/5),
               max(plot_df$PCo_2)+as.numeric(dist(range(plot_df$PCo_2))/5)) +
          theme_minimal_plots() + ggtitle(paste0(dm, " - PCoA plot")) + guide_small_legend_cont()

        umap_plot <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, col = lh_numbers)) +
          geom_point(size = size_image) +
          scale_colour_gradient2(low = "darkred", mid = "white", high = "darkblue", midpoint = mean(lh_numbers), name = "likelihood") +
          xlim(min(plot_df$UMAP_1)-as.numeric(dist(range(plot_df$UMAP_1))/5),
               max(plot_df$UMAP_1)+as.numeric(dist(range(plot_df$UMAP_1))/5)) +
          ylim(min(plot_df$UMAP_2)-as.numeric(dist(range(plot_df$UMAP_2))/5),
               max(plot_df$UMAP_2)+as.numeric(dist(range(plot_df$UMAP_2))/5)) +
          theme_minimal_plots() + ggtitle(paste0(dm, " - UMAP plot")) + guide_small_legend_cont()

      } else if (is.null(lh_numbers)){
        plot_df <- cbind(prcomp_dist, umap_dist)
        plot_df$index <- 1:length(tree_list)

        pcoa_plot <- ggplot(plot_df, aes(x = PCo_1, y = PCo_2, col = index)) +
          geom_point(size = size_image) +
          xlim(min(plot_df$PCo_1)-as.numeric(dist(range(plot_df$PCo_1))/5),
               max(plot_df$PCo_1)+as.numeric(dist(range(plot_df$PCo_1))/5)) +
          ylim(min(plot_df$PCo_2)-as.numeric(dist(range(plot_df$PCo_2))/5),
               max(plot_df$PCo_2)+as.numeric(dist(range(plot_df$PCo_2))/5)) +
          theme_minimal_plots() + scale_color_viridis_c() + guide_small_legend_cont() +
          ggtitle(paste0(dm, " - PCoA plot"))

        umap_plot <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, col = index)) +
          geom_point(size = size_image) + scale_color_viridis_c() +
          xlim(min(plot_df$UMAP_1)-as.numeric(dist(range(plot_df$UMAP_1))/5),
               max(plot_df$UMAP_1)+as.numeric(dist(range(plot_df$UMAP_1))/5)) +
          ylim(min(plot_df$UMAP_2)-as.numeric(dist(range(plot_df$UMAP_2))/5),
               max(plot_df$UMAP_2)+as.numeric(dist(range(plot_df$UMAP_2))/5)) +
          theme_minimal_plots() + ggtitle(paste0(dm, " - UMAP plot")) + guide_small_legend_cont()
      }

      if(!is.null(title)) umap_plot <- umap_plot + ggtitle(title)
      if(!is.null(title)) pcoa_plot <- pcoa_plot + ggtitle(title)

      if(save_patristic_dists & !bool_over5000) {
        res[[dm]] <- list(plot_df = plot_df,
                          path_dists_mat = path_dists_mat)
      } else {
        res[[dm]] <- list(plot_df = plot_df)
      }

      if(show_plots) {
        print(pcoa_plot)
        print(umap_plot)
      }

      if(save_plots) {
        if(nchar(pre)>0) paste0(pre,"_")
        try(ggsave(plot = pcoa_plot, filename = paste0(dir_plot,"/",pre,"",dm,"_treepcoa.",format),
                   width = im_width, height = im_height, dpi = im_res), silent = T)
        try(ggsave(plot = umap_plot, filename = paste0(dir_plot,"/",pre,"",dm,"_treeumap.",format),
                   width = im_width, height = im_height, dpi = im_res), silent = T)
      }
    } else {
      res[[dm]] <- "Trees have a distance of zero."
    }
    print(paste0(dm, " finished."))
  }
  names(res) <- distance_measures
  return(res)
}
