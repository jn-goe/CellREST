#' import_trees_lh
#'
#' @param dir Directory in which IQ-TREE output files `.runtrees` or `.treefile` and `.log` are stored.
#' @description
#' Imports trees and likelihoods available in `.runtrees` files or multiple `.treefile` and `.log` files.
#' @returns List containing trees as `multiPhylo` object and a numeric vector of likelihoods.
#'
#' @export
import_trees_lh <- function(dir) {

  last_char <- substr(dir, nchar(dir), nchar(dir))
  if(last_char != "/") {
    dir <- paste0(dir, "/")
  }

  dir_all <- list.files(dir)

  if(sum(grepl(".runtrees",dir_all)) == 1) {
    tree_list <- ape::read.tree(paste0(dir,dir_all[grepl(".runtrees",dir_all)]))
    lh_numbers <- likelihoods_from_logs(paste0(dir,dir_all[grepl(".runtrees",dir_all)]))

  } else {
    lh <- dir_all[grepl(".log",dir_all)]
    lh <- lh[!grepl("iqtree_parallel",lh)]
    dir_all <- dir_all[grepl("treefile",dir_all)]

    tree_list <- list()
    lh_numbers <- c()
    for(t_i in 1:length(dir_all)) {
      tree_list[[t_i]] <- ape::read.tree(paste0(dir,
                                                dir_all[t_i]))
      if(length(lh) > 0) {
        lh_numbers[t_i] <- likelihoods_from_logs(paste0(dir,
                                                        lh[t_i]))
      }
    }
    class(tree_list) <- "multiPhylo"
  }

  return(list(tree_list = tree_list,
              lh_numbers = lh_numbers))
}

#' likelihoods_from_logs
#'
#' @param dir Directory in which IQ-TREE output files `.log` is stored.
#' @description
#' Reads likelihoods in respective `.log` files.
#' @returns Numeric vector of likelihoods.
#'
#' @export
likelihoods_from_logs <- function(dir) {
  text_data <- readLines(dir)
  lh_numbers <- regmatches(text_data, gregexpr("(?<=\\[ lh=)-?[0-9]+", text_data, perl = TRUE))
  lh_numbers <- as.numeric(unlist(lh_numbers))

  if(length(lh_numbers)==0) {
    lh_numbers <- regmatches(text_data, gregexpr("(BEST SCORE FOUND : )-?[0-9]+", text_data, perl = TRUE))
    if(!length(unlist(lh_numbers)) == 0) {
      lh_numbers <- lh_numbers[unlist(lapply(lh_numbers, function(x) length(x) > 0))][[1]]
      lh_numbers <- sub("BEST SCORE FOUND : ","",lh_numbers)
      lh_numbers <- as.numeric(unlist(lh_numbers))
    } else {
      warning(paste0(text_data[length(text_data)-3],"\n",
                     "Remove or redo ", sub(".log","",strsplit(dir,"/")[[1]][length(strsplit(dir,"/")[[1]])])))
      lh_numbers <- "remove"
    }
  }
  return(lh_numbers)
}

#' import_cpu_time
#'
#' @param dir Directory, in which of `iqtree_parallel.log` file is stored.
#' @description
#' Reads combined CPU time from `iqtree_parallel.log` file which recorded multiple tree reconstructions.
#' @returns Numeric vector with CPU times for tree reconstructions.
#'
#' @export
import_cpu_time <- function(dir) {

  last_char <- substr(dir, nchar(dir), nchar(dir))
  if(last_char != "/") {
    dir <- paste0(dir, "/")
  }

  dir_all <- list.files(dir)

  cpu <- dir_all[grepl(".log",dir_all)]
  cpu <- cpu[!grepl("iqtree_parallel",cpu)]

  cpu_numbers_vec <- c()
  for(t_i in 1:length(cpu)) {
    text_data <- readLines(paste0(dir,cpu[t_i]))
    cpu_numbers <- regmatches(text_data, gregexpr("(Total CPU time used: )-?[0-9]+", text_data, perl = TRUE))
    cpu_numbers <- cpu_numbers[unlist(lapply(cpu_numbers, function(x) length(x) > 0))]
    cpu_numbers <- sub("Total CPU time used: ","",cpu_numbers)
    cpu_numbers <- as.numeric(unlist(cpu_numbers))
    cpu_numbers_vec <- c(cpu_numbers_vec, cpu_numbers)
  }

  return(cpu_numbers_vec)
}

#' rate_categories_from_iqtree
#'
#' @param dir Directory in which IQ-TREE output files `.iqtree` are stored.
#' @description
#' Reads rate categories from `.iqtree` files.
#' @returns List of gamma shape alpha and rate vectors.
#'
#' @export
rate_categories <- function(dir) {
  text_data <- readLines(dir, n = 100)

  shape <- regmatches(text_data, gregexpr("Gamma shape alpha: (.*)", text_data, perl = TRUE))
  shape <- as.numeric(sub("Gamma shape alpha: ", "", unlist(shape)))

  cat_1 <- regmatches(text_data, gregexpr("  1 (.*)", text_data, perl = TRUE))
  cat_1 <- unlist(cat_1)
  cat_1 <- suppressWarnings(as.numeric((strsplit(x = cat_1, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = cat_1, split = " ")[[1]])))])[2]

  cat_2 <- regmatches(text_data, gregexpr("  2 (.*)", text_data, perl = TRUE))
  cat_2 <- unlist(cat_2)
  cat_2 <- suppressWarnings(as.numeric((strsplit(x = cat_2, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = cat_2, split = " ")[[1]])))])[2]

  cat_3 <- regmatches(text_data, gregexpr("  3 (.*)", text_data, perl = TRUE))
  cat_3 <- unlist(cat_3)
  cat_3 <- suppressWarnings(as.numeric((strsplit(x = cat_3, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = cat_3, split = " ")[[1]])))])[2]

  cat_4 <- regmatches(text_data, gregexpr("  4 (.*)", text_data, perl = TRUE))
  cat_4 <- unlist(cat_4)
  cat_4 <- suppressWarnings(as.numeric((strsplit(x = cat_4, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = cat_4, split = " ")[[1]])))])[2]

  rates <- c(cat_1, cat_2, cat_3, cat_4)
  return(list(shape = shape,
              rates = rates))
}

#' rate_matrix_from_iqtree
#'
#' @param dir Directory in which IQ-TREE output files `.iqtree` are stored.
#' @description
#' Reads and plots rate matrix Q from `.iqtree` file as transition diagram.
#'
#' @export
rate_matrix_Q_from_iqtree <- function(dir) {
  text_data <- readLines(dir, n = 100)

  Q_row_A <- regmatches(text_data, gregexpr("  A (.*)", text_data, perl = TRUE))
  Q_row_A <- unlist(Q_row_A)
  Q_row_A <- suppressWarnings(as.numeric((strsplit(x = Q_row_A, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = Q_row_A, split = " ")[[1]])))])

  Q_row_C <- regmatches(text_data, gregexpr("  C (.*)", text_data, perl = TRUE))
  Q_row_C <- unlist(Q_row_C)
  Q_row_C <- suppressWarnings(as.numeric((strsplit(x = Q_row_C, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = Q_row_C, split = " ")[[1]])))])

  Q_row_G <- regmatches(text_data, gregexpr("  G (.*)", text_data, perl = TRUE))
  Q_row_G <- unlist(Q_row_G)
  Q_row_G <- suppressWarnings(as.numeric((strsplit(x = Q_row_G, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = Q_row_G, split = " ")[[1]])))])

  Q_row_T <- regmatches(text_data, gregexpr("  T (.*)", text_data, perl = TRUE))
  Q_row_T <- unlist(Q_row_T)
  Q_row_T <- suppressWarnings(as.numeric((strsplit(x = Q_row_T, split = " ")[[1]]))[!is.na(as.numeric((strsplit(x = Q_row_T, split = " ")[[1]])))])

  ratematrix <- matrix(NA, ncol = 4, nrow = 4)
  colnames(ratematrix) <- rownames(ratematrix) <- c("A","C","G","T")

  ratematrix["A",] <- Q_row_A
  ratematrix["C",] <- Q_row_C
  ratematrix["G",] <- Q_row_G
  ratematrix["T",] <- Q_row_T

  g <- igraph::graph_from_adjacency_matrix(ratematrix,
                                           mode = "directed",
                                           weighted = T,
                                           diag = F)
  layout <- data.frame("x" = c(1,1,-1,-1),
                       "y" = c(-1,1,-1,1))
  rownames(layout) <- c("C","G","A","T")
  curves <- igraph::autocurve.edges(g)

  plot(g,
       vertex.size = 80,
       vertex.color = c("A" = "#09529cff",
                        "C" = "#4c81b9ff",
                        "G" = "#a3bedeff",
                        "T" = "#edf1feff"),
       vertex.label.color = c("A" = "white",
                              "C" = "white",
                              "G" = "black",
                              "T" = "black"),
       vertex.label.cex = 2,
       vertex.shape="square",
       #vertex.label.family = "Arial",
       label.family = "Sans",
       vertex.frame.color = "transparent",
       edge.width = .1+5*E(g)$weight/max(E(g)$weight),
       edge.label = E(g)$weight,
       edge.label.color = "grey20",
       #edge.label.family = "Arial",
       edge.arrow.size = 1,
       edge.color = "grey80",
       layout = as.matrix(layout)[V(g)$name, ],
       rescale = T,
       edge.curved = .3)
}

#' rate_parameter_R_from_log
#'
#' @param dir Directory in which IQ-TREE output files `.log` are stored.
#' @description
#' Reads and plots rate parameter R from `.log` as transition diagram.
#'
#' @export
rate_parameter_R_from_log <- function(dir) {
  text_data <- readLines(dir)
  rates <- regmatches(text_data, gregexpr("Rate parameters:  (.*)", text_data, perl = TRUE))
  rates <- unlist(rates)

  if(length(rates)==0) {
    warning(paste0("No A-C-G-T rates in ", sub(".log","",strsplit(dir,"/")[[1]][length(strsplit(dir,"/")[[1]])])))
    lh_numbers <- "remove"
  } else if (length(rates) == 2 | length(rates) == 1){
    ratematrix <- matrix(NA, ncol = 4, nrow = 4)
    colnames(ratematrix) <- rownames(ratematrix) <- c("A","C","G","T")

    A_C <- regmatches(text_data, gregexpr("(A-C: )-?[0-9]+.-?[0-9]+", text_data, perl = TRUE))
    A_C <- (unlist(A_C))[length(unlist(A_C))]
    A_C <- as.numeric(sub("A-C: ","",A_C))
    ratematrix["A","C"] <- ratematrix["C","A"] <- A_C

    A_G <- regmatches(text_data, gregexpr("(A-G: )-?[0-9]+.-?[0-9]+", text_data, perl = TRUE))
    A_G <- (unlist(A_G))[length(unlist(A_G))]
    A_G <- as.numeric(sub("A-G: ","",A_G))
    ratematrix["A","G"] <- ratematrix["G","A"] <- A_G

    A_T <- regmatches(text_data, gregexpr("(A-T: )-?[0-9]+.-?[0-9]+", text_data, perl = TRUE))
    A_T <- (unlist(A_T))[length(unlist(A_T))]
    A_T <- as.numeric(sub("A-T: ","",A_T))
    ratematrix["A","T"] <- ratematrix["T","A"] <- A_T

    C_G <- regmatches(text_data, gregexpr("(C-G: )-?[0-9]+.-?[0-9]+", text_data, perl = TRUE))
    C_G <- (unlist(C_G))[length(unlist(C_G))]
    C_G <- as.numeric(sub("C-G: ","",C_G))
    ratematrix["C","G"] <- ratematrix["G","C"] <- C_G

    C_T <- regmatches(text_data, gregexpr("(C-T: )-?[0-9]+.-?[0-9]+", text_data, perl = TRUE))
    C_T <- (unlist(C_T))[length(unlist(C_T))]
    C_T <- as.numeric(sub("C-T: ","",C_T))
    ratematrix["C","T"] <- ratematrix["T","C"] <- C_T

    G_T <- regmatches(text_data, gregexpr("(G-T: )-?[0-9]+.-?[0-9]+", text_data, perl = TRUE))
    G_T <- (unlist(G_T))[length(unlist(G_T))]
    G_T <- as.numeric(sub("G-T: ","",G_T))
    ratematrix["G","T"] <- ratematrix["T","G"] <- G_T

    g <- igraph::graph_from_adjacency_matrix(ratematrix,
                                             mode = "undirected",
                                             weighted = T,
                                             diag = F)
    layout <- data.frame("x" = c(1,1,-1,-1),
                         "y" = c(-1,1,-1,1))
    rownames(layout) <- c("C","G","A","T")

    plot(g,
         vertex.size = 80,
         vertex.color = c("A" = "#09529cff",
                          "C" = "#4c81b9ff",
                          "G" = "#a3bedeff",
                          "T" = "#edf1feff"),
         vertex.label.color = c("A" = "white",
                                "C" = "white",
                                "G" = "black",
                                "T" = "black"),
         vertex.label.cex = 2,
         vertex.shape="square",
         #vertex.label.family = "Arial",
         label.family = "Sans",
         vertex.frame.color = "transparent",
         edge.width = .1+5*E(g)$weight/max(E(g)$weight),
         edge.label = E(g)$weight,
         edge.label.color = "grey20",
         #edge.label.family = "Arial",
         edge.arrow.size = 0,
         edge.color = "grey80",
         layout = as.matrix(layout)[V(g)$name, ],
         rescale = T)
  }
}
