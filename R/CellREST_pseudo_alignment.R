#' binning
#'
#' @param x Positive numerical vector.
#' @param type If `ACGT`, zeros are assigned to T and the rest of the positive values are uniformly binned and assigned accordingly.
#' @description
#' Maps a positive, continuous vector to discrete states T, G, C, A.
#' @returns Character vector with entries T, G, C, A.
#'
#' @export
binning <- function(x = x,
                    type = "ACGT") {

  if(any(x <= 0)) warning("x contains negative values; those will be ignored.")

  if (type == "ACGT") {
    max_x <- max(x)
    min_x <- min(x)
    x_ACGT <- x
    x_ACGT[which(x == 0)] <- "T"
    x_ACGT[which(x > 0 & x <= (1/3 * (max_x - min_x) + min_x))] <- "G"
    x_ACGT[which(x > (1/3 * (max_x - min_x) + min_x) & x <=
                   (2/3 * (max_x - min_x) + min_x))] <- "C"
    x_ACGT[which(x > (2/3 * (max_x - min_x) + min_x))] <- "A"
  }

  return(x_ACGT)
}

#' seurat_to_fasta
#'
#' @param obj Seurat object.
#' @param varf Character vector of highly variable genes names.
#' @param pre Prefix for the FASTA file name.
#' @param binning_method Passed to function `binning()`.
#' @param assay_in Name of Seurat object assay to use.
#' @param dir Directory for storing FASTA file.
#' @description
#' Creates and writes a pseudo-alignment from normalized counts (stored in `obj@assays[[assay_in]]$data`). Genes are subset to genes specified in `varf`.
#'
#' @export
seurat_to_fasta <- function(obj,
                            varf = NULL,
                            pre = unique(obj$orig.ident)[1],
                            binning_method = c("ACGT")[1],
                            assay_in = DefaultAssay(obj),
                            dir = getwd()) {

  if(is.null(varf)) {
    varf <- VariableFeatures(obj)
  }

  countmat <- t(as.matrix(obj@assays[[assay_in]]$data)[varf,])
  countmat_ACGT <- apply(countmat, 2, FUN = binning, list(binning = binning_method))

  countmatrix_concat <- apply(unlist(countmat_ACGT), 1, paste, collapse = "")
  MSA <- data.frame(names(countmatrix_concat), countmatrix_concat)
  colnames(MSA) <- c("cell id", "sequences")
  seqs = as.list(dplyr::pull(MSA, sequences))
  names = dplyr::pull(MSA, 'cell id')

  seqinr::write.fasta(seqs, names,
                      paste0(dir,"/",pre,"_",assay_in,"_",
                             length(varf), "hvg_",
                             binning_method,"binning.fasta"), open = "w", as.string = FALSE)
}
