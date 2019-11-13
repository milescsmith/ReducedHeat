#' @title ReducedHeat
#' @description Compute Heatmap using dimensional reduction
#'
#' @param object Object or expression matrix to visualize. Rownames should be gene names
#' @param features Features to plot
#' @param dim_embedding Cell embedding within the dimension of interest
#' @param cell_labels Labels for each cell. While this is generally the cluster identity for a cell, it can be any discrete metadata variable
#' @param enrich For every item in \code{features}, find a number of features equal to \code{enrich} that are close to that feature in the distance induced by the 1D UMAP
#' @param breaks Number of bins
#' @param slope For better visualization, transform the values with a logistic function. This is the slope of that function.
#' @param intercept For better visualization, transform the values with a logistic function. This is the intercept of that function.
#' @param assay Assay to use. Default: the object's current default assay.
#' @param slot Slot to use. Default: "data"
#' @param ident Object variable to use when determining cell_labels
#' @param reduction Reduction to visualize. Default: 'umap'
#' @param dimension The dimension to visualize. Default: 1
#' @param dendrogram See heatmaply help. Default: "row"
#' @param titleX See heatmaply help. Default: FALSE.
#' @param RowV See heatmaply help. Default: TRUE
#' @param show_legend See heatmaply help. Default: FALSE
#' @param hide_colorbar See heatmaply help. Default: FALSE
#' @param fontsize_row See heatmaply help. Default: 7
#' @param margins See heatmaply help. Default: c(70, 50, NA, 0)
#' @param ... Additional parameters to pass to heatmaply
#'
#' @importFrom Matrix colSums rowSums t
#' @importFrom Hmisc %nin%
#' @importFrom furrr future_map
#' @importFrom purrr is_empty
#' @importFrom heatmaply heatmaply
#' @importFrom pdist pdist
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom glue glue
#'
#' @return
#' @export
#'
#' @examples
ReducedHeatmap <- function(object, ...){
  UseMethod("ReducedHeatmap")
}

#' @rdname ReducedHeatmap
#' @method ReducedHeatmap default
#' @importFrom methods as
#'
#' @return
#' @export
ReducedHeatmap.default <- function(object,
                           features,
                           dim_embedding,
                           cell_labels,
                           enrich = 0,
                           breaks = 100,
                           slope = 50,
                           intercept = 0.05,
                           dendrogram = "row",
                           titleX = FALSE,
                           RowV = TRUE,
                           show_legend = FALSE,
                           hide_colorbar = FALSE,
                           fontsize_row = 7,
                           margins = c(70, 50, NA, 0),
                           ...) {
  if (class(object) != "dgCMatrix") {
    object <- as(object, "dgCMatrix")
  }

  notinrows <- (features %nin% rownames(object))
  if (sum(notinrows) > 0) {
    print(glue("The following rows are not in the matrix: {features[notinrows]}"))
  }

  features <- features[!notinrows]
  if (enrich == 0) {
    object <- object[features, ]
  }

  dim_bins <- cut(dim_embedding,
                   breaks = breaks)
  bin_counts <- aggregate.Matrix(t(object), dim_bins)
  empty_bins <- levels(dim_bins)[rowSums(bin_counts) == 0]

  if (!is_empty(empty_bins)){
    empty_bins_matrix <- matrix(0,
                                nrow = length(empty_bins),
                                ncol = ncol(bin_counts))
    rownames(empty_bins_matrix) <- empty_bins
    bin_counts <- rbind(bin_counts,
                        empty_bins_matrix)
  }

  bin_counts <- bin_counts[levels(dim_bins)[which(levels(dim_bins) %in% rownames(bin_counts))], ]
  bin_counts_s <- t(t(bin_counts) / rowSums(t(bin_counts)))

  if (enrich > 0) {
    enriched_features_list <- list()
    message("Determining enrichment")
    pdisttest <- pdist(t(bin_counts_s),
                       indices.A = features,
                       indices.B = 1:ncol(bin_counts_s))
    sortedpdist <- t(apply(as.matrix(pdisttest), 1, order, decreasing = FALSE))

    enriched_features <- sortedpdist[, 1:(enrich + 1)] %>%
      t() %>%
      as.vector() %>%
      unique()
    future_map(.x = 1:length(features),
               .progress = TRUE,
               .f = function(featuresi){
                 enriched_features_list[[features[featuresi]]] <-
                   rownames(object)[sortedpdist[featuresi, 2:(enrich + 1)]]
                 })
    bin_counts_s <- bin_counts_s[, enriched_features]
  }

  # assign a label to each column based on which of the cell_labels is the most common
  cluster_reduction <- data.frame(x = dim_embedding,
                                  y = cell_labels)
  cluster_group <- split(cluster_reduction,
                        cut(cluster_reduction$x,
                            breaks = breaks))

  group_labels <- future_map(cluster_group, function(x) {
      ux <- unique(x[, 2])
      ux[which.max(tabulate(match(x[, 2], ux)))]
    }) %>%
    unlist()

  rownames(bin_counts_s) <- glue("bin:{rownames(bin_counts_s)}, label: {unlist(group_labels[!is.na(group_labels)])}")

  group_labels <- group_labels[!is.na(group_labels)]
  group_labels <- data.frame(Group = group_labels)
  if (enrich > 0) {
    gene_colors <- rep(0, length(enriched_features))
  } else {
    gene_colors <- c()
  }
  gene_colors[colnames(bin_counts_s) %in% features] <- 1
  gene_colors <- data.frame(enriched_features = gene_colors)

  my_palette <- colorRampPalette(c("white", "red"))(n = 1000)
  row_color_palette <- colorRampPalette(c("white", "blue"))
  col_color_palette <- colorRampPalette(brewer.pal(n = 11, "Spectral"))
  toplot <- t(bin_counts_s)

  toplot2 <- 1 / (1 + exp(slope * intercept - slope * toplot))

  heatmap <- heatmaply(as.matrix(toplot2),
                       col_side_colors = group_labels,
                       row_side_colors = gene_colors,
                       row_side_palette = row_color_palette,
                       showticklabels = c(FALSE, TRUE),
                       col = my_palette,
                       dendrogram = dendrogram,
                       titleX = titleX,
                       RowV = RowV,
                       show_legend = show_legend,
                       hide_colorbar = hide_colorbar,
                       fontsize_row = fontsize_row,
                       margins = margins,
                       ...)
  return(heatmap)
}

#' @rdname ReducedHeatmap
#' @method ReducedHeatmap Seurat
#' @importFrom Seurat FetchData GetAssayData Embeddings Idents
#' @return
#' @export
ReducedHeatmap.Seurat <- function(object,
                               assay = NULL,
                               slot = "data",
                               features = NULL,
                               ident = NULL,
                               reduction = 'tsne',
                               dimension = 1,
                               ...){
  
  assay <- assay %||% DefaultAssay(object)
  
  if (is.null(ident)){
    ident <- Idents(object)
  } else {
    ident_df <- FetchData(object,
                          vars = ident)
    ident <- ident_df[,1]
    names(ident) <- rownames(ident_df)
  }
  if (is.null(features)){
    stop("You must provide a set of features to map.")
  }

  exprs <- GetAssayData(object = object,
                        assay = assay,
                        slot = slot)
  
  if(reduction %nin% names(object)){
    stop(glue("{reduction} has not been performed on this object."))
  }
  
  dim_embed <- Embeddings(object = object,
                           reduction = reduction)[,dimension]
  ReducedHeatmap.default(
    object = exprs,
    features = features,
    dim_embedding = dim_embed,
    cell_labels = ident,
    ...)
}

#' @rdname ReducedHeatmap
#' @method ReducedHeatmap SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @importFrom SummarizedExperiment colData assay
#' @importFrom dplyr select_if
#' @return
#' @export
ReducedHeatmap.SingleCellExperiment <- function(object,
                                  assay = "logcounts",
                                  features = NULL,
                                  ident = NULL,
                                  reduction = 'tsne',
                                  dimension = 1,
                                  ...){
  
  if (is.null(ident)){
    if ("ident" %in% names(colData(object))){
      ident <- colData(object)[["ident"]]
    } else if ("cluster" %in% names(colData(object))){
      ident <- colData(object)[["cluster"]]
    } else {
      # Since SingleCellExperiment objects do not have a default
      # cluster identity, guess at one by taking the first column
      # containing factor variables
      ident_tbl <- colData(object) %>%
        as.data.frame %>%
        select_if(is.factor)
      # of course, there may be no such columns...
      if(ncol(ident_tbl) > 0){
        ident <- ident_tbl[,1]
      } else {
            stop("No identity variable given and unable to guess a identity variable.
                 Please provide one.")
          }
    }
  } else {
    ident <- colData(object)[[ident]]
    names(ident) <- rownames(colData(object))
  }
  if (is.null(features)){
    stop("You must provide a set of genes to map.")
  }

  exprs <- assay(x = object, i = assay)

  if(toupper(reduction) %nin% reducedDimNames(object)){
    stop(glue("{reduction} has not been performed on this object."))
  }
  
  dim_embed <- reducedDim(
    x = object,
    type = toupper(reduction))[,dimension]
  
  ReducedHeatmap.default(
    object = exprs,
    features = features,
    dim_embedding = dim_embed,
    cell_labels = ident,
    ...)
}
