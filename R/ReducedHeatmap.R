#' @title ReducedHeat
#' @description Compute Heatmap using dimensional reduction
#'
#' @param expression_matrix Full expression matrix genes (rows) vs. cells (columns). Rownames should be gene names
#' @param genes_of_interest genes_of_interest Genes of Interest
#' @param dim_embedding Cell embedding within the dimension of interest
#' @param cell_labels Labels for each cell. While this is generally the cluster identity for a cell, it can be any discrete metadata variable
#' @param enrich For every gene in genes of interest, find enrich number of genes that are close to that gene in the distance induced by the 1D t-SNE
#' @param breaks Number of bins
#' @param slope For better visualization, transform the values with a logistic function. This is the slope of that function.
#' @param intercept For better visualization, transform the values with a logistic function. This is the intercept of that function.
#'
#' @param object Object to visualize
#' @param assay Assay to use. Default: "RNA"
#' @param slot Slot to use. Default: "data"
#' @param ident Object variable to use when determining cell_labels
#' @param reduction Reduction to visualize. Default: 'tsne'
#' @param dimension The dimension to visualize. Default: 1
#' @param expression_matrix A gene by cell matrix of expression values.
#' If passing an object, this is retrieved from that.
#' @param dendrogram See heatmaply help. Default: "row"
#' @param titleX See heatmaply help. Default: FALSE.
#' @param RowV See heatmaply help. Default: TRUE
#' @param show_legend See heatmaply help. Default: FALSE
#' @param hide_colorbar See heatmaply help. Default: FALSE
#' @param fontsize_row See heatmaply help. Default: 7
#' @param margins See heatmaply help. Default: c(70, 50, NA, 0)
#' @param ... Additional parameters to pass to heatmaply
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
#' @importFrom Matrix colSums rowSums t
#' @importFrom Hmisc %nin%
#' @importFrom furrr future_map
#' @importFrom purrr is_empty
#' @importFrom heatmaply heatmaply
#' @importFrom pdist pdist
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom Matrix.utils aggregate.Matrix
#'
#' @return
#' @export
ReducedHeatmap.default <- function(expression_matrix,
                           genes_of_interest,
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
  if (class(expression_matrix) != "dgCMatrix") {
    expression_matrix <- as(expression_matrix, "dgCMatrix")
  }

  notinrows <- (genes_of_interest %nin% rownames(expression_matrix))
  if (sum(notinrows) > 0) {
    print(glue("The following rows are not in the matrix: {genes_of_interest[notinrows]}"))
  }

  genes_of_interest <- genes_of_interest[!notinrows]
  if (enrich == 0) {
    expression_matrix <- expression_matrix[genes_of_interest, ]
  }


  dim_bins <- cut(dim_embedding,
                   breaks = breaks)
  bin_counts <- aggregate.Matrix(t(expression_matrix), dim_bins)
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
    enriched_genes_list <- list()
    message("Determining enrichment")
    pdisttest <- pdist(t(bin_counts_s),
                       indices.A = genes_of_interest,
                       indices.B = 1:ncol(bin_counts_s))
    sortedpdist <- t(apply(as.matrix(pdisttest), 1, order, decreasing = FALSE))

    enriched_genes <- sortedpdist[, 1:(enrich + 1)] %>%
      t() %>%
      as.vector() %>%
      unique()
    future_map(.x = 1:length(genes_of_interest),
               .progress = TRUE,
               .f = function(genes_of_interesti){
      enriched_genes_list[[genes_of_interest[genes_of_interesti]]] <- rownames(expression_matrix)[sortedpdist[genes_of_interesti, 2:(enrich + 1)]]
    })
    bin_counts_s <- bin_counts_s[, enriched_genes]
  }

  # assign a label to each column based on which of the cell_labels is the most common
  cluster_reduction <- data.frame(x = dim_embedding,
                                  y = cell_labels)
  cluster_group <- split(cluster_reduction,
                        cut(cluster_reduction$x,
                            breaks = breaks))

  group_labels <- map(cluster_group, function(x) {
      ux <- unique(x[, 2])
      ux[which.max(tabulate(match(x[, 2], ux)))]
    }) %>%
    unlist()

  rownames(bin_counts_s) <- glue("bin:{rownames(bin_counts_s)}, label: {unlist(group_labels[!is.na(group_labels)])}")

  group_labels <- group_labels[!is.na(group_labels)]
  group_labels <- data.frame(Group = group_labels)
  if (enrich > 0) {
    gene_colors <- rep(0, length(enriched_genes))
  } else {
    gene_colors <- c()
  }
  gene_colors[colnames(bin_counts_s) %in% genes_of_interest] <- 1
  gene_colors <- data.frame(Enriched_genes = gene_colors)

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
                               assay = "RNA",
                               slot = "data",
                               genes_of_interest = NULL,
                               ident = NULL,
                               reduction = 'tsne',
                               dimension = 1,
                               ...){
  if (is.null(ident)){
    ident <- Idents(object)
  } else {
    ident_df <- FetchData(object, vars = ident)
    ident <- ident_df[,1]
    names(ident) <- rownames(ident_df)
  }
  if (is.null(genes_of_interest)){
    stop("You must provide a set of genes to map.")
  }

  exprs <- GetAssayData(object = object, assay = assay, slot = slot) %>% as.matrix()
  if(reduction %nin% names(object)){
    stop(glue("{reduction} has not been performed on this object."))
  }
  dim_embed <- Embeddings(object = object,
                           reduction = reduction)[,dimension]
  ReducedHeatmap.default(expression_matrix = exprs,
                 genes_of_interest = genes_of_interest,
                 dim_embedding = dim_embed,
                 cell_labels = ident,
                 ...)
  }
