% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReducedHeatmap.R
\name{ReducedHeatmap}
\alias{ReducedHeatmap}
\alias{ReducedHeatmap.default}
\alias{ReducedHeatmap.Seurat}
\alias{ReducedHeatmap.SingleCellExperiment}
\title{ReducedHeat}
\usage{
ReducedHeatmap(object, ...)

\method{ReducedHeatmap}{default}(object, features, dim_embedding,
  cell_labels, enrich = 0, breaks = 100, slope = 50,
  intercept = 0.05, dendrogram = "row", titleX = FALSE,
  RowV = TRUE, show_legend = FALSE, hide_colorbar = FALSE,
  fontsize_row = 7, margins = c(70, 50, NA, 0), ...)

\method{ReducedHeatmap}{Seurat}(object, assay = NULL, slot = "data",
  features = NULL, ident = NULL, reduction = "tsne", dimension = 1,
  ...)

\method{ReducedHeatmap}{SingleCellExperiment}(object,
  assay = "logcounts", features = NULL, ident = NULL,
  reduction = "tsne", dimension = 1, ...)
}
\arguments{
\item{object}{Object or expression matrix to visualize. Rownames should be gene names}

\item{...}{Additional parameters to pass to heatmaply}

\item{features}{Features to plot}

\item{dim_embedding}{Cell embedding within the dimension of interest}

\item{cell_labels}{Labels for each cell. While this is generally the cluster identity for a cell, it can be any discrete metadata variable}

\item{enrich}{For every item in \code{features}, find a number of features equal to \code{enrich} that are close to that feature in the distance induced by the 1D UMAP}

\item{breaks}{Number of bins}

\item{slope}{For better visualization, transform the values with a logistic function. This is the slope of that function.}

\item{intercept}{For better visualization, transform the values with a logistic function. This is the intercept of that function.}

\item{dendrogram}{See heatmaply help. Default: "row"}

\item{titleX}{See heatmaply help. Default: FALSE.}

\item{RowV}{See heatmaply help. Default: TRUE}

\item{show_legend}{See heatmaply help. Default: FALSE}

\item{hide_colorbar}{See heatmaply help. Default: FALSE}

\item{fontsize_row}{See heatmaply help. Default: 7}

\item{margins}{See heatmaply help. Default: c(70, 50, NA, 0)}

\item{assay}{Assay to use. Default: the object's current default assay.}

\item{slot}{Slot to use. Default: "data"}

\item{ident}{Object variable to use when determining cell_labels}

\item{reduction}{Reduction to visualize. Default: 'umap'}

\item{dimension}{The dimension to visualize. Default: 1}
}
\description{
Compute Heatmap using dimensional reduction
}
