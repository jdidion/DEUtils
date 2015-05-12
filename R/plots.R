plot.heatmap.DESeq <- function(dds, sample.names=NULL) {
  if (is.null(sample.names)) {
    sample.names <- rownames(dds)
  }
  sampleDists <- dist( t( assay(dds) ) )
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- sample.names
  hc <- hclust(sampleDists)
  heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=rev(white.blue.pal),
           margins=c(2,10), labCol=FALSE)
}

plot.pca.DESeq <- function(dds, PCs=c(1,2), intgroup="condition", ntop=500, 
                           returnData=FALSE, plotData=TRUE, returnFig=FALSE, ...) {
    if (!all(intgroup %in% names(colData(dds)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    rv <- rowVars(assay(dds))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(dds)[select, ]))
    pct.var <- percent.var(pca)
    intgroup.df <- as.data.frame(colData(dds)[, intgroup, drop = FALSE], row.names=colnames(dds))
    if (plotData) {
        fig <- plot.pca(pca, intgroup.df, PCs, pct.var, plotData, ...)
    }
    if (returnData || returnFig) {
        result <- list()
        if (returnData) {
            result$data <- pca
            result$pct.var <- pct.var
        }
        if (returnFig) {
            result$fig <- fig
        }
        invisible(result)
    }
}
