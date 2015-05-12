run.DESeq2 <- function(dds, annot, max.padj=0.05, threads=6) {
    dds <- DESeq(dds, parallel=T, BPPARAM=MulticoreParam(workers=threads))
    res <- results(dds)
    retval <- summarize.DESeq2.result(res, annot, max.padj)
    retval$dds <- dds
    retval
}

summarize.DESeq2.result <- function(res, annot, max.padj=0.05) {
    # remove rows with a p-value of NA
    keep <- !is.na(res$padj)
    retval <- list(
        result=res,
        keep=keep,
        sig=keep & res$padj < max.padj,
        upreg=keep & res$log2FoldChange > 0,
        downreg=keep & res$log2FoldChange < 0,
        max.padj=max.padj
    )
    retval$num.genes <- sum(retval$keep)
    retval$num.sig <- sum(retval$sig)
    retval$num.upreg <- sum(retval$sig & retval$upreg)
    retval$num.downreg <- sum(retval$sig & retval$downreg)
    retval$pct.sig <- retval$num.sig / retval$num.genes
    retval$pct.upreg <- retval$num.upgreg / retval$num.genes
    retval$pct.downreg <- retval$num.downreg / retval$num.genes
    res.annotated <- cbind(res, annot)
    retval$upreg.table <- res.annotated[retval$sig & retval$upreg,]
    retval$upreg.table <- retval$upreg.table[order(retval$upreg.table$padj, 1000 - retval$upreg.table$log2FoldChange),]
    rownames(retval$upreg.table) <- NULL
    retval$downreg.table <- res.annotated[retval$sig & retval$downreg,]
    retval$downreg.table <- retval$downreg.table[order(retval$downreg.table$padj, retval$downreg.table$log2FoldChange),]
    rownames(retval$downreg.table) <- NULL
    retval
}
