topGO <- function(p) {
    l <- lapply(c("BP","MF","CC"), function(onto) topGOtest(p, onto))
    names(l) <- c("BP","MF","CC")
    l
}

topGO.test <- function(p, onto, selection.cutoff=0.01, p.value.cutoff=0.05) {
    onto <- match.arg(onto, c("BP","MF","CC"))
    GOdata <- new(Class='topGOdata', ontology=onto, allGenes=p, 
        annot=annFUN.org, mapping='org.Hs.eg.db',
        ID="ensembl", geneSelectionFun=function(p) p < selection.cutoff)
    resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
    results.table <- GenTable(GOdata, classicFisher=resultFisher, topNodes=length(resultFisher@score))
    results.table$FDR <- p.adjust(results.table[,"classicFisher"],method="BH")
    results.table[which(results.table$FDR <= p.value.cutoff),]
}

# Tests for GO enrichment of an arbitrary set of genes
goseq.test.geneset <- function(genes, all.genes) {
    p <- rep(1, length(all.genes))
    names(p) <- all.genes
    p[genes] <- 0
    goseq.test(p)
}

goseq.test.DESeq2 <- function(res, annot, p.value.threshold=0.05) {
    w <- !is.na(res$padj)
    res <- res[w,]
    annot <- annot[w,]
    upreg <- res$log2FoldChange > 0
    downreg <- res$log2FoldChange < 0
    
    p <- res$padj
    names(p) <- annot$GeneID
    
    upreg.p <- p
    upreg.p[!upreg] <- 1
    
    downreg.p <- p
    downreg.p[!downreg] <- 1
    
    list(upreg=goseq.test(upreg.p, p.value.threshold), downreg=goseq.test(downreg.p, p.value.threshold))
}

goseq.test <- function(p, p.value.threshold=0.05) {
    genes <- as.integer(p <= p.value.threshold)
    names(genes) <- names(p)
    pwf <- nullp(genes,"hg19","ensGene")
    GO.wall <- goseq(pwf,"hg19","ensGene")
    enriched <- GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH") < p.value.threshold,]
    enriched <- enriched[order(enriched$over_represented_pvalue),]
    rownames(enriched) <- NULL
    missing <- GO.wall[p.adjust(GO.wall$under_represented_pvalue, method="BH") < p.value.threshold,]
    missing <- missing[order(missing$under_represented_pvalue),]
    rownames(missing) <- NULL
    list(enriched=enriched, missing=missing)
}