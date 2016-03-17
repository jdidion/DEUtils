run.topGO.DESeq2 <- function(res) {
    p <- res$padj
    names(p) <- substr(rownames(res), 1, 15)
    run.topGO(p)
}

run.topGO <- function(p, ...) {
    l <- lapply(c("BP","MF","CC"), function(onto) topGO.test(p, onto, ...))
    names(l) <- c("BP","MF","CC")
    l
}

topGO.test <- function(p, onto, selection.cutoff=0.01, p.value.cutoff=0.05, db='org.Hs.eg.db', ID.type="ensembl") {
    onto <- match.arg(onto, c("BP","MF","CC"))
    GOdata <- new(Class='topGOdata', ontology=onto, allGenes=p, 
                  annot=annFUN.org, mapping=db,
                  ID=ID.type, geneSelectionFun=function(p) p < selection.cutoff)
    resultFisher <- runTest(GOdata, algorithm="classic", statistic="fisher")
    results.table <- GenTable(GOdata, classicFisher=resultFisher, topNodes=length(resultFisher@score))
    results.table$FDR <- p.adjust(results.table[,"classicFisher"],method="BH")
    retval <- results.table[which(results.table$FDR <= p.value.cutoff),]
    list(data=GOdata, results=retval)
}

# Tests for GO enrichment of an arbitrary set of genes
goseq.test.geneset <- function(genes, all.genes) {
    p <- rep(1, length(all.genes))
    names(p) <- all.genes
    p[genes] <- 0
    goseq.test(p)
}

goseq.test.DESeq2 <- function(res, annot, p.value.threshold=0.05, ...) {
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
    
    list(upreg=goseq.test(upreg.p, p.value.threshold, ...), downreg=goseq.test(downreg.p, p.value.threshold))
}

goseq.test <- function(p, p.value.threshold=0.05, db="hg19", ID.type="ensGene") {
    genes <- as.integer(p <= p.value.threshold)
    names(genes) <- names(p)
    pwf <- nullp(genes, db, ID.type)
    GO.wall <- goseq(pwf, db, ID.type)
    enriched <- GO.wall[p.adjust(GO.wall$over_represented_pvalue, method="BH") < p.value.threshold,]
    enriched <- enriched[order(enriched$over_represented_pvalue),]
    rownames(enriched) <- NULL
    missing <- GO.wall[p.adjust(GO.wall$under_represented_pvalue, method="BH") < p.value.threshold,]
    missing <- missing[order(missing$under_represented_pvalue),]
    rownames(missing) <- NULL
    list(enriched=enriched, missing=missing)
}

get.all.Ids <- function(db="Hs", type="ENSEMBL") {
    library(paste0("org.", db, ".eg.db"), character.only=TRUE)
    tab <- get(paste0("org.", db, ".eg", type))
    unique(unlist(as.list(tab[mappedkeys(tab)])))    
}

ids.to.entrez.map <- function(db="Hs", type="ENSEMBL") {
    library(paste0("org.", db, ".eg.db"), character.only=TRUE)
    tab <- get(paste0("org.", db, ".eg", type))
    map <- as.list(tab[mappedkeys(tab)])
    do.call(rbind, lapply(names(map), function(n) {
        data.frame(src.id=map[[n]], entrez.id=n)
    }))
}

ids.to.entrez <- function(ids, map, na.rm=T) {
    if (is.null(map)) {
        return(ids)
    }
    m <- match(ids, map$src.id)
    if (na.rm) {
        map[m[!is.na(m)], 'entrez.id']
    }
    else {
        map[m, 'entrez.id']
    }
}

# p - vector of p-values with names as gene names (can be in any format)
# map - a mapping of gene names to entrez IDs
dose.test <- function(p, map, p.value.threshold=0.05) {
    universe <- as.character(as.integer(as.character(ids.to.entrez(names(p), map))))
    ids <- universe[!is.na(p) & p <= p.value.threshold]
    summary(enrichDO(ids, universe=universe))
}
