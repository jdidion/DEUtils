parse.gff <- function(flattenedfile) {
    aggregates <- read.delim(flattenedfile, stringsAsFactors = FALSE, header = FALSE)
    colnames(aggregates) <- c("chr", "source", "class", "start",  "end", "ex", "strand", "ex2", "attr")
    aggregates$strand <- gsub("\\.", "*", aggregates$strand)
    aggregates <- aggregates[which(aggregates$class == "exonic_part"),]
    aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
    aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", aggregates$attr)
    transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", aggregates$attr)
    transcripts <- strsplit(transcripts, "\\+")
    exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", aggregates$attr)
    exoninfo <- GRanges(as.character(aggregates$chr), 
                        IRanges(start = aggregates$start, 
                                end = aggregates$end), 
                        strand = aggregates$strand)
    names(exoninfo) <- paste(aggregates$gene_id, exonids, sep = ":E")
    names(transcripts) <- rownames(exoninfo)
    list(exoninfo=exoninfo, transcripts=transcripts)
}

run.DEXSeq <- function(counts, samples, full.model, reduced.model, fitExpToVar, gff, threads=6) {
    features <- sub(":", ":E", rownames(counts))
    rownames(counts) <- features
    splitted <- strsplit(features, ":")
    genesrle <- sapply(splitted, "[[", 1)
    exons <- sapply(splitted, "[[", 2)
    
    if (!all(rownames(counts) %in% names(gff$exoninfo))) {
        stop("Count files do not correspond to the flattened annotation file")
    }
    matching <- match(rownames(counts), names(gff$exoninfo))
    stopifnot(all(names(gff$exoninfo[matching]) == rownames(counts)))
    stopifnot(all(names(gff$transcripts[matching]) == rownames(counts)))
    
    dxd <- DEXSeqDataSet(counts, samples, full.model, exons, genesrle, gff$exoninfo[matching], gff$transcripts[matching])
    DEXSeq(dxd, reducedModel=reduced.model, BPPARAM=MulticoreParam(workers=threads), fitExpToVar=fitExpToVar)
}

my.plotDEXSeq <- function (object, geneID, FDR = 0.1, fitExpToVar = "condition",
                           norCounts = FALSE, expression = TRUE, splicing = FALSE, displayTranscripts = FALSE,
                           names = FALSE, legend = FALSE, color = NULL, color.samples = NULL, ...)
{
    stopifnot(is(object, "DEXSeqResults"))
    op <- sum(c(expression, splicing, norCounts))
    if (op == 0) {
        stop("Please indicate what would you like to plot\n")
    }
    sampleData <- object@sampleData
    genomicData <- object$genomicData
    rt <- which(object$groupID == geneID)
    count <- t(t(object$countData[rt, ])/sampleData$sizeFactor)
    if (sum(count) == 0) {
        warning("No read counts falling in this gene, there is nothing to plot.")
        return()
    }
    if (FDR > 1 | FDR < 0) {
        stop("FDR has to be a numeric value between 0 - 1")
    }
    rango <- seq(along = rt)
    intervals <- (0:nrow(count))/nrow(count)
    numcond <- length(unique(sampleData[[fitExpToVar]]))
    numexons <- nrow(count)
    each <- object$padj[rt]
    exoncol <- ifelse(each <= FDR, "#F219ED", "#CCCCCC")
    exoncol[is.na(exoncol)] <- "white"
    colorlines <- ifelse(each <= FDR, "#F219ED60", "#B3B3B360")
    colorlines[is.na(colorlines)] <- "#B3B3B360"
    colorlinesB <- ifelse(each <= FDR, "#9E109B", "#666666")
    colorlinesB[is.na(colorlinesB)] <- "#666666"
    if (length(start(unlist(genomicData))) > 0) {
        sub <- data.frame(start = start(genomicData[rt]), end = end(genomicData[rt]),
                          chr = as.character(seqnames(genomicData[rt])), strand = as.character(strand(genomicData[rt])))
        rel <- (data.frame(sub$start, sub$end)) - min(sub$start)
        rel <- rel/max(rel[, 2])
        if (displayTranscripts == TRUE & !is.null(unlist(object$transcripts[rt]))) {
            transcripts <- object$transcripts[rt]
            trans <- Reduce(union, transcripts)
            if (length(trans) > 40) {
                warning("This gene contains more than 40 transcripts annotated, only the first 40 will be plotted\n")
            }
            mat <- seq_len(3 + min(length(trans), 40))
            hei <- c(8, 1, 1.5, rep(1.5, min(length(trans), 40)))
        }
        else {
            mat <- 1:3
            hei <- c(5, 1, 1.5)
        }
        if (op > 1) {
            hei <- c(rep(hei[1], op - 1), hei)
            mat <- c(mat, length(mat) + seq(along = op))
        }
        hei <- c(hei, 0.2)
        mat <- c(mat, length(mat) + 1)
        layout(matrix(mat), heights = hei)
        par(mar = c(2, 4, 4, 2))
    }
    else if (op > 1) {
        par(mfrow = c(op, 1))
    }
    if (is.null(color)) {
        if (numcond < 10) {
            color <- suppressWarnings(brewer.pal(numcond, "Set1")[seq_len(numcond)])
        }
        else {
            color <- rgb(colorRamp(brewer.pal(5, "Set1"))(seq(0,
                                                              1, length.out = numcond)), maxColorValue = 255,
                         alpha = 175)
        }
    }
    names(color) <- sort(levels(sampleData[[fitExpToVar]]))
    mf <- object@modelFrameBM
    mf <- mf[as.vector(sapply(split(seq_len(nrow(mf)), mf$sample),
                              "[", seq_len(numexons))), ]
    featuresInGene <- object$featureID[rt]
    mf$exon <- factor(rep(featuresInGene, nrow(sampleData)))
    counts <- object$countData[rt, ]
    rownames(counts) <- gsub("\\S+:", "", rownames(counts))
    dispersions <- object$dispersion[rt]
    dispersions[is.na(dispersions)] <- 1e-08
    names(dispersions) <- object$featureID[rt]
    for (i in seq_len(nrow(mf))) {
        mf[i, "dispersion"] <- dispersions[as.character(mf[i,
                                                           "exon"])]
        mf[i, "count"] <- counts[as.character(mf[i, "exon"]),
                                 as.character(mf[i, "sample"])]
    }
    mf <- droplevels(mf)
    if (expression) {
        es <- DEXSeq:::fitAndArrangeCoefs(frm = as.formula(paste("count ~",
                                                                 fitExpToVar, "* exon")), balanceExons = TRUE, mf = mf, maxRowsMF=20000)
        if (is.null(es)) {
            warning(sprintf("glm fit failed for gene %s", geneID))
            return()
        }
        coeff <- as.matrix(t(DEXSeq:::getEffectsForPlotting(es, averageOutExpression = FALSE,
                                                            groupingVar = fitExpToVar))[featuresInGene, ])
        coeff <- exp(coeff)
        ylimn <- c(0, max(coeff, na.rm = TRUE))
        coeff <- DEXSeq:::vst(coeff, object)
        DEXSeq:::drawPlot(matr = coeff, ylimn, object, intervals, rango,
                          textAxis = "Expression", rt = rt, color = rep(color[colnames(coeff)],
                                                                        each = numexons), colorlines = colorlines, ...)
    }
    if (splicing) {
        es <- DEXSeq:::fitAndArrangeCoefs(frm = as.formula(paste("count ~",
                                                                 fitExpToVar, "* exon")), balanceExons = TRUE, mf = mf, maxRowsMF=20000)
        if (is.null(es)) {
            warning(sprintf("glm fit failed for gene %s", geneID))
            return()
        }
        coeff <- as.matrix(t(DEXSeq:::getEffectsForPlotting(es, averageOutExpression = TRUE,
                                                            groupingVar = fitExpToVar))[featuresInGene, ])
        coeff <- exp(coeff)
        ylimn <- c(0, max(coeff, na.rm = TRUE))
        coeff <- DEXSeq:::vst(coeff, object)
        DEXSeq:::drawPlot(matr = coeff, ylimn, object, intervals, rango,
                          textAxis = "Exon usage", rt = rt, color = rep(color[colnames(coeff)],
                                                                        each = numexons), colorlines = colorlines, ...)
    }
    if (norCounts) {
        ylimn <- c(0, max(count, na.rm = TRUE))
        count <- DEXSeq:::vst(count, object)
        if (is.null(color.samples)) {
            colorcounts <- rep(color[as.character(sampleData[[fitExpToVar]])],
                               each = numexons)
        }
        else {
            colorcounts <- rep(color.samples, each = numexons)
        }
        DEXSeq:::drawPlot(matr = count, ylimn, object, intervals, rango,
                          textAxis = "Normalized counts", rt = rt, color = colorcounts,
                          colorlines = colorlines, ...)
    }
    if (length(start(unlist(genomicData))) > 0) {
        par(mar = c(0, 4, 0, 2))
        plot.new()
        segments(apply((rbind(rel[rango, 2], rel[rango, 1])),
                       2, median), 0, apply(rbind(intervals[rango], intervals[rango +
                                                                                  1] - ((intervals[rango + 1] - intervals[rango]) *
                                                                                            0.2)), 2, median), 1, col = colorlinesB)
        par(mar = c(1.5, 4, 0, 2))
        DEXSeq:::drawGene(min(sub$start), max(sub$end), tr = sub, rango,
                          exoncol = exoncol, names, trName = "Gene model",
                          cex = 0.8)
        if (length(unlist(object$transcripts[rt])) > 0) {
            if (displayTranscripts) {
                for (i in seq_len(min(length(trans), 40))) {
                    logicexons <- sapply(transcripts, function(x) {
                        length(which(x == trans[i]))
                    })
                    tr <- as.data.frame(reduce(IRanges(sub$start[logicexons ==
                                                                     1], sub$end[logicexons == 1])))[, c("start",
                                                                                                         "end")]
                    DEXSeq:::drawGene(min(sub$start), max(sub$end), tr = tr,
                                      rango, exoncol = NULL, names, trName = trans[i],
                                      cex = 0.8)
                }
            }
        }
        axis(1, at = round(seq(min(sub$start), max(sub$end),
                               length.out = 10)), labels = round(seq(min(sub$start),
                                                                     max(sub$end), length.out = 10)), pos = 0, lwd.ticks = 0.2,
             padj = -0.7, ...)
    }
    if (legend) {
        mtext(paste(geneID, unique(sub$strand)), side = 3, adj = 0.25,
              padj = 1.5, line = 0, outer = TRUE, cex = 1.5)
        posforlegend <- seq(0.7, 0.9, length.out = numcond)
        for (i in seq(along = color)) {
            mtext(names(color[i]), side = 3, adj = posforlegend[i],
                  padj = 1.5, line = 0, outer = TRUE, col = color[i],
                  ...)
        }
    }
    else {
        mtext(paste(geneID, unique(sub$strand)), side = 3, adj = 0.5,
              padj = 1.5, line = 0, outer = TRUE, cex = 1.5)
    }
}
