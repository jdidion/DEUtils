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