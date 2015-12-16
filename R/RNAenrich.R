library(annotate)
library(stats)
library(mgcv)
library(rms)
library(utils)

###########################################################################################
#  Function for RNA-Enrich: A cut-off free functional enrichment testing method 
#  for RNA-seq with improved detection power
#  Written by: Chee Lee, University of Michigan, 2015
###########################################################################################
#'
#'  This function uses modified random sets method with weights calculated from a spline to 
#'  adjust for any relationship between read count and differential expression p-value to test for 
#'  enriched (or depleted) biological categories in gene expression data. 
#'
#'  Please acknowledge your use of LRpath in publications by referencing:
#'  Lee C, Patil S, Sartor MA. RNA-Enrich: A cut-off free functional enrichment testing 
#'  method for RNA-seq with improved detection power. 
#'
#'  Inputs:
#'  data :: data.frame with three or four columns: geneid (Entrez gene IDs), p.value (p-values), 
#'  readcount (), and (optionally) direction (direction of effect for directional testing
#'  values should be negative if down-regulated and positive if up-regulated (e.g. -1 or 1). 
#'  Only the sign is used).
#'  species :: Used only for KEGG, GO, and Cytoband human="hsa", mouse="mmu", rat="rno", 
#'  "dme" for fly, "dre" for zebrafish, "cel" for C.elegans, "sce" for yeast etc. (no default)
#'  min.g :: The minimum number of unique gene IDs analyzed in category to be tested (default = 10)
#'  max.g :: The maximum number of unique gene IDs analyzed in category to be tested
#'  (default = NA (99999))
#'  sig.cutoff" Entrez gene IDs in each category with p-values<sig.cutoff will be returned (default
#'  = 0.05)
#'  database :: database to be tested - choices are "GO", "KEGG", "Cytoband", "custom"
#'  (default "GOBP"). "custom" is for user provided gene sets. When "custom" is selected, 
#'  conceptList is a list of custom gene set IDs, nullsetList is list of genes connected to
#'  conceptList, should be a 1 gene set ID to 1 gene ID relationship between conceptList and
#'  nullsetList.
#'  read_lim :: Read limit of genes that will be included in analysis (default = 5)
#'  conceptList :: List of concepts or gene sets, use only when database="custom"
#'  nullsetList :: List of genes, use only when database="custom." conceptList and nullsetList
#'  should be a 1 gene set ID to 1 gene ID relationship between conceptList and nullsetList.
#'  plot_file :: Name of file of output plot (default = "RNA-Enrich_plot.jpg")
#'  plot_height :: Height in pixels of plot (default = 480 pixels)
#'  plot_width :: Width in pixels of plot (default = 480 pixels)
#'  results_file :: Name of file of tab-delimited text output of results (default =
#'  "RNA-Enrich_results.txt")
#'
#'  Outputs:
#'  object is dataframe with the following columns:
#'  1) GO, KEGG, or Concept ID  
#'  2) GO, KEGG, or Concept term - name of category
#'  3) database - Database of gene sets that was tested
#'  4) n.genes - Number of unique Entrez Gene IDs in category
#'  5) coeff - coefficient of slope (positive values indicate enrichment or up-regulation 
#'  negative values indicate depletion or down-regulation)
#'  6) odds.ratio -  Odds ratio, as measure of strenth of enrichment (or depletion)
#'  7) status - "Enriched" or "Depleted", or if direction is used "Up-regulated" or "Down-regulated"
#'  8) p.value  -  P-value that slope does not equal zero, i.e. that term is enriched (or depleted)
#'  9) FDR  - False Discovery Rate (Benjamini & Hochberg, 1995)
#'  10) sig.genes   - comma separated Entrez gene ids in category with p-value<"sig.cutoff"
#'
#' @ex
#'    pvalues <- IBMT.results$IBMT.p[,1]
#'    entrez.geneid <- unlist(mget(featureNames(affy.eset),env=hgu133aENTREZID))
#'    GO.results <- LRpath(sigvals=pvalues, geneids=entrez.geneid, species="hsa", database="GO")
#'    KEGG.results <- LRpath(sigvals=pvalues, geneids=entrez.geneid, database="KEGG",
#'      species="hsa")
#'
##########################################################################################
rna_enrich <- function(data, species, min.g=10, max.g=NA, sig.cutoff=0.05, database="GOBP", 
                       read_lim=5, conceptList=NULL, nullsetList=NULL, 
                       results_file="RNA-Enrich_results.txt") {

    has.direction <- ncol(data) > 3
    
    if (has.direction) {
        data <- stats::aggregate(cbind(p.value, readcount, direction) ~ geneid, data, mean) 
        data$direction <- ifelse(data$direction > 0, 1, -1)
    }
    else {  
        data <- stats::aggregate(cbind(p.value, readcount) ~ geneid, data, mean)
    }
    
    data <- na.omit(data)
    data <- subset(data, geneid !='')
    
    if (!species %in% c('dme','sce')) {
        data$geneid <- suppressWarnings(as.numeric(as.character(data$geneid)))
        data <- na.omit(data)
    }
    
    # Get entire gene list using GO.db, for RNA-Enrich, must test each branch of 
    # GO separately for later weighting step to work correctly
    if (database %in% c('GOBP','GOCC','GOMF')) {
        library("GO.db")
        
        GO.data <- list(
            hsa=list(lib="org.Hs.eg.db", db="org.Hs.egGO2ALLEGS"),
            mmu=list(lib="org.Mm.eg.db", db="org.Mm.egGO2ALLEGS"),
            rno=list(lib="org.Rn.eg.db", db="org.Rn.egGO2ALLEGS"),
            dme=list(lib="org.Dm.eg.db", db="org.Dm.egGO2ALLEGS"),
            dre=list(lib="org.Dr.eg.db", db="org.Dr.egGO2ALLEGS"),
            cel=list(lib="org.Ce.eg.db", db="org.Ce.egGO2ALLEGS"),
            sce=list(lib="org.Sc.sgd.db", db="org.Sc.sgdGO2ALLORFS")
        )[[species]]
        
        library(GO.data$lib)
        genes <- as.list(get(gGO.data$db))
        ontols <- Ontology(names(genes))
        
        if (database == 'GOBP') {
            genes <- genes[ontols=='BP']
        }
        else if (database == 'GOCC') {
            genes <- genes[ontols=='CC']
        }
        else if (database == 'GOMF') {
            genes <- genes[ontols=='MF']
        }
    
        ##Remove GO ids that aren't mapped to any Gene id
        genes <- genes[!is.na(genes)]

        ## Determine which Entrez Gene IDs are annotated anywhere in GO
        ## Gene ids not annotated in GO are left out of the environment
        if (species %in% c("sce")) {   
            nullset <- unique(unlist(genes))
        } 
        else {
            nullset <- as.numeric(unique(unlist(genes)))
        }
        
        concept_descriptions <- Term(names(genes))
    }
    #########  KEGG works for human, mouse, rat, and any other species that uses Entrez geneID
    else if (database=="KEGG") {
        library(KEGG.db)
        genes <- as.list(KEGGPATHID2EXTID)
        genes <- genes[!is.na(genes)]
        genes <- genes[grep(species, names(genes))]

        ## Determine which Entrez Gene IDs are annotated anywhere in KEGG
        if (species %in% c("dme","sce")) {
           nullset <- names(as.list(KEGGEXTID2PATHID))
        }
        else {
            nullset <- as.numeric(names(as.list(KEGGEXTID2PATHID)))
        }
        nullset <- nullset[!is.na(nullset)]
        
        KEGGterms <- as.list(KEGGPATHNAME2ID)
        kterms <- as.vector(names(KEGGterms))
        k.ids <- paste(species, unlist(KEGGterms), sep="")
        keggrows <- match(names(genes), k.ids)
        concept_descriptions <- kterms[keggrows]
    }
    #######  For this, we'll use data in species specific R packages, e.g. org.Hs.eg.db for human
    else if (database=="Cytoband") {
        if (species != "hsa") {
            stop("Only human is supported in Cytoband")
        }
        
        library("org.Hs.eg.db")
        ###  Get GO list of Entrez genes
        temp <- as.list(org.Hs.egMAP)
        ## Currently, list names are Entrez IDs, and values are cytobands.  Need to switch.
        temp <- temp[!is.na(temp)]
        ## names(xx1) are the entrez IDs to match this
        cyto <- as.vector(sapply(temp, "[", 1))
        bands <- unique(cyto)
        genes <- lapply(bands, function(x) names(temp)[cyto == x])
        names(genes) <- as.character(bands)
        ## Get null set vector for Cytoband
        nullset <- as.numeric(names(temp))
        nullset <- nullset[!is.na(nullset)]
        concept_descriptions <- NA
    }
    else if (database=="custom"){
        if (is.null(conceptList) || is.null(nullsetList)) {
            stop("Need conceptList and nullsetList to proceed!")
        }
        genes <- by(nullsetList, conceptList, list)
        genes <- genes[!is.na(genes)]
        nullset <- unique(as.numeric(unlist(genes)))
        concept_descriptions <- names(genes)
    }
    
    catsizes <- sapply(genes, length)
    if (is.na(max.g)) { 
        max.g <- max(catsizes)
    }
    genes <- genes[catsizes >= min.g & catsizes <= max.g]
    ncats <- length(genes)

    #make empty d dataframe to initiate gpw creation
    allgenes <- data.frame(
        geneid=nullset[which(!nullset %in% data$geneid)],
        p.value=NA, readcount=0
    )
    if (has.direction) {
        allgenes$direction <- NA
        allgenes$direction <- as.numeric(allgenes$direction)
    }
    gpw <- rbind(data, allgenes)

    # Add log10 p-values to gpw, add min p[!=0] to all pvals so can log transform later
    gpw$log10_p.value <- (-1) * log10(gpw$p.value + min(gpw$p.value[gpw$p.value !=0], na.rm=T))

    # remove genes with read lim readcounts (default is 5) and log readcount
    if (read_lim==0) { 
        gpw <- subset(gpw, readcount > read_lim)
    }
    else {
        gpw <- subset(gpw, readcount >= read_lim)
    }
    gpw$log10_readcount <- log10(gpw$readcount)

    # Sort by readcount 
    gpw <- gpw[order(gpw$readcount, decreasing=T),]
    gpw <- gpw[!is.na(gpw$log10_p.value),]
    
    # Create model - adjust for readcount
    model <- "log10_p.value ~ s(log10_readcount,bs='cr')"

    # Compute binomial spline fit.
    fit <- gam(as.formula(model), data=gpw, family="gaussian")

    # Compute weights for each gene, based on the predicted prob(DE) for each gene. 
    pP <- fitted(fit)
    w0 <- 1 / (pP/mean(gpw$log10_p.value, na.rm=T))
    w0 <- w0 / mean(w0, na.rm=T)

    gpw$weight <- w0
    gpw$prob_P <- pP
    gpw$resid.dev <- resid(fit, type="deviance")

    # change log pvals if directional test
    if (has.direction) {
        gpw$log10_p.value <- gpw$log10_p.value * gpw$direction
    }

    gpw <- subset(gpw, select=c(
        "geneid","readcount","log10_readcount","p.value",
        "log10_p.value","weight","prob_P","resid.dev"
    ))

    # now run test   

    lrm.fast <- function(x,y) {
        fit <- lrm.fit(x, y)
        vv <- diag(fit$var)
        cof <- fit$coef
        z <- cof / sqrt(vv)
        pval <- pchisq(z^2,1, lower.tail=F)
        c(cof[2], pval[2])
    }
    
    # Restrict our genes/weights/peaks to only those genes in the genesets. 
    gpw <- subset(gpw,geneid %in% nullset)
    
    # Re-normalize weights so that mean is 1 (if any genes were removed from previous gpw step)
    gpw$weight <- gpw$weight / mean(gpw$weight)
    
    N <- length(genes)
    
    result <- data.frame(
        "Concept ID"=names(genes),
        "Concept name"=concept_descriptions,
        "database"=database,
        "n.genes"=integer(N),
        "coeff"=numeric(N),
        "odds.ratio"=numeric(N),
        "status"=character(N),
        "p.value"=numeric(N),
        "FDR"=numeric(N),
        "sig.genes"=integer(N),
        stringsAsFactors=F
    )
    result$coeff <- NA
    result$p.value <- NA
    result$status <- NA
    
    pb <- txtProgressBar(min = 0, max = N, title="Running tests")
    
    # TODO: parallelize?
    
    for (i in 1:N) {
        setTxtProgressBar(pb, i)
    
        if (species %in% c('dme','sce')) {
            go_genes <- as.character(genes[[i]])
        } 
        else{
            go_genes <- as.numeric(as.character(genes[[i]]))
        }
        
        # Eliminate GO genes that aren't in the gpw. 
        go_genes <- go_genes[go_genes %in% gpw$geneid]
        result[i, "n.genes"] <- length(go_genes)
        
        # A boolean vector of gene membership of all geneids in go_genes
        b_genes <- gpw$geneid %in% go_genes
        sg_go <- gpw$log10_p.value[b_genes]
        wg_go <- gpw$weight[b_genes]
        
        go_genes_DE <- gpw$geneid[gpw$p.value<=sig.cutoff & b_genes==TRUE] 
        if (length(go_genes_DE)==0) {
            go_genes_DE=NA
        }
        result[i, "sig.genes"] <- paste(go_genes_DE, collapse=", ")
        
        testm <- cbind(y=as.numeric(b_genes), x=gpw$weight * gpw$log10_p.value)
        
        if (sum(testm[,'y']==1) != nrow(gpw) && sum(testm[,'y']) != 0) {
            result[i, c("coeff", "p.value")] <- lrm.fast(testm[,"x"], testm[,"y"])
        }
    }
    
    close(pb)
    
    result$FDR <- p.adjust(result$p.value, method="BH")
    result$odds.ration <- exp(result$coeff)
    
    is_depleted <- result$coeff < 0 
    is_enriched <- result$coeff > 0
    if (has.direction) {
        result$status[is_depleted] <- "down"
        result$status[is_enriched] <- "up"
    }
    else{
        result$status[is_depleted] <- "depleted"
        result$status[is_enriched] <- "enriched"
    }
    
    result <- result[!is.na(result$p.value),]
    result <- result[order(result$FDR, result$p.value),]
    
    write.table(result, file=results_file, sep="\t",row.names=F, quote=F)
    
    result
}

rna_enrich_plot <- function(data, filename="rna_enrich_plot.jpg", plot_height=480, plot_width=480,
                            min.readcount=5, group.size=25) {
    data <- na.omit(data)
    data <- data[data$readcount >= min.readcount,]
    data <- data[order(data$p.value),]
        
    rownames(data) <- NULL
    data$group <- ceiling(as.numeric(rownames(data)) / group.size)
    gavgrc <- as.numeric(by(data$readcount, data$group, mean))
    gavgp <- as.numeric(by(data$p.value, data$group, mean))
    log10_gavgrc <- log10(gavgrc)
    log10_gavgp <- -log10(gavgp + min(gavgp[gavgp!=0]))

    data$log10_p.value <- (-1)*log10(data$p.value + min(data$p.value[data$p.value !=0],na.rm=T))
    data$log10_readcount <- log10(data$readcount)
    data <- data[!is.na(data$log10_p.value),]
    model <- "log10_p.value ~ s(log10_readcount,bs='cr')"
    
    # Compute binomial spline fit.  
    fit <- gam(as.formula(model),data=data,family="gaussian")
    # Compute weights for each gene, based on the predicted prob(peak) for each gene. 
    pP <- fitted(fit)
     
    jpeg(filename, height=plot_height, width=plot_width)
    par(mar= c(5, 5, 4, 2) + 0.1)
      
    plot(log10_gavgp ~ log10_gavgrc, pch=20, cex.axis=1.5, cex.lab=1.5
        ylab='-log10(p-value) - binned 25 genes', 
        xlab='log10(avg read count) - binned 25 genes')
    
    lines(pP[order(data$log10_readcount)]~sort(data$log10_readcount),col='red',lwd=2)
    
    dev.off()
}
