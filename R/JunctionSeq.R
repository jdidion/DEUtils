library(JunctionSeq)
# Run each step of the analysis separately, and
# save the results of each step.
run.JunctionSeq <- function(sample.files, sample.names, condition, flat.gff.file, nCores,
        data.file="jscs.RData", table.prefix="result") {

    analysis.type <- "junctionsAndExons"
    use.junctions <- TRUE
    use.exons <- TRUE
    use.novel.junctions <- TRUE
    meanCountTestableThreshold <- 7.5
    test.formula0 <- formula(~sample + countbin)
    test.formula1 <- formula(~sample + countbin + condition:countbin)
    effect.formula <- formula(~condition + countbin + condition:countbin)
    geneLevel.formula <- formula(~condition)
    gzip.output <- TRUE
    test.aggregated.genes <- TRUE
    fitDispersionsForExonsAndJunctionsSeparately <- TRUE
    use.alternate.method <- TRUE
    keep.debug.model.data <- TRUE
    gtf.format <- TRUE
    verbose <- TRUE
            
    if (verbose) {
        message("> STARTING runJunctionSeqAnalyses (", date(),
            ")")
        message(paste("> rJSA: sample.files: ", paste0(sample.files,
            collapse = ", ")))
        message(paste("> rJSA: sample.names: ", paste0(sample.names,
            collapse = ", ")))
        message(paste("> rJSA: condition: ", paste0(condition,
            collapse = ", ")))
        message(paste("> rJSA: analysis.type: ", analysis.type))
        message(paste("> rJSA: use.junctions: ", use.junctions))
        message(paste("> rJSA: use.novel.junctions: ", use.novel.junctions))
        message(paste("> rJSA: use.exons: ", use.exons))
        message(paste("> rJSA: nCores: ", nCores))
        message(paste("> rJSA: test.formula0: ", paste0(test.formula0,
            collapse = " ")))
        message(paste("> rJSA: test.formula1: ", paste0(test.formula1,
            collapse = " ")))
        message(paste("> rJSA: gzip.output: ", gzip.output))
        message(paste("> rJSA: use.alternate.method: ", use.alternate.method))
        message(paste("> rJSA: test.aggregated.genes: ", test.aggregated.genes))
    }
    design <- data.frame(condition = condition)
    row.names(design) <- sample.names

    if (verbose) {
        message(paste0("> rJSA: Reading Count files...", " ",
            date(), "."))
    }

    jscs <- readJunctionSeqCounts(countfiles = as.character(sample.files),
            samplenames = sample.names, design = design, flat.gff.file = flat.gff.file,
            verbose = verbose, use.junctions = use.junctions, 
            use.novel.junctions = use.novel.junctions,
            use.exons = use.exons, nCores = nCores)
    save(jscs, file=data.file)

    if (verbose) {
        message(paste0("> rJSA: Count files read.", " ", date(),
            "."))
        message(paste0("> rJSA: Estimating Size Factors...",
            " ", date(), "."))
    }

    jscs <- estimateSizeFactors(jscs)
    save(jscs, file=data.file)

    if (verbose) {
        message(paste0("> rJSA: Size Factors Done. Size Factors are:",
            "."))
        message(paste0("> rJSA: ", paste0(names(sizeFactors(jscs)),
            collapse = ",")))
        message(paste0("> rJSA: ", paste0(sizeFactors(jscs),
            collapse = ",")))
        message(paste0("> rJSA: Estimating Dispersions...", " ",
            date(), "."))
    }

    jscs <- estimateJunctionSeqDispersions(jscs, nCores = nCores,
        test.formula1 = test.formula1, meanCountTestableThreshold = meanCountTestableThreshold,
        use.alternate.method = use.alternate.method, test.aggregated.genes = test.aggregated.genes,
        verbose = verbose)
    save(jscs, file=data.file)

    if (verbose) {
        message(paste0("> rJSA: Dispersions estimated.", " ",
            date(), "."))
        message(paste0("> rJSA: Fitting Dispersion Fcn...", " ",
            date(), "."))
    }

    jscs <- fitDispersionFunction(jscs, verbose = verbose, fitDispersionsForExonsAndJunctionsSeparately = fitDispersionsForExonsAndJunctionsSeparately)
    save(jscs, file=data.file)

    if (verbose) {
        message(paste0("> rJSA: Dispersions Fcn Fitted.", " ",
            date(), "."))
        message(paste0("> rJSA: Testing for DEU...", " ", date(),
            "."))
    }

    jscs <- testForDiffUsage(jscs, nCores = nCores, test.formula0 = test.formula0,
        test.formula1 = test.formula1, use.alternate.method = use.alternate.method,
        verbose = verbose)
        save(jscs, file=data.file)

    if (verbose) {
        message(paste0("> rJSA: DEU tests complete.", " ", date(),
            "."))
        message(paste0("> rJSA: Estimating effect sizes using effects models...",
            " ", date(), "."))
    }

    jscs <- estimateEffectSizes(jscs, effect.formula = effect.formula,
        geneLevel.formula = geneLevel.formula, nCores = nCores)
    sampleNames(jscs) = as.character(sample.names)
    save(jscs, file=data.file)
    
    if (verbose) {
        message(paste0("> rJSA: Effect Sizes estimated.", "."))
        message(paste0("> rJSA: Generating results table...",
            " ", date(), "."))
    }

    writeCompleteResults(jscs, outfile.prefix=table.prefix, save.jscs = TRUE)
}
