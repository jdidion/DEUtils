#' Perform differential expression analysis using sleuth.
#' requires a data table with one row per sample. The only
#' required column is `sample`, which gives the sample name.
#' There should also be at least one additional column that
#' specifies a condition that divides the samples into at
#' at least two groups. There can also be a `path` column to
#' specify that directory that contains the kallisto output 
#' for the sample, otherwise the path is assumed to be
#' <data.dir>/<sample>.
#'
#' The default values for the full model and beta value
#' (which you should almost certainly change) assume 
#' that the model is testing a single response variable 
#' (condition) assumed to have two possible values 
#' ('case' and 'control').
#'
#' @param samples data.frame, data.table, or path to a csv,
#' tsv, or Rdata file containing the data table.
#' @param models model(s), specified as formula. If a single
#' value, taken to be the full model, otherwise a list of
#' model names and formulas. If a model is named 'full' that
#' will be used as the full model, otherwise the first model
#' in the list.
#' @param beta response variable(s) to test for differential
#' expression. Each is a concatanation of the variable name
#' and the factor relative to which beta values should be
#' calculated. If a vector, each beta will be tested for each
#' model; otherwise, a list mapping model name to beta.
#' @param data.dir parent directory containing kallisto 
#' output directories.
#' @param outfile RData file in which to save final sleuth
#' object. Defaults to "sleuth_results.RData" in `data.dir`.
#' @param db data.frame with a target_id column containing
#' the transcript IDs used in the Kallisto index and mapping
#' them to other features (e.g. gene IDs/names). Uses the
#' bioMart Ensembl database by default. Use NA to avoid using
#' transript->gene mappings.
run.sleuth <- function(samples, models=~condition, betas="conditioncontrol", 
                       data.dir=".", outfile=file.path(data.dir, "sleuth_result.RData"), 
                       db=NULL, max.bootstrap=NULL) {
    if (is.character(samples)) {
        ext <- file.ext(samples) == 'RData'
        if (ext) {
            samples <- my.load(samples)
        }
        else {
            sep <- ifelse(ext == "csv", ",", "\t")
            read.table(samples, sep=sep, header=TRUE, stringsAsFactors=FALSE)
        }
    }
    if (!('path' %in% colnames(samples))) {
        samples$path <- file.path(data.dir, samples$sample)
    }
    
    if (length(models) == 1) {
        names(models) <- "full"
    }
    else {
        n <- names(models)
        if (is.null(n)) {
            n <- c("full", paste0("partial", 1:(length(models)-1)))
        }
        else if (!("full" %in% names(models))) {
            n[1] <- "full"
        }
        names(models) <- n
    }
    
    if (is.null(db)) {
        message("Loading Ensembl tx->gene mapping")
        db <- get.ensembl.db()
    }
    
    message("Preparing data")
    data <- sleuth::sleuth_prep(samples, full.model, target_mapping=db, max_bootstrap=max.bootstrap)
    
    for (mod in names(partial.models)) {
        message(paste("Fitting model", mod))
        data <- sleuth::sleuth_fit(data, formula=models[[mod]], fit_name=mod)
    }
    
    if (is.list(betas)) {
        for (mod in names(betas)) {
            for (b in betas[[mod]]) {
                message(paste("Testing", b, "for model", mod))
                data <- sleuth::sleuth_wt(data, which_beta=beta, which_mod=mod)
            }
        }
    }
    else {
        for (mod in names(models)) {
            for (b in betas) {
                message(paste("Testing", b, "for model", mod))
                data <- sleuth::sleuth_wt(data, which_beta=beta, which_mod=mod)
            }
        }
    }
    
    if (!is.null(outfile)) {
        save(data, file=outfile)
    }
    
    invisible(data)
}

get.ensembl.db <- function() {
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
        "external_gene_name"), mart = mart)
    t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
        ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    t2g
}

cluster <- function(sleuth, ..., group.vars=NULL) {
    mat <- sleuth_to_matrix(sleuth, "obs_norm", "tpm")
    mat <- mat$data[sleuth$filter_bool,]
    s <- sleuth$sample_to_covariates
    if (!is.null(group.vars)) {
        s <- s[,group.vars]
    }
    d <- call.with(feature.distance, mat=mat, ...)
    call.with(plot.mds, d=d, group.df=group.vars, ...)
}
