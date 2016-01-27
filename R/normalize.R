# Copyright (c) 2015, Universitat Rovira i Virgili (Spain), University of 
# Aveiro (Portugal) & Aarhus University (Denmark)
# 
# Written by Carlos P. Roca
# as Research funded by the European Union
# for the research paper by Roca, Gomes, Amorim & Scott-Fordsmand: "Variation-
# preserving normalization unveils blind spots in gene expression profiling".
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# Normalize expression data using one of two condition-decomposition methods,
# median  condition-decomposition normalization or standard-vector 
# condition-decomposition normalization (see Roca et al for details).
#
# Args:
#     expr.data: matrix of expression data, with variables (probes or genes)
#                in rows and samples in columns
#     expr.cond: vector of sample conditions (character or factor). This
#                must be the same length as the number of columns in expr.data,
#                with NA values for all samples to be excluded.
#     method: "median" or "stdvec"
#
# Returns: list, which includes 'data' (the normalized expression matrix)
# and additional details of the normalization procedure.
condec.norm <- function(expr.data, expr.cond, method=c("median","stdvec")) {
    method <- match.arg(method)
    expr.data <- as.matrix(expr.data)
    stopifnot(length(expr.cond) == ncol(expr.data))
    if (method == "median") {
        normalize.median.condec(expr.data, expr.cond)
    }
    else {
        normalize.stdvec.condec(expr.data, expr.cond)
    }
}


# Implements median condition-decomposition normalization


normalize.median.condec <- function( expression.data, expression.condition, 
    norm.probability=0.5, verbose=FALSE, p.value.graph=NULL )
{
    # identify samples and conditions to normalize
    normalize.sample <- 
        colnames( expression.data )[ ! is.na( expression.condition ) ]
    normalize.sample.condition <- 
        as.character( na.omit( expression.condition ) )
    normalize.condition <- unique( normalize.sample.condition )
    
    if ( length( normalize.condition ) == 0 )
        stop( "No condition to normalize" )
    else if ( min( table( normalize.sample.condition ) ) < 2 )
        stop( "There must be 2 or more samples in each condition" )
    
    # select probes for normalization
    # no missing values in any normalization sample
    expression.probe <- rownames( expression.data )
    normalize.probe.idx <- which( rowSums( 
        is.na( expression.data[ , normalize.sample ] ) ) == 0 )

    # normalize within conditions

    normalize.expr.data <- matrix( nrow = length( expression.probe ), 
        ncol = length( normalize.sample ) )
    rownames( normalize.expr.data ) <- expression.probe
    colnames( normalize.expr.data ) <- normalize.sample
    
    normalize.within.cond.offset <- rep( 0, length( normalize.sample ) )
    names( normalize.within.cond.offset ) <- normalize.sample
    
    normalize.within.cond.convergence <- vector( "list", 
        length( normalize.condition ) )
    names( normalize.within.cond.convergence ) <- normalize.condition
    
    condition.sample.idx <- split( 1 : length( expression.condition ), 
        expression.condition )
    
    for ( sample.idx in condition.sample.idx )
    {
        condition <- as.character( expression.condition[ sample.idx[ 1 ] ] )
        norm.sample.idx <- which( condition == normalize.sample.condition )
        
        within.cond.norm.result <- normalize.median.within.condition(
            expression.data[ , sample.idx ], normalize.probe.idx, condition, 
            norm.probability, verbose )
                
        normalize.expr.data[ , norm.sample.idx ] <- 
            within.cond.norm.result$data
        
        normalize.within.cond.offset[ norm.sample.idx ] <- 
            within.cond.norm.result$offset
    }
    
    # normalize between conditions
    
    normalize.between.cond.offset <- rep( 0, length( normalize.condition ) )
    names( normalize.between.cond.offset ) <- normalize.condition
    
    normalize.between.cond.h0.probe <- NULL
    normalize.between.cond.convergence <- NULL
    
    if ( length( normalize.condition ) > 1 )
    {
        between.cond.norm.result <- 
            normalize.median.between.condition( normalize.expr.data, 
                normalize.probe.idx, normalize.condition, 
                normalize.sample.condition, norm.probability, verbose, 
                p.value.graph )
        
        normalize.expr.data <- between.cond.norm.result$data
        
        normalize.between.cond.offset <- between.cond.norm.result$offset
        
        normalize.between.cond.h0.probe <- between.cond.norm.result$h0.probe
        
        normalize.between.cond.convergence <-
            between.cond.norm.result$convergence
    }
    
    normalize.offset <- normalize.within.cond.offset[ normalize.sample ] + 
        normalize.between.cond.offset[ normalize.sample.condition ]
    
    list( data = normalize.expr.data, offset = normalize.offset, 
        within.condition.offset = normalize.within.cond.offset, 
        between.condition.offset = normalize.between.cond.offset, 
        between.condition.h0.probe = normalize.between.cond.h0.probe, 
        between.condition.convergence = normalize.between.cond.convergence )
}


normalize.median.within.condition <- function( edata, norm.probe.idx, 
    condition, norm.prob, verbose )
{
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    # find median offset
    norm.median.offset <- apply( edata[ norm.probe.idx, ], 2, quantile, 
        probs=norm.prob )
    norm.median.offset <- norm.median.offset - mean( norm.median.offset )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.median.offset )
    
    list( data = edata, offset = norm.median.offset )
}


normalize.median.between.condition <- function( edata, norm.probe.idx, 
    norm.cond, norm.sample.cond, norm.prob, verbose, p.value.graph )
{
    # identify samples available per condition
    within.cond.n <- as.vector( table( norm.sample.cond )[ norm.cond ] )
    names( within.cond.n ) <- norm.cond
    
    # calculate balanced within-condition means
    bal.mean.n <- min( within.cond.n )
    
    expr.bal.mean.data <- sapply( norm.cond, function( cond ) 
        if ( within.cond.n[ cond ] == bal.mean.n )
            rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] )
        else  # within.cond.n[ cond ] > bal.mean.n
            rowMeans( t( apply( 
                t( edata[ norm.probe.idx, cond == norm.sample.cond ] ), 
                2, function( ed ) sample( ed, bal.mean.n ) ) ) )
    )
    
    # for F-statistic
    # calculate within-condition means
    expr.mean.data <- sapply( norm.cond, function( cond ) 
        rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] ) )
    
    # calculate within-condition variances
    within.cond.var <- sweep( sapply( norm.cond, function( cond ) 
        apply( edata[ norm.probe.idx, cond == norm.sample.cond ], 1, var ) ), 
        2, within.cond.n - 1, "*" )
    within.cond.var <- rowSums( within.cond.var ) / 
        ( sum( within.cond.n ) - length( norm.cond ) )
    
    # calculate normalization
    norm.median.result <- normalize.median.selection( expr.bal.mean.data, 
        expr.mean.data, within.cond.var, within.cond.n, norm.prob, verbose, 
        p.value.graph )
    
    # normalize each condition with its offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <-
            edata[ , cond == norm.sample.cond ] - 
            norm.median.result$offset[ cond ]
    
    list( data = edata, offset = norm.median.result$offset, 
        h0.probe = norm.median.result$h0.probe, 
        convergence = norm.median.result$convergence )
}


normalize.median.selection <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, verbose, p.value.graph )
{
    iter.max <- 100
    median.offset.accum.step.max <- 10
    median.offset.accum.threshold <- 0.1
    median.offset.single.threshold <- 0.01
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    last.median.offset <- vector( "numeric" )
    last.median.h0.probe <- vector( "character" )
    
    norm.median.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.median.offset.ratio <- vector( "numeric" )
    norm.median.offset.accum.step <- vector( "numeric" )
    norm.median.h0.probe.num <- vector( "numeric" )

    median.offset <- rep( 0, ncol( edata ) )
    
    if ( verbose )
        cat( paste0( "between.condition\n" ) )
    
    iter <- 0
    median.offset.accum.step <- 0
    median.offset.ratio <- 1
    
    while ( iter < iter.max && 
        median.offset.accum.step < median.offset.accum.step.max &&
        median.offset.ratio >= median.offset.single.threshold )
    {
        iter <- iter + 1
        
        # obtain next step of median offset     
        median.offset.step <- calculate.median.offset( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, norm.prob, 
            p.value.graph, iter )
        
        median.offset.delta <- median.offset.step$value
        median.h0.probe <- median.offset.step$h0.probe
        
        median.h0.probe.num <- length( median.h0.probe )
        
        # check errors
        if ( any( is.nan( median.offset.delta ) ) )
            stop( "NaN error in normalize.median.condec" )
        
        # update total median offset
        median.offset <- median.offset + median.offset.delta
        median.offset <- median.offset - mean( median.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, median.offset )
        if ( ! is.null( edata.fstat ) )
            last.norm.data.fstat <- sweep( edata.fstat, 2, median.offset )
        
        # check convergence
        median.offset.sd <- sd( median.offset )
        median.offset.delta.sd <- sd( median.offset.delta )
        
        median.offset.ratio <- median.offset.delta.sd / median.offset.sd
        
        if ( median.offset.ratio < median.offset.accum.threshold )
            median.offset.accum.step <- median.offset.accum.step + 1
        else
            median.offset.accum.step <- 0
        
        # store last results
        last.median.offset <- median.offset
        last.median.h0.probe <- median.h0.probe
        
        # store step results
        norm.median.offset <- rbind( norm.median.offset, median.offset )
        norm.median.offset.ratio <- c( norm.median.offset.ratio, 
            median.offset.ratio )
        norm.median.offset.accum.step <- c( norm.median.offset.accum.step, 
            median.offset.accum.step )
        norm.median.h0.probe.num <- c( norm.median.h0.probe.num, 
            median.h0.probe.num )
        
        if ( verbose )
            cat( sprintf( "  %2d %g %g %02d %d\n", iter, median.offset.sd, 
                median.offset.ratio, median.offset.accum.step, 
                median.h0.probe.num ) )
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( median.offset.accum.step < median.offset.accum.step.max && 
            median.offset.ratio >= median.offset.single.threshold )
        stop( "no convergence in normalize.median.condec" )
    
    dimnames( norm.median.offset ) <- NULL
    
    norm.median.convergence <- list( offset = norm.median.offset, 
        offset.ratio = norm.median.offset.ratio, 
        offset.accum.step = norm.median.offset.accum.step, 
        h0.probe.num = norm.median.h0.probe.num )
    
    list( offset = last.median.offset, h0.probe = last.median.h0.probe, 
        convergence = norm.median.convergence )
}


calculate.median.offset <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, norm.prob, p.value.graph, iter )
{
    ks.test.alpha <- 1e-3
    
    # calculate f statistics for each probe
    expr.k <- length( within.cond.n )
    expr.n <- sum( within.cond.n )
    
    expr.grand.mean <- apply( edata.fstat, 1, function( ef ) 
        sum( ef * within.cond.n ) ) / expr.n
    
    between.cond.var <- apply( ( edata.fstat - expr.grand.mean )^2, 1, 
        function( ef2 ) sum( ef2 * within.cond.n ) ) / ( expr.k - 1 )
    
    expr.f <- between.cond.var / within.cond.var
    
    expr.f <- na.omit( expr.f )  # in case of 0/0
    attr( expr.f, "na.action" ) <- NULL
    
    expr.p.value <- pf( expr.f, df1 = expr.k-1, df2 = expr.n-expr.k, 
        lower.tail = FALSE )
    
    # identify H0 probes with one-sided up Kolmogorov-Smirnov test
    ks.test.d <- sqrt( - log( ks.test.alpha ) / 2 )
    
    epv <- sort( expr.p.value )
    epv.n <- length( expr.p.value )
    
    epv.i <- 1
    ks.test.D.up <- 1:epv.n / epv.n - epv
    ks.test.reject <- any( ks.test.D.up > ks.test.d / sqrt( epv.n ) )
    
    while ( ks.test.reject )
    {
        epv.i <- epv.i + 1
        
        ks.test.D.up <- 
            ( epv.i : epv.n - epv.i + 1 ) / ( epv.n - epv.i + 1 ) - 
            ( epv[ epv.i : epv.n ] - epv[ epv.i - 1 ] ) / 
            ( 1 - epv[ epv.i - 1 ] )
        
        ks.test.reject <- 
            any( ks.test.D.up > ks.test.d / sqrt( epv.n - epv.i + 1 ) )
    }
    
    epv.h0.i <- epv.i
    epv.h0.n <- epv.n - epv.i + 1
    epv.h0.p <-  ifelse( epv.i == 1, 0, epv[ epv.i - 1 ] )
    epv.h0.q <- ( epv.i - 1 ) / epv.n
    
    h0.probe.idx <- order( expr.p.value, decreasing=TRUE )[ 1 : epv.h0.n ]
    h0.probe <- names( expr.p.value )[ h0.probe.idx ]
    
    pi0.est <- ( 1 - epv.h0.q ) / ( 1 - epv.h0.p )
    
    # plot graph of p-values
    if ( ! is.null( p.value.graph ) )
    {
        if ( p.value.graph != "" )
        {
            if ( ! file.exists( p.value.graph ) )
                dir.create( p.value.graph, recursive=TRUE )
            
            png.filename <- sprintf( "%s/pvalue_iter%02d.png", 
                p.value.graph, iter )
            
            png( png.filename, width=640, height=360 )
        }
        
        par.default <- par( no.readonly=TRUE )
        
        par( mfrow = c(1,2), pty="s", xpd=FALSE, 
            mar = c( 3.4, 2.2, 2.2, 2.2 ), oma = c( 0, 3.8, 0, 0 ) )
        
        epv.quant <- 1:epv.n / epv.n
        epv.h0.idx <- epv.h0.i : epv.n
        
        epv.x0 <- epv.h0.p
        epv.y0 <- epv.h0.q
        epv.yd <- ks.test.d * ( sqrt( epv.h0.n ) / epv.n )
        
        xylim <- list( c(0,0), c( epv.h0.p, epv.h0.q ) )
        
        for ( i in 1:2  )
        {
            plot( 0, type="n", 
                xlim = c( xylim[[i]][1], 1 ), ylim = c( xylim[[i]][2], 1 ), 
                xlab="p-value", ylab="", cex.axis=1.3, cex.lab=1.5 )
            
            segments( x0 = c( 0, epv[ - epv.h0.idx ] ), 
                y0 = c( 0, epv.quant[ - epv.h0.idx ] ), 
                x1 = c( epv[ - epv.h0.idx ], epv[ epv.h0.idx ][ 1 ] ), 
                lwd=2 )
            segments( x0 = epv[ epv.h0.idx ], y0 = epv.quant[ epv.h0.idx ], 
                x1 = c( epv[ epv.h0.idx ][ -1 ], 1 ), lwd=2, col="red" )
            points( epv.h0.p, epv.h0.q, pch=20, cex=2.5 )
            
            segments( 0, 1 - pi0.est, 1, 1, col="blue", lwd=1.5, lty=2 )
            segments( 0, 1 - pi0.est + epv.yd, 1, 1 + epv.yd, col="blue", 
                lwd=1.5, lty=3 )
            
            if ( i==1 )
                graph.title <- sprintf( "iter=%02d", iter )
            else
                graph.title <- substitute( paste( "#", H[0], "=", h0.n ), 
                    list( h0.n = sprintf( "%5d", epv.h0.n ) ) )
            
            title( main=graph.title, line = ifelse( i==1, 1.4, 1.75 ), 
                font.main=1, cex.main=1.7 )
        }
        
        mtext( "F( p-value )", side=2, line=1.3, outer=TRUE, cex=1.5 )
        
        if ( p.value.graph != "" )
            dev.off()
        else
            par( par.default )
    }
    
    # identify probes for normalization
    median.probe <- h0.probe
    
    # calculate offset
    median.offset <- apply( edata[ median.probe, ], 2, quantile, 
        probs=norm.prob )
    median.offset <- median.offset - mean( median.offset )
    
    list( value = median.offset, h0.probe = h0.probe )
}

# Implements standard-vector condition-decompositon normalization


normalize.stdvec.condec <- function( expression.data, expression.condition, 
    verbose=FALSE, vector.graph=NULL, p.value.graph=NULL )
{
    # identify samples and conditions to normalize
    normalize.sample <- 
        colnames( expression.data )[ ! is.na( expression.condition ) ]
    normalize.sample.condition <- 
        as.character( na.omit( expression.condition ) )
    normalize.condition <- unique( normalize.sample.condition )
    
    if ( length( normalize.condition ) == 0 )
        stop( "No condition to normalize" )
    else if ( min( table( normalize.sample.condition ) ) < 2 )
        stop( "There must be 2 or more samples in each condition" )
    
    # select probes for normalization
    # no missing values in any normalization sample
    expression.probe <- rownames( expression.data )
    normalize.probe.idx <- which( rowSums( 
        is.na( expression.data[ , normalize.sample ] ) ) == 0 )
    
    # normalize within conditions
    
    normalize.expr.data <- matrix( nrow = length( expression.probe ), 
        ncol = length( normalize.sample ) )
    rownames( normalize.expr.data ) <- expression.probe
    colnames( normalize.expr.data ) <- normalize.sample
    
    normalize.within.cond.offset <- rep( 0, length( normalize.sample ) )
    names( normalize.within.cond.offset ) <- normalize.sample
    
    normalize.within.cond.convergence <- vector( "list", 
        length( normalize.condition ) )
    names( normalize.within.cond.convergence ) <- normalize.condition
    
    condition.sample.idx <- split( 1 : length( expression.condition ), 
        expression.condition )
    
    for ( sample.idx in condition.sample.idx )
    {
        condition <- as.character( expression.condition[ sample.idx[ 1 ] ] )
        norm.sample.idx <- which( condition == normalize.sample.condition )
        
        within.cond.norm.result <- normalize.stdvec.within.condition( 
            expression.data[ , sample.idx ], normalize.probe.idx, condition, 
            verbose, vector.graph )
        
        normalize.expr.data[ , norm.sample.idx ] <- 
            within.cond.norm.result$data
        
        normalize.within.cond.offset[ norm.sample.idx ] <- 
            within.cond.norm.result$offset
        
        normalize.within.cond.convergence[[ condition ]] <- 
            within.cond.norm.result$convergence
    }
    
    # normalize between conditions
    
    normalize.between.cond.offset <- rep( 0, length( normalize.condition ) )
    names( normalize.between.cond.offset ) <- normalize.condition
    
    normalize.between.cond.h0.probe <- NULL
    normalize.between.cond.convergence <- NULL
    
    if ( length( normalize.condition ) > 1 )
    {
        between.cond.norm.result <- 
            normalize.stdvec.between.condition( normalize.expr.data, 
                normalize.probe.idx, normalize.condition, 
                normalize.sample.condition, verbose, vector.graph, 
                p.value.graph )
        
        normalize.expr.data <- between.cond.norm.result$data
        
        normalize.between.cond.offset <- between.cond.norm.result$offset
        
        normalize.between.cond.h0.probe <- between.cond.norm.result$h0.probe
        
        normalize.between.cond.convergence <-
            between.cond.norm.result$convergence
    }
    
    normalize.offset <- normalize.within.cond.offset[ normalize.sample ] + 
        normalize.between.cond.offset[ normalize.sample.condition ]
    
    list( data = normalize.expr.data, offset = normalize.offset, 
        within.condition.offset = normalize.within.cond.offset, 
        between.condition.offset = normalize.between.cond.offset, 
        between.condition.h0.probe = normalize.between.cond.h0.probe, 
        within.condition.convergence = normalize.within.cond.convergence, 
        between.condition.convergence = normalize.between.cond.convergence )
}


normalize.stdvec.within.condition <- function( edata, norm.probe.idx, 
    condition, verbose, vector.graph )
{
    # calculate normalization
    norm.stdvec.result <- normalize.standard.vector( edata[ norm.probe.idx, ], 
        NULL, NULL, NULL, condition, verbose, vector.graph, NULL )
    
    # normalize with obtained offset
    edata <- sweep( edata, 2, norm.stdvec.result$offset )
    
    list( data = edata, offset = norm.stdvec.result$offset, 
        convergence = norm.stdvec.result$convergence )
}


normalize.stdvec.between.condition <- function( edata, norm.probe.idx, 
    norm.cond, norm.sample.cond, verbose, vector.graph, p.value.graph )
{
    # identify samples available per condition
    within.cond.n <- as.vector( table( norm.sample.cond )[ norm.cond ] )
    names( within.cond.n ) <- norm.cond
    
    # calculate balanced within-condition means
    bal.mean.n <- min( within.cond.n )
    
    expr.bal.mean.data <- sapply( norm.cond, function( cond ) 
        if ( within.cond.n[ cond ] == bal.mean.n )
            rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] )
        else  # within.cond.n[ cond ] > bal.mean.n
            rowMeans( t( apply( 
                t( edata[ norm.probe.idx, cond == norm.sample.cond ] ), 
                2, function( ed ) sample( ed, bal.mean.n ) ) ) )
    )
    
    # for F-statistic
    # calculate within-condition means
    expr.mean.data <- sapply( norm.cond, function( cond ) 
        rowMeans( edata[ norm.probe.idx, cond == norm.sample.cond ] ) )
    
    # calculate within-condition variances
    within.cond.var <- sweep( sapply( norm.cond, function( cond ) 
        apply( edata[ norm.probe.idx, cond == norm.sample.cond ], 1, var ) ), 
        2, within.cond.n - 1, "*" )
    within.cond.var <- rowSums( within.cond.var ) / 
        ( sum( within.cond.n ) - length( norm.cond ) )
    
    # calculate normalization
    norm.stdvec.result <- normalize.standard.vector( expr.bal.mean.data, 
        expr.mean.data, within.cond.var, within.cond.n, "between.condition", 
        verbose, vector.graph, p.value.graph )
    
    # normalize each condition with its offset
    for( cond in norm.cond )
        edata[ , cond == norm.sample.cond ] <-
            edata[ , cond == norm.sample.cond ] - 
            norm.stdvec.result$offset[ cond ]
    
    list( data = edata, offset = norm.stdvec.result$offset, 
        h0.probe = norm.stdvec.result$h0.probe, 
        convergence = norm.stdvec.result$convergence )
}


normalize.standard.vector <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, condition, verbose, vector.graph, p.value.graph )
{
    iter.max <- 100
    stdvec.offset.accum.step.max <- 10
    stdvec.offset.accum.threshold <- 0.1
    stdvec.offset.single.threshold <- 0.01
    vector.graph.probe.num <- 5000
    
    last.norm.data <- edata
    last.norm.data.fstat <- edata.fstat
    
    last.stdvec.offset <- vector( "numeric" )
    last.stdvec.h0.probe <- vector( "character" )
    
    norm.stdvec.offset <- matrix( nrow=0, ncol = ncol( edata ) )
    norm.stdvec.offset.sd <- vector( "numeric" )
    norm.stdvec.offset.stderr <- vector( "numeric" )
    norm.stdvec.offset.delta.sd <- vector( "numeric" )
    norm.stdvec.offset.accum.step <- vector( "numeric" )
    norm.stdvec.numerical.demand <- vector( "numeric" )
    norm.stdvec.watson.u2 <- matrix( nrow=0, 
        ncol = ( ncol( edata ) - 1 ) %/% 3 + 1 )
    norm.stdvec.h0.probe.num <- vector( "numeric" )
    
    stdvec.offset <- rep( 0, ncol( edata ) )
    
    vector.graph.probe <- NULL
    
    if ( ! is.null( vector.graph ) )
    {
        # select a sample of probes for plotting standardized sample vectors
        edata.probe <- rownames( edata )
        edata.probe.num <- length( edata.probe )
        
        if ( edata.probe.num <= vector.graph.probe.num )
            vector.graph.probe <- edata.probe
        else
            vector.graph.probe <- sample( edata.probe, vector.graph.probe.num )
    }
    
    if ( verbose )
        cat( paste0( condition, "\n" ) )
    
    iter <- 0
    stdvec.offset.accum.step <- 0
    stdvec.offset.ratio <- 1
    
    while ( iter < iter.max && 
        stdvec.offset.accum.step < stdvec.offset.accum.step.max &&
        stdvec.offset.ratio >= stdvec.offset.single.threshold )
    {
        iter <- iter + 1
        
        # obtain next step of standard vector offset     
        stdvec.offset.step <- calculate.stdvec.offset( last.norm.data, 
            last.norm.data.fstat, within.cond.var, within.cond.n, vector.graph, 
            vector.graph.probe, p.value.graph, condition, iter )
        
        stdvec.offset.delta <- stdvec.offset.step$value
        stdvec.offset.stderr <- stdvec.offset.step$stderr
        stdvec.numerical.demand <- stdvec.offset.step$numerical.demand
        stdvec.watson.u2 <- stdvec.offset.step$watson.u2
        stdvec.h0.probe <- stdvec.offset.step$h0.probe
        
        stdvec.h0.probe.num <- NULL
        if ( ! is.null( stdvec.h0.probe ) )
            stdvec.h0.probe.num <- length( stdvec.h0.probe )
        
        # check errors
        if ( any( is.nan( stdvec.offset.delta ) ) || 
                is.nan( stdvec.offset.stderr ) )
            stop( "NaN error in normalize.stdvec.condec" )
        
        if ( stdvec.numerical.demand < .Machine$double.eps * 10^3 )
            stop( "numerical error in normalize.stdvec.condec" )
        
        # update total standard vector offset
        stdvec.offset <- stdvec.offset + stdvec.offset.delta
        stdvec.offset <- stdvec.offset - mean( stdvec.offset )
        
        # update data at once with total offset
        last.norm.data <- sweep( edata, 2, stdvec.offset )
        if ( ! is.null( edata.fstat ) )
            last.norm.data.fstat <- sweep( edata.fstat, 2, stdvec.offset )
        
        # check convergence
        stdvec.offset.sd <- sd( stdvec.offset )
        stdvec.offset.delta.sd <- sd( stdvec.offset.delta )
        
        stdvec.offset.stderr.ratio <- 
            stdvec.offset.stderr / stdvec.offset.sd
        stdvec.offset.delta.sd.ratio <- 
            stdvec.offset.delta.sd / stdvec.offset.sd
        
        stdvec.offset.ratio <- stdvec.offset.delta.sd / stdvec.offset.stderr
        
        if ( stdvec.offset.ratio < stdvec.offset.accum.threshold )
            stdvec.offset.accum.step <- stdvec.offset.accum.step + 1
        else
            stdvec.offset.accum.step <- 0
        
        # store last results
        last.stdvec.offset <- stdvec.offset
        last.stdvec.h0.probe <- stdvec.h0.probe
        
        # store step results
        norm.stdvec.offset <- rbind( norm.stdvec.offset, stdvec.offset )
        norm.stdvec.offset.sd <- c( norm.stdvec.offset.sd, stdvec.offset.sd )
        norm.stdvec.offset.stderr <- c( norm.stdvec.offset.stderr, 
            stdvec.offset.stderr )
        norm.stdvec.offset.delta.sd <- c( norm.stdvec.offset.delta.sd, 
            stdvec.offset.delta.sd )
        norm.stdvec.offset.accum.step <- c( norm.stdvec.offset.accum.step, 
            stdvec.offset.accum.step )
        norm.stdvec.numerical.demand <- c( norm.stdvec.numerical.demand, 
            stdvec.numerical.demand )
        if ( ! is.null( stdvec.watson.u2 ) )
            norm.stdvec.watson.u2 <- rbind( norm.stdvec.watson.u2, 
                stdvec.watson.u2 )
        if ( ! is.null( stdvec.h0.probe.num ) )
            norm.stdvec.h0.probe.num <- c( norm.stdvec.h0.probe.num, 
                stdvec.h0.probe.num )
        
        if ( verbose )
        {
            if ( ! is.null( stdvec.watson.u2 ) )
                stdvec.watson.u2.char <- paste0( signif( stdvec.watson.u2, 6 ), 
                    collapse=" " )
            else
                stdvec.watson.u2.char <- ""
            
            cat( sprintf( "  %2d %g %g %g %02d %g [%s]", iter, 
                stdvec.offset.sd, stdvec.offset.stderr.ratio, 
                stdvec.offset.delta.sd.ratio, stdvec.offset.accum.step, 
                stdvec.numerical.demand, stdvec.watson.u2.char ), 
                stdvec.h0.probe.num, "\n" )
        }
    }
    
    if ( verbose )
        cat( "\n" )
    
    if ( stdvec.offset.accum.step < stdvec.offset.accum.step.max && 
        stdvec.offset.ratio >= stdvec.offset.single.threshold )
        stop( "no convergence in normalize.stdvec.condec" )
    
    dimnames( norm.stdvec.offset ) <- NULL
    dimnames( norm.stdvec.watson.u2 ) <- NULL
    
    norm.stdvec.convergence <- list( offset = norm.stdvec.offset, 
        offset.sd = norm.stdvec.offset.sd, 
        offset.stderr = norm.stdvec.offset.stderr, 
        offset.delta.sd = norm.stdvec.offset.delta.sd, 
        offset.accum.step = norm.stdvec.offset.accum.step, 
        numerical.demand = norm.stdvec.numerical.demand, 
        watson.u2 = norm.stdvec.watson.u2 )
    
    if ( length( norm.stdvec.h0.probe.num ) > 0 )
        norm.stdvec.convergence$h0.probe.num <- norm.stdvec.h0.probe.num
    
    if ( is.null( last.stdvec.h0.probe ) ) {
        list( offset = last.stdvec.offset, 
            convergence = norm.stdvec.convergence )
    } else {
        list( offset = last.stdvec.offset, h0.probe = last.stdvec.h0.probe, 
            convergence = norm.stdvec.convergence )
    }
}


calculate.stdvec.offset <- function( edata, edata.fstat, within.cond.var, 
    within.cond.n, vector.graph, vector.graph.probe, p.value.graph, condition, 
    iter )
{
    stdvec.trim <- 0.01
    ks.test.alpha <- 1e-3
    
    h0.probe <- NULL
    
    if ( is.null( edata.fstat ) )
    {
        # identify probes for normalization
        expr.var <- apply( edata, 1, var )
        expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
        
        stdvec.probe <- names( which( 
            expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
            expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ) )
    }
    else
    {
        # calculate f statistics for each probe
        expr.k <- length( within.cond.n )
        expr.n <- sum( within.cond.n )
        
        expr.grand.mean <- apply( edata.fstat, 1, function( ef ) 
            sum( ef * within.cond.n ) ) / expr.n
        
        between.cond.var <- apply( ( edata.fstat - expr.grand.mean )^2, 1, 
            function( ef2 ) sum( ef2 * within.cond.n ) ) / ( expr.k - 1 )
        
        expr.f <- between.cond.var / within.cond.var
        
        expr.f <- na.omit( expr.f )  # in case of 0/0
        attr( expr.f, "na.action" ) <- NULL
        
        expr.p.value <- pf( expr.f, df1 = expr.k-1, df2 = expr.n-expr.k, 
            lower.tail = FALSE )
        
        # identify H0 probes with one-sided up Kolmogorov-Smirnov test
        ks.test.d <- sqrt( - log( ks.test.alpha ) / 2 )
        
        epv <- sort( expr.p.value )
        epv.n <- length( expr.p.value )
        
        epv.i <- 1
        ks.test.D.up <- 1:epv.n / epv.n - epv
        ks.test.reject <- any( ks.test.D.up > ks.test.d / sqrt( epv.n ) )
        
        while ( ks.test.reject )
        {
            epv.i <- epv.i + 1
            
            ks.test.D.up <- 
                ( epv.i : epv.n - epv.i + 1 ) / ( epv.n - epv.i + 1 ) - 
                ( epv[ epv.i : epv.n ] - epv[ epv.i - 1 ] ) / 
                    ( 1 - epv[ epv.i - 1 ] )
            
            ks.test.reject <- 
                any( ks.test.D.up > ks.test.d / sqrt( epv.n - epv.i + 1 ) )
        }
        
        epv.h0.i <- epv.i
        epv.h0.n <- epv.n - epv.i + 1
        epv.h0.p <-  ifelse( epv.i == 1, 0, epv[ epv.i - 1 ] )
        epv.h0.q <- ( epv.i - 1 ) / epv.n
        
        h0.probe.idx <- order( expr.p.value, decreasing=TRUE )[ 1 : epv.h0.n ]
        h0.probe <- names( expr.p.value )[ h0.probe.idx ]
        
        pi0.est <- ( 1 - epv.h0.q ) / ( 1 - epv.h0.p )
                
        # plot graph of p-values
        if ( ! is.null( p.value.graph ) )
        {
            if ( p.value.graph != "" )
            {
                if ( ! file.exists( p.value.graph ) )
                    dir.create( p.value.graph, recursive=TRUE )
                
                png.filename <- sprintf( "%s/pvalue_iter%02d.png", 
                    p.value.graph, iter )
                
                png( png.filename, width=640, height=360 )
            }
            
            par.default <- par( no.readonly=TRUE )
            
            par( mfrow = c(1,2), pty="s", xpd=FALSE, 
                mar = c( 3.2, 2.2, 1.6, 2.2 ), oma = c( 0, 3.8, 0, 0 ) )
            
            epv.quant <- 1:epv.n / epv.n
            epv.h0.idx <- epv.h0.i : epv.n
            
            epv.x0 <- epv.h0.p
            epv.y0 <- epv.h0.q
            epv.yd <- ks.test.d * ( sqrt( epv.h0.n ) / epv.n )
            
            xylim <- list( c(0,0), c( epv.h0.p, epv.h0.q ) )
            
            for ( i in 1:2  )
            {
                plot( 0, type="n", 
                    xlim = c( xylim[[i]][1], 1 ), ylim = c( xylim[[i]][2], 1 ), 
                    xlab="p-value", ylab="", cex.axis=1.3, cex.lab=1.5 )
                
                segments( x0 = c( 0, epv[ - epv.h0.idx ] ), 
                    y0 = c( 0, epv.quant[ - epv.h0.idx ] ), 
                    x1 = c( epv[ - epv.h0.idx ], epv[ epv.h0.idx ][ 1 ] ), 
                    lwd=2 )
                segments( x0 = epv[ epv.h0.idx ], y0 = epv.quant[ epv.h0.idx ], 
                    x1 = c( epv[ epv.h0.idx ][ -1 ], 1 ), lwd=2, col="red" )
                points( epv.h0.p, epv.h0.q, pch=20, cex=2.5 )
                
                segments( 0, 1 - pi0.est, 1, 1, col="blue", lwd=1.5, lty=2 )
                segments( 0, 1 - pi0.est + epv.yd, 1, 1 + epv.yd, col="blue", 
                    lwd=1.5, lty=3 )
                
                if ( i==1 )
                    graph.title <- sprintf( "iter=%02d", iter )
                else
                    graph.title <- substitute( paste( "#", H[0], "=", h0.n ), 
                        list( h0.n = sprintf( "%5d", epv.h0.n ) ) )
                
                title( main=graph.title, line = ifelse( i==1, 1.4, 1.75 ), 
                    font.main=1, cex.main=1.7 )
            }
            
            mtext( "F( p-value )", side=2, line=1.3, outer=TRUE, cex=1.5 )
            
            if ( p.value.graph != "" )
                dev.off()
            else
                par( par.default )
        }
        
        # identify probes for normalization
        expr.var <- apply( edata.fstat[ h0.probe, ], 1, var )
        expr.var[ expr.var < max( expr.var ) * .Machine$double.eps  ] <- NA
        
        stdvec.probe <- h0.probe[ 
            expr.var > quantile( expr.var, stdvec.trim/2, na.rm=TRUE ) & 
            expr.var < quantile( expr.var, 1 - stdvec.trim/2, na.rm=TRUE ) ]
    }
    
    print(paste(length(stdvec.probe), nrow(edata)))
    print(all(stdvec.probe %in% rownames(edata)))
    if (!all(stdvec.probe %in% rownames(edata))) {
        print(head(stdvec.probe))
        print(head(rownames(edata)))
        print(setdiff(stdvec.probe, rownames(edata)))
    }
    
    # center and scale expression data
    expr.mean <- rowMeans( edata[ stdvec.probe, ] )
    expr.centered <- sweep( edata[ stdvec.probe, ], 1, expr.mean, "-" )
    
    expr.sd.inv <- 1 / apply( expr.centered, 1, sd )
    expr.scaled <- sweep( expr.centered, 1, expr.sd.inv, "*" )
    
    # calculate offset
    expr.sd.inv.sum <- sum( expr.sd.inv )
    stdvec.offset <- apply( expr.scaled, 2, sum ) / expr.sd.inv.sum
    stdvec.offset <- stdvec.offset - mean( stdvec.offset )
    
    # estimate error and numerical demand
    stdvec.offset.stderr <- sqrt( length( stdvec.probe ) ) / expr.sd.inv.sum
    expr.sd.inv.min <- min( expr.sd.inv )
    stdvec.numerical.demand <- expr.sd.inv.min / expr.sd.inv.sum
    
    # calculate and plot density distribution of standard vector angles
    
    dimension.num <- ncol( expr.scaled )
    
    if ( dimension.num < 3 )
    {
        theta.watson.u2 = NULL
    }
    else
    {
        # identify condition groups
        dimension.group.num <- ( dimension.num - 1 ) %/% 3 + 1
        
        dimension.group <- lapply( 1 : dimension.group.num, function ( g )
            if ( g < dimension.group.num )
                ( 3*g - 2 ) : ( 3*g )
            else
                ( dimension.num - 2 ) : dimension.num )
        
        theta.watson.u2 <- vector( "numeric", dimension.group.num )
        
        uv <- matrix( c( 0, -1/sqrt(2), 1/sqrt(2), 
            2/sqrt(6), -1/sqrt(6), -1/sqrt(6) ), nrow=3 )
        
        for ( dim.group.idx in 1:dimension.group.num )
        {
            # select expression values for each condition group
            expr.dim.group <- expr.scaled[ , 
                dimension.group[[ dim.group.idx ]] ]
            
            if ( dimension.group.num > 1 )
            {
                # re-standardize again for this group
                expr.dim.group.mean <- rowMeans( expr.dim.group )
                expr.dim.group.centered <- sweep( expr.dim.group, 1, 
                    expr.dim.group.mean, "-" )
                
                expr.dim.group.sd <- apply( expr.dim.group.centered, 1, sd )
                expr.dim.group.sel <- expr.dim.group.sd != 0
                expr.dim.group <- sweep( 
                    expr.dim.group.centered[ expr.dim.group.sel, ], 1, 
                    expr.dim.group.sd[ expr.dim.group.sel ], "/" )
            }
            
            expr.uv <- expr.dim.group %*% uv
            expr.u <- expr.uv[ , 1 ]
            expr.v <- expr.uv[ , 2 ]
            
            # calculate density distribution of angles
            theta.density.n <- 2^11
            theta.density.adjust <- 0.5
            
            expr.theta <- atan2( expr.v, expr.u )
            expr.theta.ex <- c( expr.theta, 
                expr.theta + ifelse( expr.theta > 0, -2*pi, 2*pi ) )
            
            expr.theta.density <- density( expr.theta.ex, 
                adjust=theta.density.adjust, n=theta.density.n )
            expr.theta.density.sel <- expr.theta.density$x > -pi & 
                expr.theta.density$x <= pi
            
            # calculate density distribution of angles after permutations
            expr.theta.permu <- cbind( expr.theta.ex, 
                expr.theta.ex + (2*pi)/3, 
                expr.theta.ex - (2*pi)/3, 
                - expr.theta.ex + pi, 
                - expr.theta.ex + pi/3, 
                - expr.theta.ex - pi/3 )
            
            invar.theta <- expr.theta.permu[ expr.theta.permu > -pi &
                expr.theta.permu <= pi ]
            invar.theta.ex <- c( invar.theta, 
                invar.theta + ifelse( invar.theta > 0, -2*pi, 2*pi ) )
            
            invar.theta.density <- density( invar.theta.ex, 
                adjust=theta.density.adjust, n=theta.density.n )
            invar.theta.density.sel <- invar.theta.density$x > -pi & 
                invar.theta.density$x <= pi
            
            # calculate Watson U2 statistic
            theta.watson.u2[ dim.group.idx ] <- 
                watson.u2( expr.theta, sample( invar.theta, 
                    length( expr.theta ), replace=TRUE ) )
            
            if ( ! is.null( vector.graph ) && require( plotrix, quietly=TRUE ) )
            {
                if ( vector.graph != "" )
                {
                    if ( ! file.exists( vector.graph ) )
                        dir.create( vector.graph, recursive=TRUE )
                    
                    png.filename <- sprintf( 
                        "%s/stdvec_%s%s_iter%02d.png", vector.graph, condition, 
                        ifelse( dimension.group.num > 1, 
                            sprintf( "_dg%02d", dim.group.idx ), "" ), 
                        iter )
                    
                    png( png.filename, width=640, height=360 )
                }
                
                par.default <- par( no.readonly=TRUE )
                
                par( mfrow = c(1,2), pty="s", xpd=FALSE, 
                    mar = c( 1.2, 1.2, 2.2, 1.2 ), oma = c( 0, 0, 0, 0 ) )
                
                # select offset values for condition group
                stdvec.offset.uv <- stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ] %*% uv
                stdvec.offset.u <- stdvec.offset.uv[ , 1 ]
                stdvec.offset.v <- stdvec.offset.uv[ , 2 ]

                # plot a sample of standardized sample vectors
                uv.lim <- c( -1.5, 1.5 )
                plot( 0, type="n", xlim=uv.lim, ylim=uv.lim, axes=FALSE, 
                    ann=FALSE, frame.plot=FALSE, asp=1 )
                
                expr.uv.probe <- names( expr.u )
                if ( length( expr.uv.probe ) > length( vector.graph.probe ) ) {
                    expr.uv.probe <- intersect( expr.uv.probe, 
                        vector.graph.probe )
                }
                expr.uv.probe.num <- length( expr.uv.probe )
                
                expr.uv.factor <- 1.05
                expr.uv.color <- gray( 0.3 )
                expr.uv.width <- ifelse( expr.uv.probe.num > 1000, 0.1, 0.2 )
                segments( 0, 0, expr.uv.factor * expr.u[ expr.uv.probe ], 
                    expr.uv.factor * expr.v[ expr.uv.probe ], 
                    col=expr.uv.color, lwd=expr.uv.width )
                
                grid.pos.x <- c( 0, -sqrt(3/4), sqrt(3/4) )
                grid.pos.y <- c( 1, -1/2, -1/2 )
                
                grid.length <- expr.uv.factor * sqrt( 2 )
                segments( 0, 0, grid.length * grid.pos.x, 
                    grid.length * grid.pos.y, lwd=2 )
                
                stdvec.offset.uv.factor <- 10
                stdvec.offset.uv.color <- "red"
                segments( 0, 0, stdvec.offset.uv.factor * stdvec.offset.u, 
                    stdvec.offset.uv.factor * stdvec.offset.v, lwd=2, 
                    col=stdvec.offset.uv.color )
                
                par( xpd=TRUE )
                
                grid.label.length <- 1.63
                grid.labels <- c( "s1", "s2", "s3" )
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels, cex=1.5 )
                
                par( xpd=FALSE )
                
                stdvec.offset.mag <- sqrt( sum( stdvec.offset[ 
                    dimension.group[[ dim.group.idx ]] ]^2 ) )
                mtext( paste0( "||offset|| = ", sprintf( "%.3e", 
                    stdvec.offset.mag ) ), side=1, line=0.2, 
                    cex=1.5 )
                
                # plot polar distributions of standard vector angles
                polar.expr.theta <- expr.theta.density$x[ 
                    expr.theta.density.sel ]
                polar.expr.rho <- 2 * expr.theta.density$y[ 
                    expr.theta.density.sel ]

                polar.invar.theta <- invar.theta.density$x[ 
                    invar.theta.density.sel ]
                polar.invar.rho <- 2 * invar.theta.density$y[ 
                    invar.theta.density.sel ]
                
                rho.grid <- seq( 0, 3/(2*pi), length.out=4 )
                rho.labels <- c( "", "", expression(1/pi), "" )
                theta.labels <- c( "", "", "" )
                
                radial.plot( polar.expr.rho, polar.expr.theta - pi/2, 
                    start=pi/2, rp.type="p", radial.lim = rho.grid, 
                    radial.labels = rho.labels, 
                    show.grid.labels = length( theta.labels ), 
                    labels = theta.labels, mar = par( "mar" ) )
                
                radial.plot( polar.invar.rho, polar.invar.theta - pi/2, 
                    start=pi/2, rp.type="p", radial.lim = rho.grid, lty=2, 
                    line.col="blue", add=TRUE )
                
                grid.label.length <- 0.523
                text( grid.label.length * grid.pos.x, 
                    grid.label.length * grid.pos.y, grid.labels, cex=1.5 )
                
                mtext( substitute( paste( "Watson U"^"2", " = ", wu2 ), 
                    list( wu2 = sprintf( "%.3e", 
                        theta.watson.u2[ dim.group.idx ] ) ) ), 
                    side=1, line=0.2, cex=1.5 )
                
                graph.title <- sprintf( "%s%s    -    iter=%02d", condition, 
                    ifelse( dimension.group.num > 1, 
                        sprintf( ":%02d", dim.group.idx ), "" ), 
                    iter )
                title( main=graph.title, outer=TRUE, line=-1.8, font.main=1, 
                    cex.main=1.7 )
                
                if ( vector.graph != "" )
                    dev.off()
                else
                    par( par.default )
            }
        }
    }
    
    list( value = stdvec.offset, stderr = stdvec.offset.stderr, 
        numerical.demand = stdvec.numerical.demand, 
        watson.u2 = theta.watson.u2, 
        h0.probe = h0.probe )
}


watson.u2 <- function( x, y )
{
    # see section 6.5 of Durbin, Distribution Theory for Tests Based on the 
    # Sample Distribution Function, SIAM, Philadelphia (1973)
    
    n <- length( x )
    m <- length( y )
    
    r <- c( sort( x ), sort( y ) )
    r.rank <- rank( r, ties.method="average" )
    
    z <- ( r.rank[ 1:n ] - 1:n ) / m - ( 1:n - 1/2 ) / n 
    
    ( m / (n+m) ) * sum( ( z - mean( z ) )^2 ) + 
        ( m*(m+2*n) ) / ( 12*n*m*(n+m) ) 
}