# This file contains functions written for the carcinosarcoma project to deal with sparse matrices.

# Number of nonzero elements per column:

col_nnz <- function(spmat, use_names = TRUE) {
    if(class(spmat) == 'dgTMatrix' || class(spmat) == 'dgRMatrix') {
        if(use_names) {
            setNames(tabulate(spmat@j + 1, nbins = spmat@Dim[2]), colnames(spmat))
        } else {
            tabulate(spmat@j + 1, nbins = spmat@Dim[2])
        }
    } else if(class(spmat) == 'dgCMatrix') {
        if(use_names) {
            setNames(diff(spmat@p), colnames(spmat))
        } else {
            diff(spmat@p)
        }
    } else {
        stop('Matrix should have class dgTMatrix, dgCMatrix or dgRMatrix.')
    }
}

# Number of nonzero elements per row:

row_nnz <- function(spmat, use_names = TRUE) {
    if(class(spmat) == 'dgTMatrix' || class(spmat) == 'dgCMatrix') {
        if(use_names) {
            setNames(tabulate(spmat@i + 1, nbins = spmat@Dim[1]), rownames(spmat))
        } else {
            tabulate(spmat@i + 1, nbins = spmat@Dim[1])
        }
    } else if(class(spmat) == 'dgRMatrix') {
        if(use_names) {
            setNames(diff(spmat@p), rownames(spmat))
        } else {
            diff(spmat@p)
        }
    } else {
        stop('Matrix should have class dgTMatrix, dgCMatrix or dgRMatrix.')
    }
}

# There already exists a function nnzero() in Matrix to find the number of nonzero elements, but
# microbenchmark shows the following is faster:

nnz <- function(spmat) {
    length(spmat@x)
    # if(class(spmat) == 'dgTMatrix' || class(spmat) == 'dgCMatrix') {
    #     length(spmat@i)
    # } else if(class(spmat) == 'dgRMatrix') {
    #     length(spmat@j)
    # } else {
    #     stop('Matrix should have class dgTMatrix, dgCMatrix or dgRMatrix.')
    # }
}

# Convert rows or columns to fractions of row/column totals:

to_frac <- function(spmat, MARGIN = 2) {
    
    MARGIN <- match.arg(as.character(MARGIN), c('1', '2'))
    
    if(MARGIN == '1') {
        
        if(class(spmat) == 'dgTMatrix' || class(spmat) == 'dgCMatrix') {
            out <- spmat
            out@x <- unname(
                unlist(
                    unname(
                        tapply( # Add names to x so that we can preserve the order of x
                            setNames(spmat@x, as.character(1:length(spmat@x))),
                            spmat@i,
                            function(y) {y/sum(y)},
                            simplify = FALSE
                        )
                    )
                )[as.character(1:length(spmat@x))]
            )
            out
        } else if(class(spmat) == 'dgRMatrix') { # Then spmat doesn't have an i attribute
            out <- spmat
            out@x <- unlist(
                unname(
                    tapply(
                        spmat@x, # Don't need to use names here, as p is ordered with respect to x
                        cut(1:length(spmat@x), spmat@p), # Will this cause problems if there are all-zero rows at the bottom of the matrix?
                        function(y) {y/sum(y)},
                        simplify = FALSE
                    )
                )
            )
            out
        } else {
            stop('Matrix should have class dgTMatrix, dgCMatrix or dgRMatrix.')
        }
        
    } else { # Then MARGIN == '2'
        
        if(class(spmat) == 'dgTMatrix' || class(spmat) == 'dgRMatrix') {
            out <- spmat
            out@x <- unname(
                unlist(
                    unname(
                        tapply(
                            setNames(spmat@x, as.character(1:length(spmat@x))),
                            spmat@j,
                            function(y) {y/sum(y)},
                            simplify = FALSE
                        )
                    )
                )[as.character(1:length(spmat@x))]
            )
            out
        } else if(class(spmat) == 'dgCMatrix') { # Then spmat doesn't have a j attribute
            out <- spmat
            out@x <- unlist(
                unname(
                    tapply(
                        spmat@x,
                        cut(1:length(spmat@x), spmat@p), # Will this cause problems if there are all-zero columns at the right of the matrix?
                        function(y) {y/sum(y)},
                        simplify = FALSE
                    )
                )
            )
            out
        } else {
            stop('Matrix should have class dgTMatrix, dgCMatrix or dgRMatrix.')
        }
        
    }
    
}

# Log transform without changing all the zeros (as '+ 1' usually does):

log_transform <- function(spmat, base = 2) {
    spmat@x <- log(spmat@x + 1, base = base)
    spmat
}
