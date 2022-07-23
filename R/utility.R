#--------------------------------------#
#'@noRd
assert <- function (expr, error) {
  if (!expr) stop(error, call. = FALSE)
}

#' Expresses a square numeric matrix as a list of position indices and values
#'
#' @description
#' Converts a square matrix or numeric data.frame to a the list of all
#' values described with respect to the position of the indices.
#' If are present present rows and cols names the information are adds like additional columns.
#'
#' @param X matrix or data.frame
#' @param diag includes or not include the diagonal
#' @param sparse removes zero values
#' @param symmetric if true select only the upper triangular values
#'@export
#--------------------------------------#
unroll <- function(X, diag=FALSE, sparse=TRUE, symmetric=TRUE){
  #Checks
  #---------------------------#
  assert(class(X)[1]=="matrix" || class(X)[1]=="data.frame", "X must be a matrix or a data.frame")
  assert(nrow(X)==ncol(X),"X must be square")
  assert(is.logical(sparse),"sparse must be logical")
  assert(is.logical(symmetric),"symmetric must be logical")
  assert(is.logical(diag),"diag must be logical")
  assert(all(is.numeric(X)), "Found non numeric elements in X")
  #---------------------------#

  # get row and columns number; checking the diagonal or only the triangular upper
  if(symmetric==TRUE){
    index <- as.data.frame(which(upper.tri(X,diag),arr.ind=TRUE))
  } else {
    index <- as.data.frame(which(upper.tri(X,diag),arr.ind=TRUE))
    index <- rbind(index, as.data.frame(which(upper.tri(X,FALSE),arr.ind=TRUE)))
  }
  # Add the values of matrix
  values <- X[as.matrix(index)]
  # Prepare the output form
  res <- data.frame(index[,1],index[,2],values)
  colnames(res) <- c("row","col","value")
  # If are presente the rows and cols names add them to the results
  if(!is.null(colnames(X)) && !is.null(rownames(X))){
    res$row.names <- rownames(X)[index[,1]]
    res$col.names <- colnames(X)[index[,2]]
  }
  # Remove the zero values if request
  if(sparse==TRUE && length(which(res$value==0))!=0){
    res <- res[-which(res$value==0),]
  }
  # reorder the result rows similar to igraph
  res <- res[order(res$row),]
  return(res)
}


#' Upper triangular matrix
#'
#' @description
#' Returns the upper triangular values of a square matrix as one-dimensional array.
#'
#' @param X matrix or data.frame
#' @param diag includes or not the diagonal
#'@export
triu <- function(X, diag=FALSE){
  assert(class(X)[1]=="matrix" || class(X)[1]=="data.frame", "X must be a matrix or a data.frame")
  assert(nrow(X)==ncol(X),"X must be square")

  return(X[upper.tri(X,diag=diag)])
}


#' Lower triangular matrix
#'
#' @description
#' Returns the lower triangular values of a square matrix as one-dimensional array.
#'
#' @param X matrix or data.frame
#' @param diag includes or not the diagonal
#'@export
tril <- function(X, diag=FALSE){
  assert(class(X)[1]=="matrix" || class(X)[1]=="data.frame", "X must be a matrix or a data.frame")
  assert(nrow(X)==ncol(X),"X must be square")

  return(X[lower.tri(X,diag=diag)])
}


#' Apply function to columns/elements pairwise
#' Literally copied from the package AkselA/R-ymse: Ymse!!!
#'
#' Pairwise application of a function to the columns of a matrix/data.frame or
#' elements of a list
#'
#' @param x a matrix or data.frame
#' @param FUN any function that takes two vectors as input and returns a single
#' value
#' @param ... further arguments passed to FUN
#' @param comm logical; is FUN commutative? If true, only the lower
#' triangle, including the diagonal, is computed
#'
#' @return An \eqn{n\times n}{n×n} square matrix with \eqn{n} the number of columns
#' of \code{x}.
#'
#' @seealso \code{\link{similarity}} for a few more examples
#'
#' @examples
#' dtf <- data.frame(aa=c(1, 1, 2, 2, 3, 2, 4),
#'                   bb=c(1, 1, 2, 3, 3, 3, 4),
#'                   cc=c(3, 3, 2, 1, 1, 1, 1),
#'                   dd=c(1, 2, 2, 2, 1, 1, 2))
#'
#' # Root Mean Square Deviation
#' pairwise(dtf, function(x, y) sqrt(mean((x-y)^2)))
#'
#' # using with cor.test() to accompany cor()
#' pv <- pairwise(dtf, function(x, y) cor.test(x, y)$p.val)
#' pvn <- 6^(1.1-pv)-5
#' pvn[pvn<1] <- 1
#'
#' set_mar(1, 1, 1, 1)
#' plot(0, xlim=c(0.5, 4.5), ylim=c(0.5, 4.5), cex=0, ann=FALSE, xaxt="n", yaxt="n")
#' text(rep(1:4, 4), rep(4:1, each=4), t(round(cor(dtf), 2)), cex=pvn,
#'   col=c("black", "darkgrey")[(pv>0.1)+1])
#'
#'\url{https://rdrr.io/github/AkselA/R-ymse/src/R/pairwise.R}
#'
#' @export
pairwise <- function(x, FUN, ..., comm=FALSE) {

  nc <- ncol(x)
  cnames <- colnames(x)
  FUN <- match.fun(FUN)

  if (comm) {
    cb <- t(combn(nc, 2))
    nr <- nrow(cb)
    v <- vector(length=nr)
    for (i in 1:nr) {
      cc <- FUN(x[, cb[i,1]], x[, cb[i,2]], ...)
      v[i] <- cc
    }
    m <- matrix(1, nc, nc)
    m[lower.tri(m)] <- v
    m <- m*t(m)
    diag(m) <- sapply(1:nc, function(j) FUN(x[,j], x[,j], ...))
  } else {
    eg <- expand.grid(1:nc, 1:nc)
    nr <- nrow(eg)
    v <- vector(length=nr)
    for (i in 1:nr) {
      cc <- FUN(x[,eg[i,1]], x[,eg[i,2]], ...)
      v[i] <- cc
    }
    m <- matrix(v, nc, byrow=TRUE)
  }
  dimnames(m) <- list(cnames, cnames)
  m
}

