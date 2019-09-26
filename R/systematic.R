#' systematic_ra function assigns units in ordered sequence
#'
#' Binary treatment assigned
#'
#' @param p assignment probability of length n
#' @param n Number of units
#' @param m Number of units to be assigned to treatment
#'
#' @export
#'
#' @examples
#' systematic_ra(n = 10)
#' systematic_ra(p = c(.3, .4, .9))
systematic_ra <- function(p = .5, n = NULL) {
  # Housekeeping
  if(is.null(n)) n <-  length(p)
  if(length(p)!=n) {if(length(p)!=1) stop("Incompatible lengths for n and p")
    p <- rep(p, n)}

  # Dealing with p vector that does not sum to an integer
  if((sum(p) %% 1) < 0.0000001) {
    m <- sum(p)
    tag = 0
    } else {
    m <- ceiling(sum(p))
    p <- c(p, ceiling(sum(p)) - sum(p))
    n <- n+1
    tag <- 1
    }

  s   <- (cumsum(p) +m*runif(1))%%m
  e   <- s - floor(s)
  out <- 1*(e < c(e[n], e[-n]))
  if(tag==1) out <- out[-n]
  out
  }

#' prob_ra function for single treatments
#'
#' Binary treatment assigned in blocks b with probability p
#'
#' @param p assignment probability of length n
#' @param b block indicator, of length n
#' @param n Number of units
#' @param tolerance parameter for rounding errors
#'
#' @examples
#' prob_ra_internal(n = 10)
#'

prob_ra_internal <- function(
  p = .5,
  b = NULL,
  n = NULL,
  tol = 10){

  # Housekeeping

  if(is.null(n)) {if(!is.null(b)) n <- length(b)
  if( is.null(b) & length(p)>1)   n <- length(p)}
  if(length(p) == 1) p <- rep(p, n)
  if(is.null(b))     b <- rep(1, n)
  p   <- round(p, tol)

  if(ceiling(sum(p)) == 0) return(rep(0, length(p)))

  # Figure out if we have to deal with a random total

  base <- p - p%%1
  p    <- p - base

  # randomly order blocks then reorder within blocks
  b_names   <- unique(b)
  k         <- length(b_names)
  seq1      <- rep(NA, length(b))
  b_shuffle <- sample(1:k)
  for(j in 1:k) seq1[b==b_names[j]] <- b_shuffle[j]

  seq2      <- rank(seq1 + runif(n))
  p[seq2]   <- p

  # Now  do systematic assignment
  out <- systematic_ra(p)


  out <- out[seq2]
  out + base
}



#' prob_ra
#'
#' Multivalued treatment assigned
#'
#' @param p assignment probability of length n
#' @param b block indicator of length n
#' @param n Number of units
#'
#' @export
#'
#' @examples
#' prob_ra(n = 10)
#'
prob_ra <- function(p = .5,
                    b = NULL,
                    n = NULL){

  if(is.null(ncol(p))) {Z <- prob_ra_internal(p, b, n)
  } else {
    Z <- matrix(NA, nrow(p), ncol(p))
    Z[,1] <- prob_ra_internal(p[,1],b,n)

    for(j in 2:ncol(p)){
      q <- p[,j]
      q[apply(Z, 1, sum, na.rm = TRUE)==1] <- 0
      q <- q/(1-apply(as.matrix(p[,1:(j-1)]), 1, sum, na.rm = TRUE))
      q[is.nan(q)] <- 0
      Z[,j] <- prob_ra_internal(as.vector(q), b, n)
    }
    Z <-Z%*%matrix(1:ncol(p),ncol(p))
  }
  Z
  }

