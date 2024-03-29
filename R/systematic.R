#' systematic_ra function assigns units in ordered sequence
#'
#' Usual use is to sample (or assign to a binary treatment) a set of units given possibly heterogeneous sampling (assignment) probabilities.
#' Systematic sampling uses the ordering of units provided, takes an initial random starting point, and selects at fixed intervals over the list, with each unit selected according to its specified probability.
#' For fully random assignment the order of lists should be randomized first, as done for example by \code{prob_ra}.
#' In cases where probabilities exceed 1, all units are first assigned their "guaranteed" amounts and then residual amounts are randomized across units. This might be used for example for assigning
#' targets to blocks rather than individuals to a treatment.
#'
#'
#' @param prob assignment probability of length n
#' @param n Number of units
#' @param m Number of units to be assigned to treatment
#'
#' @export
#'
#' @examples
#' # Simple aaplication: note fixed gap between selected units
#' systematic_ra(n = 20)
#' systematic_ra(prob = rep(1/3, 20))
#' systematic_ra(prob = rep(2/7, 20))
#'
#' # Nonuniform probabilities
#' # For an assignment with prob = c(.3, .4, .9) we need
#' # (a) individual assignments according to prob and
#' # (b) the total number assigned should be either 1 or 2
#' systematic_ra(prob = c(.3, .4, .9))
#' reps <- replicate(1000, systematic_ra(prob = c(.3, .4, .9)))
#' apply(reps, 1, mean)
#' table(apply(reps, 2, sum))
#'
#' # Probabilities > 1 allowed
#' systematic_ra(prob = c(.5, 1.5))
#' reps <- replicate(500, systematic_ra(prob = c(.5, 1.5)))
#' apply(reps, 1, mean)

systematic_ra <- function(prob = .5, n = NULL) {

  # Housekeeping
  if(is.null(n)) n <-  length(prob)
  if(length(prob)!=n) {if(length(prob)!=1) stop("Incompatible lengths for n and prob")
    prob <- rep(prob, n)}

  # Probabilities > 1 given base (implies assignments of more than 1 to unit)
  base <- prob - prob%%1
  prob <- prob - base
  if(sum(prob) == 0) return(base)

  # Dealing with prob vector that does not sum to an integer
  if((sum(prob) %% 1) < 0.0000001) {
    m <- sum(prob)
    tag = 0
    } else {
    m <- ceiling(sum(prob))
    prob <- c(prob, ceiling(sum(prob)) - sum(prob))
    n <- n+1
    tag <- 1
    }

  s   <- (cumsum(prob) +m*runif(1))%%m
  e   <- s - floor(s)
  out <- 1*(e < c(e[n], e[-n]))
  if(tag==1) out <- out[-n]

  out + base
  }

