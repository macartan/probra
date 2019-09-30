#' prob_ra function for single treatments
#'
#' Binary treatment assigned in blocks blocks with probability prob
#'
#' @param prob assignment probability of length n
#' @param blocks block indicator, of length n
#' @param n Number of units
#' @param tolerance parameter for rounding errors
#' @param shuffle_blocks logicial defaulting to TRUE randomly reorders blocks before applying sequential_ra
#' @param shuffle_within_blocks logicial defaulting to TRUE randomly reorders units within blocks before applying sequential_ra
#'
#' @examples
#' prob_ra_internal(n = 10)
#' # For an assignment with prob = .5 and blocks = c(1,1,1,2,2,2) we require exactly 3 units
#' # chosen each time, though sometimes 1 from block 1 and 2 from block 2 and sometimes the opposite
#'
#' reps <- replicate(1000, prob_ra_internal(blocks = c(1,1,1,2,2,2)))
#' apply(reps, 1, mean)              # .5 probabilities for each unit
#' table(apply(reps, 2, sum))        # exactly 3 selected each time
#' table(apply(reps[1:3,], 2, sum))  # 1 or 2 in block 1 each time
#'
#' # Same principles hold for non uniform probabilities
#' reps <- replicate(1000, prob_ra_internal(blocks = c(1,1,1,2,2,2), prob = (1:6)/6))
#' apply(reps, 1, mean)              # .5 probabilities for each unit
#' table(apply(reps, 2, sum))        # 3 or 4 selected each time
#' table(apply(reps[1:3,], 2, sum))  # 1 in block 1 each time
#'
#'
#' # Probabilities > 1 allowed
#' prob_ra_internal(blocks = c(1,1,1,2,2,2), prob = (1:6))    # No randomization needed
#' prob_ra_internal(blocks = c(1,1,1,2,2,2), prob = (1:6)/2)  # Randomization needed

prob_ra_internal <- function(
  prob = .5,
  blocks = NULL,
  n = NULL,
  shuffle_blocks = TRUE,
  shuffle_within_blocks = TRUE,
  tol = 10){

  # Housekeeping
  if(is.null(n)) {if(!is.null(blocks)) n <- length(blocks)
  if( is.null(blocks) & length(prob)>1)   n <- length(prob)}
  if(length(prob) == 1) prob <- rep(prob, n)
  if(is.null(blocks))     blocks <- rep(1, n)
  prob   <- round(prob, tol)

  if(ceiling(sum(prob)) == 0) return(rep(0, length(prob)))

  if(shuffle_blocks){
    # randomly order blocks then reorder within blocks
    b_names   <- unique(blocks)
    k         <- length(b_names)
    seq1      <- rep(NA, length(blocks))
    b_shuffle <- sample(1:k)
    for(j in 1:k) seq1[blocks==b_names[j]] <- b_shuffle[j]
    } else {
    seq1 = 1
    }

  if(shuffle_within_blocks){
    seq2       <- rank(seq1 + runif(n))
  } else {
    seq2       <- rank(seq1)
  }

  prob[seq2] <- prob

  # Now  do systematic assignment
  out <- systematic_ra(prob)

  # Reorder before return
  out[seq2]
}

#' prob_ra
#'
#' Complete random assignment with heterogeneous assignment probabilities and blocks.
#'
#' User provides a scalar or vector of probabilities (for sampling or assignment of a single unit) or a n*k matrix of assignment probabilities for n units to k treatments, and, optionally, indicators of block membership.
#' Assignment respects these probabilities with  minimal variance in the numbers assigned to a treatment. Though often achievable, minimal variance is not guaranteed for k>2. See vignette for more discussion.
#'
#' @param prob assignment probability that is either (a) a scalar, the same for all units (b) a vector of length n (c) a n * k matrix with probabilities for k treatments
#' @param blocks block indicator of length n
#' @param n Number of units
#' @param shuffle_blocks logicial defaulting to TRUE randomly reorders blocks before applying sequential_ra
#' @param shuffle_within_blocks logicial defaulting to TRUE randomly reorders units within blocks before applying sequential_ra
#'
#' @export
#'
#' @examples
#' prob_ra(n = 10)
#' # For an assignment with prob = .5 and blocks = c(1,1,1,2,2,2) we require exactly 3 units
#' # chosen each time, though sometimes 1 from block 1 and 2 from block 2 and sometimes the opposite
#'
#' reps <- replicate(1000, prob_ra(blocks = c(1,1,1,2,2,2)))
#' apply(reps, 1, mean)              # .5 probabilities for each unit
#' table(apply(reps, 2, sum))        # exactly 3 selected each time
#' table(apply(reps[1:3,], 2, sum))  # 1 or 2 in block 1 each time
#'
#' # Illustration with three treatments
#' # Simple case
#'  prob = matrix(1/3, 12, 3)
#'  table(prob_ra(prob = prob))
#'
#' # Uneven probabilities
#' p1 <- runif(20)
#' p2 <- runif(20)*(1-p1)
#' prob <- cbind(p1,p2, 1 - p1 - p2)
#' table(prob_ra(prob = prob))
#'
#' # If row probabilities do not sume to 1 a residual condition is defined
#' p1 <- runif(20)
#' p2 <- runif(20)*(1-p1)
#' prob <- cbind(p1,p2)
#' table(prob_ra(prob = prob))
#'
#' # Problematic cases.
#' sims <- 10000
#' prob <- t(matrix(c(.15,.65,.2, .47, .48, .05), 3,2))
#' # Results are demonstrably not as tight as possible since ideally there would be
#' # at least one unit in treatment 2 in each draw, but sometimes none here….
#' set.seed(17)
#' prob_ra(prob = prob)
#'
#' # Though assignment probabilities are still correct.
#' runs <- sapply(1:sims, function(j) prob_ra(prob = prob))
#' round(sapply(1:3, function(j) apply(runs==j, 1, sum)), 2)/sims

prob_ra <- function(prob = .5,
                    blocks = NULL,
                    n = NULL,
                    shuffle_blocks = TRUE,
                    shuffle_within_blocks = TRUE){
  # Simple case where prob is a scalar or vector
  if(is.null(ncol(prob))){
  return(prob_ra_internal(prob, blocks, n,shuffle_blocks= shuffle_blocks,  shuffle_within_blocks = shuffle_within_blocks))
  } else {

  # prob is a matrix with one column per treatment

    Z <- matrix(NA, nrow(prob), ncol(prob))
    Z[,1] <- prob_ra_internal(prob[,1], blocks,
                              n,
                              shuffle_blocks= shuffle_blocks,
                              shuffle_within_blocks = shuffle_within_blocks)

    for(j in 2:ncol(prob)){
      q <- prob[,j]
      q[apply(Z, 1, sum, na.rm = TRUE)==1] <- 0
      q <- q/(1-apply(as.matrix(prob[,1:(j-1)]), 1, sum, na.rm = TRUE))
      q[is.nan(q)] <- 0
      Z[,j] <- prob_ra_internal(as.vector(q), blocks, n)
    }
    Z <-Z%*%matrix(1:ncol(prob),ncol(prob))
  }
  Z
  }
