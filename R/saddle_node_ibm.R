# @file: ind_based_models.R
# @author: Carl Boettiger, <cboettig@gmail.com>
# @section DESCRIPTION wrapper to the C code containing the gillespie simulation for individual-based models.  

## Some example parameters to try:

## In the C code
# Pars = {n, e, a, K, h, i, Da, Dt} 
# inits[8] = {572, .5, 160, 1000, 200, 0, 1, 100};

# in the R code
#pars = c(Xo = 730, e = 0.5, a = 100, K = 1000, h = 200,
#         i = 0, Da = .09, Dt = 0, p = 2)



#' @title Individual based simulation of a saddle node bifurcation
#'
#' Uses the Gillespie algorithm to simulate an individual birth-death process
#' for the saddle node bifurcation.  The system gradually approaches the bifurcation
#' at the specified rate.  
#' 
#' @param pars a list of parameters for the saddle-node model
#' @param times a sequence of times at which we sample the system state 
#' (note that the dynamics of the system itself are continuous and independent
#' of this sampling process.  
#' @param reps how many replicate simulations should we run?
#' @return a list with the observed values in the matrix as "x1" (columns are reps, if desired),
#' mean values "m1" across replicates, "v1" variance across replicates, "parameters" is 
#' the input list of pars, and "time" is the input list of times. 
#' 
#' @details 
#'
#' If given replicates, this produces some summary statistics of the replicates as well.  
#' 
#' pars contains a named list of these elements, in this order:
#' Xo is initial population size
#' e is the natural per-capita death-rate
#' a is the environmental toxin level
#' K scales the birth rate (hence the equilibrium size)
#' h is the half-max growth rate
#' i is a place-holder for an internal counter, not real parameter
#' Da is the rate of environmental degradation
#' Dt is the time at which environmental degradation begins
#'
#' @examples
#' pars = c(Xo = 730, e = 0.5, a = 100, K = 1000, h = 200, i = 0, Da = .09, Dt = 0, p = 2)
#' time=seq(0,500, length=500)
#' sn <- saddle_node_ibm(pars,time)
#' X <- data.frame(time=time, value=sn$x1)
#'
#' @useDynLib populationdynamics 
#' @export
saddle_node_ibm <- function(
    pars=c("Xo" = 570, "e" = .5, "a" = 160, "K" = 1000, "h" = 200,
           "i" = 0, "Da" = 1, "Dt" = 100, "p"=2),
    times = seq(0,150,length=50), 
    reps=1)
{
    samples <- length(times)
    N <- reps*samples
    maxtime <- max(times)

    o <- .C("saddle_node_direct", double(N), as.double(pars), as.integer(samples), as.integer(reps), as.double(maxtime) )

    x1 <- matrix(o[[1]], samples, reps)
    m1 <- sapply(1:samples, function(i) mean(x1[i,]))
    v1 <- sapply(1:samples, function(i) var(x1[i,]))

    list(x1 = x1,  m1=m1, v1=v1, parameters = pars, Xo = pars[1], time=times)
}



