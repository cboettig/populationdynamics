\name{saddle_node_ibm}
\alias{saddle_node_ibm}
\title{Individual based simulation of a saddle node bifurcation

Uses the Gillespie algorithm to simulate an individual birth-death process
for the saddle node bifurcation.  The system gradually approaches the bifurcation
at the specified rate.}
\usage{
  saddle_node_ibm(pars = c(Xo = 570, e = 0.5, a = 160, K = 1000, h = 200, i = 0, Da = 1, Dt = 100, p = 2),
    times = seq(0, 150, length = 50), reps = 1)
}
\arguments{
  \item{pars}{a list of parameters for the saddle-node
  model}

  \item{times}{a sequence of times at which we sample the
  system state (note that the dynamics of the system itself
  are continuous and independent of this sampling process.}

  \item{reps}{how many replicate simulations should we
  run?}
}
\value{
  a list with the observed values in the matrix as "x1"
  (columns are reps, if desired), mean values "m1" across
  replicates, "v1" variance across replicates, "parameters"
  is the input list of pars, and "time" is the input list
  of times.
}
\description{
  Individual based simulation of a saddle node bifurcation

  Uses the Gillespie algorithm to simulate an individual
  birth-death process for the saddle node bifurcation.  The
  system gradually approaches the bifurcation at the
  specified rate.
}
\details{
  If given replicates, this produces some summary
  statistics of the replicates as well.

  pars contains a named list of these elements, in this
  order: Xo is initial population size e is the natural
  per-capita death-rate a is the environmental toxin level
  K scales the birth rate (hence the equilibrium size) h is
  the half-max growth rate i is a place-holder for an
  internal counter, not real parameter Da is the rate of
  environmental degradation Dt is the time at which
  environmental degradation begins
}
\examples{
pars = c(Xo = 730, e = 0.5, a = 100, K = 1000, h = 200, i = 0, Da = .09, Dt = 0, p = 2)
time=seq(0,500, length=500)
sn <- saddle_node_ibm(pars,time)
X <- data.frame(time=time, value=sn$x1)
}

