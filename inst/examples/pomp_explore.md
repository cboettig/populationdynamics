



```r
require(pomp)
gompertz.proc.sim <- function(x, t, params, delta.t, ...) {
    r <- params["r"]
    K <- params["K"]
    sigma <- params["sigma"]
    X <- x["X"]  # the state at time t:
    ## generate a log-normal random variable:
    eps <- exp(rnorm(n = 1, mean = 0, sd = sigma))
    ## compute the state at time t+delta.t:
    S <- exp(-r * delta.t)
    c(X = unname(K^(1 - S) * X^S * eps))
}
```







```r
gompertz.meas.sim <- function(x, t, params, ...) {
    tau <- params["tau"]
    X <- x["X"]  # state at time t:
    y <- c(Y = unname(rlnorm(n = 1, meanlog = log(X), sd = tau)))
}
```




create a container of class `pomp` to hold the model and data.



```r
gompertz <- pomp(data = data.frame(time = 1:100, Y = NA), times = "time", 
    rprocess = discrete.time.sim(step.fun = gompertz.proc.sim, delta.t = 1), 
    rmeasure = gompertz.meas.sim, t0 = 0)
```




Parameters and  inital condition:



```r
theta <- c(r = 0.1, K = 1, sigma = 0.1, tau = 0.1, X.0 = 1)
```






```r
gompertz <- simulate(gompertz, params = theta)
plot(gompertz, variables = "Y")
```

![plot of chunk simulate](http://farm9.staticflickr.com/8012/7168130004_62bc33148d_o.png) 



Likelihood using a particle filter (sequential Monte Carlo) requires the measurement density (but not the process density).  




```r
gompertz.meas.dens <- function(y, x, t, params, log, ...) {
    tau <- params["tau"]
    ## state at time t:
    X <- x["X"]
    ## observation at time t:
    Y <- y["Y"]
    ## compute the likelihood of Y|X,tau
    dlnorm(x = Y, meanlog = log(X), sdlog = tau, log = log)
}
```




Stick the new function into our container:



```r
gompertz <- pomp(gompertz, dmeasure = gompertz.meas.dens)
```




Then we get a point estimate of the likelihood,



```r
pf <- pfilter(gompertz, params = theta, Np = 1000)
logLik(pf)
```



```
[1] 49.29
```










# Iterated filtering

We will do this in the transformed variable space, so we add
our transformation method to the container:



```r
gompertz <- pomp(gompertz, parameter.transform = function(params, 
    ...) {
    exp(params)
}, parameter.inv.transform = function(params, ...) {
    log(params)
})
```





Now we're ready for an iterated-filtering run.  This is gonna be slow,
so let's set up a parallel architecture first:



```r
require(snowfall)
sfInit(parallel = TRUE, cpu = 10)
```



```
R Version:  R version 2.15.0 (2012-03-30) 

```



```r
sfLibrary(pomp)
```



```
Library pomp loaded.
```






```r
rep <- function(dummy_index) {
    estpars <- c("r", "sigma", "tau")
    theta.guess <- coef(gompertz)
    theta.guess[estpars] <- rlnorm(n = length(estpars), meanlog = log(theta.guess[estpars]), 
        sdlog = 1)
    mif(gompertz, Nmif = 100, start = theta.guess, transform = TRUE, pars = estpars, 
        rw.sd = c(r = 0.02, sigma = 0.02, tau = 0.05), Np = 2000, var.factor = 4, 
        ic.lag = 10, cooling.factor = 0.999, max.fail = 10)
}
sfExportAll()
system.time(mf <- sfSapply(1:10, rep))
```



```
   user  system elapsed 
   0.02    0.00 1037.39 
```







```r
save(list = ls(), file = "pomp_explore.rda")
```




> The key idea of iterated filtering is to replace the model we are interested in fitting—which has time invariant parameters—with a model that is just the same except that its parameters take a random walk in time. As the intensity of this random walk approaches zero, the modiﬁed model approaches the original model. 

(which it does by cooling? So convergence is already going to happen by cooling?)  


> Adding additional variability in this way has three positive effects: 
> 
> - smooths the likelihood surface, which makes optimization easier
> - it combats particle depletion, the fundamental difficulty associated with the particle filter, and 
> - the additional variability can be exploited to estimate of the gradient of the (smoothed) likelihood surface with no more computation than is required to estimate of the value of the likelihood. 

(So uncertainty estimates comes from the gradient estimate?)



Evaluate convergence 
--------------------


 just the mean of coefs across the 10 reps -- gosh wish Aaron would learn data structures



```r
theta.mif <- apply(sapply(mf, coef), 1, mean)
theta.mif
```



```
      r       K   sigma     tau     X.0 
0.04385 1.00000 0.09647 0.10495 1.00000 
```




Evaluating the log-likelihoods at the convergent parameters requires the particle filter.  This line applies the particle filter to each of the parameter estimates




```r
loglik.mif <- replicate(n = 10, logLik(pfilter(mf[[1]], params = theta.mif, 
    Np = 10000)))
bl <- mean(loglik.mif)
loglik.mif.est <- bl + log(mean(exp(loglik.mif - bl)))
loglik.mif.se <- sd(exp(loglik.mif - bl))/exp(loglik.mif.est - bl)
c(est = loglik.mif.est, se = loglik.mif.se)
```



```
     est       se 
51.45397  0.09316 
```




"True" uses the original model and coefficients rather than the estimated one, but still approximates the log.likelihood using particle filtering 



```r
loglik.true <- replicate(n = 10, logLik(pfilter(gompertz, params = coef(gompertz), 
    Np = 10000)))
loglik.true.est <- bl + log(mean(exp(loglik.true - bl)))
loglik.true.se <- sd(exp(loglik.true - bl))/exp(loglik.true.est - 
    bl)
c(est = loglik.true.est, se = loglik.true.se)
```



```
    est      se 
50.0128  0.1062 
```












