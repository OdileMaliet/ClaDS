## ------------
## Utils functions
## ------------

# Geometric Brownian motion expectation - eq. 10 in Guindon 2013 - Systematic Biology
geometricBM <- function(par, parents, times, sigma, len){
    results <- .Call("geometricExpectation", lambda=par, parents=as.integer(parents), length=as.integer(len), sigma=sigma, brlength=times)
    return(results)
}

# Geometric Brownian motion from a Gamma distribution parameterized by the moments - eq. 10, 24 in Guindon 2013 - Systematic Biology
geometricBMGibbs <- function(par, parents, times, sigma, len){
    results <- .Call("geometricExpectationGibbs", lambda=par, parents=as.integer(parents), length=as.integer(len), sigma=sigma, brlength=times)
    return(results)
}

# Arithmetic average of nodes-rates; e.g. Kishino et al. 2001
arithmeticBM <- function(par, parents, len){
    results <- .Call("arithmetic", lambda=par, parents=as.integer(parents), length=as.integer(len))
    return(results)
}

# Simulate a normal bactrian variate (Yang & Rodriguez 2013 - PNAS)

rbactrian <- function(n, m=0.95){
    mBactrian = m
    sBactrian = sqrt(1-m^2)
    z = mBactrian + rnorm(n,0,1)*sBactrian
    rdunif <- runif(n,0,1)<0.5
    sign <- ifelse(rdunif,-1,1)
    z=z*sign
    return(z)
}

# Simulate a triangular bactrian variate (Yang & Rodriguez 2013 - PNAS)

rbactrianTriangle <- function(n, m=0.95){
    mBactrian = m
    sBactrian = sqrt(1-m^2)
    # triangle variate
    u <- runif(n,0,1)
    usign <- u<0.5
    variate = ifelse(usign,sqrt(6)-2*sqrt(3*(1-u)),-sqrt(6)+2*sqrt(3*u))
    
    z = mBactrian + variate*sBactrian
    rdunif <- runif(n,0,1)<0.5
    sign <- ifelse(rdunif,-1,1)
    z=z*sign
    return(z)
}

# Simulate a laplace bactrian variate (Yang & Rodriguez 2013 - PNAS)

rbactrianLaplace <- function(n, m=0.95){
    mBactrian = m
    sBactrian = sqrt(1-m^2)
    u <- runif(n,0,1) - 0.5
    variate = log(1-2*abs(u)) * 0.70710678118654752440
    variate = ifelse(u>=0, -variate, variate)

    z = mBactrian + variate*sBactrian
    rdunif <- runif(n,0,1)<0.5
    sign <- ifelse(rdunif,-1,1)
    z=z*sign
    return(z)
}