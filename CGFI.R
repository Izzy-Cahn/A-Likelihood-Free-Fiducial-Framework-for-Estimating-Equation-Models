############################################################
# CGFI Simulation Study
# Likelihood-Free Fiducial Inference Using Gaussian Copula
############################################################

library(MASS)
library(ggplot2)

set.seed(123)



############################################################
# Moment function
############################################################

moments <- function(theta, Y, X, Z){
  
  residual <- Y - theta*X - 0.5*X^2
  
  Z * residual
  
}



############################################################
# Pre-compute moment components
#
# Since:
#
# g(theta)=Z(Y-theta*X-0.5X^2)
#
# the sample moments can be written as:
#
# gbar(theta)=ZY - theta*ZX - ZX2
#
############################################################

prepare_moments <- function(Y, X, Z){
  
  list(
    
    ZY = colMeans(Z * Y),
    
    ZX = colMeans(Z * X),
    
    ZX2 = colMeans(Z * (0.5*X^2)),
    
    n = length(Y),
    
    q = ncol(Z)
    
  )
  
}



############################################################
# Fast sample moment function
############################################################

moment_means <- function(theta, prep){
  
  prep$ZY -
    theta*prep$ZX -
    prep$ZX2
  
}



############################################################
# Nonlinear GMM estimator
############################################################

gmm_est <- function(prep){
  
  objective <- function(theta){
    
    g <- moment_means(
      theta,
      prep
    )
    
    sum(g^2)
    
  }
  
  
  optimize(
    objective,
    interval = c(0,4),
    tol = 1e-8
  )$minimum
  
}



############################################################
# Estimate covariance matrix of moments
############################################################

estimate_Omega <- function(theta_hat,
                           Y,
                           X,
                           Z){
  
  
  g_mat <- moments(
    theta_hat,
    Y,
    X,
    Z
  )
  
  
  Omega <- cov(g_mat)
  
  
  # Regularization for high-dimensional moments
  
  lambda <- 0.01
  
  Omega <-
    (1-lambda)*Omega +
    lambda*diag(ncol(Omega))
  
  
  Omega
  
}



############################################################
# Prepare Gaussian copula distribution
############################################################

prepare_copula <- function(Omega){
  
  
  scales <- sqrt(
    diag(Omega)
  )
  
  
  R <- Omega /
    outer(scales, scales)
  
  
  R <- (R+t(R))/2
  
  
  L <- chol(R)
  
  
  list(
    
    scales = scales,
    
    chol = L
    
  )
  
  
}



############################################################
# Generate Gaussian copula moment shock
############################################################

cgfi_shock <- function(copula){
  
  
  z <- copula$chol %*%
    rnorm(
      length(copula$scales)
    )
  
  
  eta <- copula$scales*z
  
  
  as.numeric(eta)
  
}

############################################################
# Fast CGFI draw generator
############################################################

cgfi_draw <- function(theta_hat,
                      prep,
                      copula){
  
  
  n <- prep$n
  
  
  # Generate Gaussian copula perturbation
  
  eta <- cgfi_shock(
    copula
  )
  
  
  target <- -eta/sqrt(n)
  
  
  # Solve perturbed estimating equation
  
  objective <- function(theta){
    
    g <- moment_means(
      theta,
      prep
    )
    
    
    sum(
      (g-target)^2
    )
    
  }
  
  
  optimize(
    objective,
    interval=c(0,4),
    tol=1e-8
  )$minimum
  
}



############################################################
# Wald standard error
############################################################

wald_interval <- function(theta_hat,
                          Omega,
                          prep){
  
  
  # numerical derivative of moments
  
  eps <- 1e-5
  
  
  G <- (
    
    moment_means(
      theta_hat + eps,
      prep
    )
    -
      moment_means(
        theta_hat - eps,
        prep
      )
    
  )/(2*eps)
  
  
  G <- matrix(
    G,
    ncol=1
  )
  
  
  W <- solve(Omega)
  
  
  V <- solve(
    t(G)%*%W%*%G
  )/prep$n
  
  
  se <- sqrt(
    as.numeric(V)
  )
  
  
  theta_hat +
    c(-1,1)*
    1.96*
    se
  
}



############################################################
# Simulation function
############################################################

simulate_q <- function(q,
                       R=100,
                       B=200,
                       n=500,
                       theta0=2){
  
  
  coverage <- matrix(
    0,
    R,
    3
  )
  
  
  lengths <- matrix(
    0,
    R,
    3
  )
  
  
  boot_time <- numeric(R)
  cgfi_time <- numeric(R)
  
  
  
  for(r in 1:R){
    
    
    ########################################################
    # Generate data
    ########################################################
    
    X <- rnorm(n)
    
    
    Z <- matrix(
      rnorm(n*q),
      nrow=n,
      ncol=q
    )
    
    
    u <- rt(
      n,
      df=5
    )
    
    
    Y <- theta0*X +
      0.5*X^2 +
      u
    
    
    ########################################################
    # Estimate GMM parameter
    ########################################################
    
    prep <- prepare_moments(
      Y,
      X,
      Z
    )
    
    
    theta_hat <- gmm_est(
      prep
    )
    
    
    ########################################################
    # Estimate covariance
    ########################################################
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    
    copula <- prepare_copula(
      Omega
    )
    
    
    
    ########################################################
    # Wald interval
    ########################################################
    
    ci <- wald_interval(
      theta_hat,
      Omega,
      prep
    )
    
    
    coverage[r,1] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    
    lengths[r,1] <-
      diff(ci)
    
    
    
    ########################################################
    # Bootstrap
    ########################################################
    
    start <- Sys.time()
    
    
    boot <- numeric(B)
    
    
    for(b in 1:B){
      
      
      id <- sample(
        1:n,
        n,
        replace=TRUE
      )
      
      
      prep_b <- prepare_moments(
        Y[id],
        X[id],
        Z[id,]
      )
      
      
      boot[b] <- gmm_est(
        prep_b
      )
      
      
    }
    
    
    boot_time[r] <-
      as.numeric(
        Sys.time()-start,
        units="secs"
      )
    
    
    ci <- quantile(
      boot,
      c(.025,.975)
    )
    
    
    coverage[r,2] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    
    lengths[r,2] <-
      diff(ci)
    
    
    
    ########################################################
    # CGFI
    ########################################################
    
    start <- Sys.time()
    
    
    cgfi <- numeric(B)
    
    
    for(b in 1:B){
      
      cgfi[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula
      )
      
    }
    
    
    cgfi_time[r] <-
      as.numeric(
        Sys.time()-start,
        units="secs"
      )
    
    
    ci <- quantile(
      cgfi,
      c(.025,.975)
    )
    
    
    coverage[r,3] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    
    lengths[r,3] <-
      diff(ci)
    
    
    
    if(r %% 10 == 0){
      
      cat(
        "Completed",
        r,
        "of",
        R,
        "replications\n"
      )
      
    }
    
    
  }
  
  
  
  data.frame(
    
    q=q,
    
    Method=c(
      "Wald",
      "Bootstrap",
      "CGFI"
    ),
    
    Coverage=c(
      mean(coverage[,1]),
      mean(coverage[,2]),
      mean(coverage[,3])
    ),
    
    Avg_Length=c(
      mean(lengths[,1]),
      mean(lengths[,2]),
      mean(lengths[,3])
    ),
    
    Avg_Time=c(
      NA,
      mean(boot_time),
      mean(cgfi_time)
    )
    
  )
  
}

############################################################
# Experiment 1:
# Increasing number of moment conditions
############################################################


q_values <- c(
  10,
  50,
  100,
  250
)


q_results <- do.call(
  rbind,
  lapply(
    q_values,
    function(q){
      
      simulate_q(
        q=q,
        n=500,
        R=100,
        B=200
      )
      
    }
  )
)


print(q_results)



############################################################
# Coverage plot
############################################################


ggplot(
  q_results,
  aes(
    x=q,
    y=Coverage,
    group=Method,
    color=Method
  )
)+
  geom_line(
    linewidth=1
  )+
  geom_point(
    size=2
  )+
  geom_hline(
    yintercept=0.95,
    linetype="dashed"
  )+
  ylim(
    0.5,
    1
  )+
  theme_minimal()+
  labs(
    x="Number of Moment Conditions (q)",
    y="Coverage Probability",
    color="Method"
  )



############################################################
# Interval length plot
############################################################


ggplot(
  q_results,
  aes(
    x=q,
    y=Avg_Length,
    group=Method,
    color=Method
  )
)+
  geom_line(
    linewidth=1
  )+
  geom_point(
    size=2
  )+
  theme_minimal()+
  labs(
    x="Number of Moment Conditions (q)",
    y="Average Confidence Interval Length",
    color="Method"
  )



############################################################
# Experiment 2:
# Increasing sample size
# Holding q/n = 0.2
############################################################


n_values <- c(
  100,
  250,
  500,
  1000
)


n_results <- do.call(
  rbind,
  lapply(
    n_values,
    function(n){
      
      simulate_q(
        q=floor(0.2*n),
        n=n,
        R=100,
        B=200
      )
      
    }
  )
)


print(n_results)



############################################################
# Runtime experiment only
############################################################


runtime_results <- q_results[
  ,
  c(
    "q",
    "Method",
    "Avg_Time"
  )
]


print(runtime_results)