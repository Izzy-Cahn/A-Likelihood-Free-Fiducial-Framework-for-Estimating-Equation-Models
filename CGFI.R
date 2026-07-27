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
  
  residual <- Y - theta * X - 0.5 * X^2
  
  Z * residual
  
}

############################################################
# Pre-compute sample moment quantities
############################################################

prepare_moments <- function(Y, X, Z){
  
  list(
    
    ZY  = colMeans(Z * Y),
    
    ZX  = colMeans(Z * X),
    
    ZX2 = colMeans(Z * (0.5 * X^2)),
    
    n = length(Y),
    
    q = ncol(Z)
    
  )
  
}

############################################################
# Mean moment vector
############################################################

moment_means <- function(theta, prep){
  
  prep$ZY -
    theta * prep$ZX -
    prep$ZX2
  
}

############################################################
# Nonlinear GMM estimator
############################################################

gmm_est <- function(prep, W = NULL){
  
  if(is.null(W))
    W <- diag(prep$q)
  
  objective <- function(theta){
    
    g <- moment_means(theta, prep)
    
    as.numeric(
      t(g) %*% W %*% g
    )
    
  }
  
  optimize(
    objective,
    interval = c(0,4),
    tol = 1e-8
  )$minimum
  
}

############################################################
# Estimate covariance of estimating equations
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
  
  # Small ridge regularization
  
  lambda <- 0.01
  
  Omega <-
    (1 - lambda) * Omega +
    lambda * diag(ncol(Omega))
  
  Omega
  
}

############################################################
# Prepare Gaussian copula
############################################################

prepare_copula <- function(Omega){
  
  scales <- sqrt(diag(Omega))
  
  R <- Omega /
    outer(scales, scales)
  
  R[!is.finite(R)] <- 0
  
  diag(R) <- 1
  
  R <- (R + t(R))/2
  
  eig <- eigen(R)
  
  eig$values[eig$values < 1e-8] <- 1e-8
  
  R <- eig$vectors %*%
    diag(eig$values) %*%
    t(eig$vectors)
  
  L <- chol(R)
  
  list(
    
    scales = scales,
    
    chol = L
    
  )
  
}

############################################################
# Generate Gaussian copula perturbation
############################################################

cgfi_shock <- function(copula){
  
  z <- copula$chol %*%
    rnorm(length(copula$scales))
  
  eta <- copula$scales * z
  
  as.numeric(eta)
  
}

############################################################
# Generate one CGFI draw
############################################################

cgfi_draw <- function(theta_hat,
                      prep,
                      copula,
                      W){
  
  eta <- cgfi_shock(copula)
  
  target <- -eta / sqrt(prep$n)
  
  objective <- function(theta){
    
    g <- moment_means(theta, prep)
    
    as.numeric(
      t(g - target) %*%
        W %*%
        (g - target)
    )
    
  }
  
  optimize(
    objective,
    interval = c(0,4),
    tol = 1e-8
  )$minimum
  
}

############################################################
# Wald confidence interval
############################################################

wald_interval <- function(theta_hat,
                          Omega,
                          prep){
  
  eps <- 1e-5
  
  G <- (
    
    moment_means(
      theta_hat + eps,
      prep
    ) -
      
      moment_means(
        theta_hat - eps,
        prep
        
      )
    
  ) / (2 * eps)
  
  G <- matrix(
    G,
    ncol = 1
  )
  
  W <- solve(Omega)
  
  V <- solve(
    t(G) %*%
      W %*%
      G
  ) / prep$n
  
  se <- sqrt(as.numeric(V))
  
  theta_hat +
    c(-1,1) *
    1.96 *
    se
  
}

############################################################
# Experiment 1:
# Increasing number of moment conditions
############################################################

simulate_q <- function(q,
                       R = 1000,
                       B = 500,
                       n = 500,
                       theta0 = 2){
  
  coverage <- matrix(0, R, 3)
  
  lengths <- matrix(0, R, 3)
  
  boot_time <- numeric(R)
  
  cgfi_time <- numeric(R)
  
  theta_store <- rep(NA, R)
  
  failures <- 0
  
  for(r in 1:R){
    
    ########################################################
    # Generate data
    ########################################################
    
    X <- rnorm(n)
    
    Z <- matrix(
      rnorm(n * q),
      nrow = n,
      ncol = q
    )
    
    u <- rt(
      n,
      df = 5
    )
    
    Y <- theta0 * X +
      0.5 * X^2 +
      u
    
    ########################################################
    # Prepare moments
    ########################################################
    
    prep <- prepare_moments(
      Y,
      X,
      Z
    )
    
    ########################################################
    # Two-step GMM
    ########################################################
    
    theta_hat <- tryCatch({
      
      theta_initial <- gmm_est(prep)
      
      Omega <- estimate_Omega(
        theta_initial,
        Y,
        X,
        Z
      )
      
      W <- solve(Omega)
      
      theta_final <- gmm_est(
        prep,
        W
      )
      
      theta_final
      
    }, error = function(e){
      
      failures <<- failures + 1
      
      NA
      
    })
    
    if(is.na(theta_hat))
      next
    
    theta_store[r] <- theta_hat
    
    ########################################################
    # Update covariance at final estimate
    ########################################################
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    W <- solve(Omega)
    
    copula <- prepare_copula(
      Omega
    )
    
    ########################################################
    # Wald
    ########################################################
    
    ci <- wald_interval(
      theta_hat,
      Omega,
      prep
    )
    
    coverage[r,1] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,1] <- diff(ci)
    
    ########################################################
    # Bootstrap
    ########################################################
    
    start <- Sys.time()
    
    boot <- numeric(B)
    
    for(b in 1:B){
      
      id <- sample(
        1:n,
        n,
        replace = TRUE
      )
      
      prep_b <- prepare_moments(
        Y[id],
        X[id],
        Z[id,]
      )
      
      theta_b0 <- gmm_est(
        prep_b
      )
      
      Omega_b <- estimate_Omega(
        theta_b0,
        Y[id],
        X[id],
        Z[id,]
      )
      
      W_b <- solve(Omega_b)
      
      boot[b] <- gmm_est(
        prep_b,
        W_b
      )
      
    }
    
    boot_time[r] <-
      as.numeric(
        Sys.time() - start,
        units = "secs"
      )
    
    ci <- quantile(
      boot,
      c(.025,.975)
    )
    
    coverage[r,2] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,2] <- diff(ci)
    
    ########################################################
    # CGFI
    ########################################################
    
    start <- Sys.time()
    
    cgfi <- numeric(B)
    
    for(b in 1:B){
      
      cgfi[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula,
        W
      )
      
    }
    
    cgfi_time[r] <-
      as.numeric(
        Sys.time() - start,
        units = "secs"
      )
    
    ci <- quantile(
      cgfi,
      c(.025,.975)
    )
    
    coverage[r,3] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,3] <- diff(ci)
    
    if(r %% 25 == 0){
      
      cat(
        "Completed",
        r,
        "of",
        R,
        "replications\n"
      )
      
    }
    
  }
  
  ##########################################################
  # Summary statistics
  ##########################################################
  
  valid <- !is.na(theta_store)
  
  theta_store <- theta_store[valid]
  
  cover <- c(
    mean(coverage[valid,1]),
    mean(coverage[valid,2]),
    mean(coverage[valid,3])
  )
  
  cover_se <- sqrt(
    cover * (1 - cover) /
      length(theta_store)
  )
  
  data.frame(
    
    q = q,
    
    Method = c(
      "Wald",
      "Bootstrap",
      "CGFI"
    ),
    
    Coverage = cover,
    
    Coverage_SE = cover_se,
    
    Avg_Length = c(
      mean(lengths[valid,1]),
      mean(lengths[valid,2]),
      mean(lengths[valid,3])
    ),
    
    Avg_Time = c(
      NA,
      mean(boot_time[valid]),
      mean(cgfi_time[valid])
    ),
    
    Speedup = c(
      NA,
      NA,
      mean(boot_time[valid]) /
        mean(cgfi_time[valid])
    ),
    
    Mean_Estimate = rep(
      mean(theta_store),
      3
    ),
    
    Bias = rep(
      mean(theta_store - theta0),
      3
    ),
    
    RMSE = rep(
      sqrt(
        mean((theta_store - theta0)^2)
      ),
      3
    ),
    
    Failures = rep(
      failures,
      3
    )
    
  )
  
}

############################################################
# Simulation function
############################################################

simulate_q <- function(q,
                       R = 100,
                       B = 200,
                       n = 500,
                       theta0 = 2){
  
  coverage <- matrix(0, R, 3)
  lengths  <- matrix(0, R, 3)
  
  theta_hat_store <- numeric(R)
  
  boot_time <- numeric(R)
  cgfi_time <- numeric(R)
  
  for(r in 1:R){
    
    ########################################################
    # Generate data
    ########################################################
    
    X <- rnorm(n)
    
    Z <- matrix(
      rnorm(n*q),
      nrow = n,
      ncol = q
    )
    
    u <- rt(n, df = 5)
    
    Y <- theta0*X +
      0.5*X^2 +
      u
    
    ########################################################
    # Two-step GMM estimator
    ########################################################
    
    prep <- prepare_moments(
      Y,
      X,
      Z
    )
    
    theta_initial <- gmm_est(prep)
    
    Omega <- estimate_Omega(
      theta_initial,
      Y,
      X,
      Z
    )
    
    W <- solve(Omega)
    
    theta_hat <- gmm_est(
      prep,
      W
    )
    
    theta_hat_store[r] <- theta_hat
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    W <- solve(Omega)
    
    copula <- prepare_copula(
      Omega
    )
    
    ########################################################
    # Wald
    ########################################################
    
    ci <- wald_interval(
      theta_hat,
      Omega,
      prep
    )
    
    coverage[r,1] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,1] <- diff(ci)
    
    ########################################################
    # Bootstrap
    ########################################################
    
    start <- Sys.time()
    
    boot <- numeric(B)
    
    for(b in 1:B){
      
      id <- sample(
        1:n,
        n,
        replace = TRUE
      )
      
      prep_b <- prepare_moments(
        Y[id],
        X[id],
        Z[id,]
      )
      
      theta_b0 <- gmm_est(prep_b)
      
      Omega_b <- estimate_Omega(
        theta_b0,
        Y[id],
        X[id],
        Z[id,]
      )
      
      W_b <- solve(Omega_b)
      
      boot[b] <- gmm_est(
        prep_b,
        W_b
      )
      
    }
    
    boot_time[r] <-
      as.numeric(
        Sys.time() - start,
        units = "secs"
      )
    
    ci <- quantile(
      boot,
      c(.025,.975)
    )
    
    coverage[r,2] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,2] <- diff(ci)
    
    ########################################################
    # CGFI
    ########################################################
    
    start <- Sys.time()
    
    cgfi <- numeric(B)
    
    for(b in 1:B){
      
      cgfi[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula,
        W
      )
      
    }
    
    cgfi_time[r] <-
      as.numeric(
        Sys.time() - start,
        units = "secs"
      )
    
    ci <- quantile(
      cgfi,
      c(.025,.975)
    )
    
    coverage[r,3] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,3] <- diff(ci)
    
    if(r %% 10 == 0){
      
      cat(
        "Completed",
        r,
        "of",
        R,
        "\n"
      )
      
    }
    
  }
  
  cover <- c(
    mean(coverage[,1]),
    mean(coverage[,2]),
    mean(coverage[,3])
  )
  
  data.frame(
    
    q = q,
    
    Method = c(
      "Wald",
      "Bootstrap",
      "CGFI"
    ),
    
    Coverage = cover,
    
    MCSE = sqrt(
      cover*(1-cover)/R
    ),
    
    Avg_Length = c(
      mean(lengths[,1]),
      mean(lengths[,2]),
      mean(lengths[,3])
    ),
    
    Avg_Time = c(
      NA,
      mean(boot_time),
      mean(cgfi_time)
    ),
    
    Bias = rep(
      mean(theta_hat_store - theta0),
      3
    ),
    
    RMSE = rep(
      sqrt(
        mean(
          (theta_hat_store-theta0)^2
        )
      ),
      3
    )
    
  )
  
}

############################################################
# Experiment 1
# Increasing number of moments
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
        q = q,
        n = 500,
        R = 100,
        B = 200
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
    x = q,
    y = Coverage,
    color = Method,
    group = Method
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = Coverage - 1.96*MCSE,
      ymax = Coverage + 1.96*MCSE
    ),
    width = 4
  ) +
  geom_hline(
    yintercept = 0.95,
    linetype = "dashed"
  ) +
  ylim(0.5,1) +
  theme_minimal() +
  labs(
    x = "Number of Moment Conditions",
    y = "Coverage Probability"
  )

############################################################
# Interval length plot
############################################################

ggplot(
  q_results,
  aes(
    x = q,
    y = Avg_Length,
    color = Method,
    group = Method
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    x = "Number of Moment Conditions",
    y = "Average Interval Length"
  )

############################################################
# Experiment 2
# Increasing sample size
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
        q = floor(0.2*n),
        n = n,
        R = 100,
        B = 200
      )
      
    }
  )
)

print(n_results)

############################################################
# Runtime summary
############################################################

runtime_results <-
  q_results[
    ,
    c(
      "q",
      "Method",
      "Avg_Time"
    )
  ]

print(runtime_results)

############################################################
############################################################
# Experiment 3:
# Robustness to dependent moment conditions
############################################################
############################################################


simulate_dependence <- function(q,
                                rho = 0.5,
                                R = 100,
                                B = 200,
                                n = 500,
                                theta0 = 2){
  
  
  coverage <- matrix(0,R,3)
  lengths  <- matrix(0,R,3)
  
  boot_time <- numeric(R)
  cgfi_time <- numeric(R)
  
  
  for(r in 1:R){
    
    
    ########################################################
    # Generate correlated instruments
    ########################################################
    
    
    Sigma_Z <- matrix(
      rho,
      q,
      q
    )
    
    diag(Sigma_Z) <- 1
    
    
    L <- chol(Sigma_Z)
    
    
    Z <- matrix(
      rnorm(n*q),
      nrow=n,
      ncol=q
    ) %*% L
    
    
    
    X <- rnorm(n)
    
    u <- rt(
      n,
      df=5
    )
    
    
    Y <- theta0*X +
      0.5*X^2 +
      u
    
    
    
    ########################################################
    # Two-step GMM
    ########################################################
    
    
    prep <- prepare_moments(
      Y,
      X,
      Z
    )
    
    
    theta_initial <- gmm_est(
      prep
    )
    
    
    Omega <- estimate_Omega(
      theta_initial,
      Y,
      X,
      Z
    )
    
    
    W <- solve(Omega)
    
    
    theta_hat <- gmm_est(
      prep,
      W
    )
    
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    
    W <- solve(Omega)
    
    
    copula <- prepare_copula(
      Omega
    )
    
    
    
    ########################################################
    # Wald
    ########################################################
    
    
    ci <- wald_interval(
      theta_hat,
      Omega,
      prep
    )
    
    
    coverage[r,1] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    
    lengths[r,1] <- diff(ci)
    
    
    
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
      
      
      theta_b0 <- gmm_est(
        prep_b
      )
      
      
      Omega_b <- estimate_Omega(
        theta_b0,
        Y[id],
        X[id],
        Z[id,]
      )
      
      
      W_b <- solve(Omega_b)
      
      
      boot[b] <- gmm_est(
        prep_b,
        W_b
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
    
    
    lengths[r,2] <- diff(ci)
    
    
    
    ########################################################
    # CGFI
    ########################################################
    
    
    start <- Sys.time()
    
    
    cgfi <- numeric(B)
    
    
    for(b in 1:B){
      
      
      cgfi[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula,
        W
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
    
    
    lengths[r,3] <- diff(ci)
    
    
    
    if(r %% 10 == 0){
      
      cat(
        "Completed",
        r,
        "of",
        R,
        "\n"
      )
      
    }
    
    
  }
  
  
  
  cover <- c(
    mean(coverage[,1]),
    mean(coverage[,2]),
    mean(coverage[,3])
  )
  
  
  data.frame(
    
    rho=rho,
    
    Method=c(
      "Wald",
      "Bootstrap",
      "CGFI"
    ),
    
    Coverage=cover,
    
    MCSE=sqrt(
      cover*(1-cover)/R
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
# Run dependence experiment
############################################################


dependence_results <- do.call(
  rbind,
  lapply(
    c(0,0.25,0.50,0.75),
    function(rho){
      
      simulate_dependence(
        q=100,
        rho=rho,
        n=500,
        R=100,
        B=200
      )
      
    }
  )
)


print(dependence_results)



############################################################
# Dependence coverage plot
############################################################


ggplot(
  dependence_results,
  aes(
    x=rho,
    y=Coverage,
    color=Method,
    group=Method
  )
)+
  geom_line(linewidth=1)+
  geom_point(size=2)+
  geom_errorbar(
    aes(
      ymin=Coverage-1.96*MCSE,
      ymax=Coverage+1.96*MCSE
    ),
    width=.03
  )+
  geom_hline(
    yintercept=.95,
    linetype="dashed"
  )+
  ylim(.5,1)+
  theme_minimal()+
  labs(
    x="Correlation Among Moment Conditions",
    y="Coverage Probability"
  )



############################################################
############################################################
# Experiment 4:
# Alternative Error Distributions
############################################################
############################################################

simulate_errors <- function(error_type,
                            q=100,
                            R=100,
                            B=200,
                            n=500,
                            theta0=2){
  
  coverage <- matrix(0, R, 3)
  lengths  <- matrix(0, R, 3)
  
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
    
    ########################################################
    # Error distributions
    ########################################################
    
    if(error_type=="Normal"){
      u <- rnorm(n)
    }
    
    if(error_type=="t5"){
      u <- rt(n, df=5)
    }
    
    if(error_type=="Heavy"){
      u <- rt(n, df=3)
    }
    
    Y <- theta0*X + 0.5*X^2 + u
    
    ########################################################
    # Two-step GMM
    ########################################################
    
    prep <- prepare_moments(Y,X,Z)
    
    theta_initial <- gmm_est(prep)
    
    Omega <- estimate_Omega(
      theta_initial,
      Y,
      X,
      Z
    )
    
    W <- solve(Omega)
    
    theta_hat <- gmm_est(
      prep,
      W
    )
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    W <- solve(Omega)
    
    copula <- prepare_copula(Omega)
    
    ########################################################
    # Wald
    ########################################################
    
    ci <- wald_interval(
      theta_hat,
      Omega,
      prep
    )
    
    coverage[r,1] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,1] <- diff(ci)
    
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
      
      theta0_b <- gmm_est(prep_b)
      
      Omega_b <- estimate_Omega(
        theta0_b,
        Y[id],
        X[id],
        Z[id,]
      )
      
      W_b <- solve(Omega_b)
      
      boot[b] <- gmm_est(
        prep_b,
        W_b
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
    
    lengths[r,2] <- diff(ci)
    
    ########################################################
    # CGFI
    ########################################################
    
    start <- Sys.time()
    
    draws <- numeric(B)
    
    for(b in 1:B){
      
      draws[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula,
        W
      )
      
    }
    
    cgfi_time[r] <-
      as.numeric(
        Sys.time()-start,
        units="secs"
      )
    
    ci <- quantile(
      draws,
      c(.025,.975)
    )
    
    coverage[r,3] <-
      theta0 >= ci[1] &
      theta0 <= ci[2]
    
    lengths[r,3] <- diff(ci)
    
    if(r %% 10 == 0){
      cat(
        error_type,
        ": completed",
        r,
        "of",
        R,
        "\n"
      )
    }
    
  }
  
  data.frame(
    
    Error=error_type,
    
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
    
    MCSE=c(
      sqrt(mean(coverage[,1])*(1-mean(coverage[,1]))/R),
      sqrt(mean(coverage[,2])*(1-mean(coverage[,2]))/R),
      sqrt(mean(coverage[,3])*(1-mean(coverage[,3]))/R)
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
# Run Experiment 4
############################################################

error_results <- do.call(
  rbind,
  lapply(
    c(
      "Normal",
      "t5",
      "Heavy"
    ),
    simulate_errors
  )
)

print(error_results)

############################################################
# Coverage Plot
############################################################

ggplot(
  error_results,
  aes(
    x=Error,
    y=Coverage,
    group=Method,
    color=Method
  )
)+
  geom_point(size=3)+
  geom_line(linewidth=1)+
  geom_hline(
    yintercept=.95,
    linetype="dashed"
  )+
  ylim(.5,1)+
  theme_minimal()+
  labs(
    y="Coverage Probability",
    x="Error Distribution"
  )


############################################################
############################################################
# Experiment 5:
# Heteroskedastic Error Robustness
############################################################
############################################################


simulate_heteroskedastic <- function(error_type,
                                     q=100,
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
    
    
    ####################################################
    # Generate regressors and instruments
    ####################################################
    
    
    X <- rnorm(n)
    
    
    Z <- matrix(
      rnorm(n*q),
      nrow=n,
      ncol=q
    )
    
    
    ####################################################
    # Heteroskedastic disturbances
    ####################################################
    
    
    if(error_type=="Homoskedastic"){
      
      sigma <- rep(1,n)
      
    }
    
    
    if(error_type=="Moderate"){
      
      sigma <- 1 + 0.5*abs(X)
      
    }
    
    
    if(error_type=="Strong"){
      
      sigma <- exp(0.5*X)
      
    }
    
    
    u <- sigma * rt(
      n,
      df=5
    )
    
    
    
    ####################################################
    # Outcome model
    ####################################################
    
    
    Y <- theta0*X +
      0.5*X^2 +
      u
    
    
    ####################################################
    # Two-step GMM
    ####################################################
    
    
    prep <- prepare_moments(
      Y,
      X,
      Z
    )
    
    
    theta_initial <- gmm_est(
      prep
    )
    
    
    Omega <- estimate_Omega(
      theta_initial,
      Y,
      X,
      Z
    )
    
    
    W <- solve(Omega)
    
    
    theta_hat <- gmm_est(
      prep,
      W
    )
    
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    
    W <- solve(Omega)
    
    
    copula <- prepare_copula(
      Omega
    )
    
    
    
    ####################################################
    # Wald
    ####################################################
    
    
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
    
    
    
    ####################################################
    # Bootstrap
    ####################################################
    
    
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
      
      
      theta_b0 <- gmm_est(
        prep_b
      )
      
      
      Omega_b <- estimate_Omega(
        theta_b0,
        Y[id],
        X[id],
        Z[id,]
      )
      
      
      W_b <- solve(Omega_b)
      
      
      boot[b] <- gmm_est(
        prep_b,
        W_b
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
    
    
    
    ####################################################
    # CGFI
    ####################################################
    
    
    start <- Sys.time()
    
    
    cgfi <- numeric(B)
    
    
    for(b in 1:B){
      
      cgfi[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula,
        W
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
        error_type,
        ": completed",
        r,
        "of",
        R,
        "\n"
      )
      
    }
    
  }
  
  
  data.frame(
    
    Error_Structure=error_type,
    
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
    
    MCSE=c(
      sqrt(mean(coverage[,1])*
             (1-mean(coverage[,1]))/R),
      sqrt(mean(coverage[,2])*
             (1-mean(coverage[,2]))/R),
      sqrt(mean(coverage[,3])*
             (1-mean(coverage[,3]))/R)
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
# Run heteroskedasticity experiment
############################################################


hetero_results <- do.call(
  rbind,
  lapply(
    c(
      "Homoskedastic",
      "Moderate",
      "Strong"
    ),
    function(e){
      
      simulate_heteroskedastic(
        error_type=e,
        q=100,
        n=500,
        R=100,
        B=200
      )
      
    }
  )
)


print(hetero_results)



############################################################
# Plot coverage
############################################################


ggplot(
  hetero_results,
  aes(
    x=Error_Structure,
    y=Coverage,
    color=Method,
    group=Method
  )
)+
  geom_point(
    size=3
  )+
  geom_line(
    linewidth=1
  )+
  geom_hline(
    yintercept=.95,
    linetype="dashed"
  )+
  ylim(.5,1)+
  theme_minimal()+
  labs(
    x="Error Structure",
    y="Coverage Probability"
  )

############################################################
############################################################
# Experiment 6:
# Weak Instrument / Weak Identification Robustness
############################################################
############################################################


simulate_weak_iv <- function(pi_strength,
                             q=100,
                             R=100,
                             B=200,
                             n=500,
                             theta0=2){
  
  
  coverage <- matrix(0,R,3)
  lengths <- matrix(0,R,3)
  
  boot_time <- numeric(R)
  cgfi_time <- numeric(R)
  
  
  for(r in 1:R){
    
    
    ####################################################
    # Instruments
    ####################################################
    
    Z <- matrix(
      rnorm(n*q),
      nrow=n,
      ncol=q
    )
    
    
    ####################################################
    # Endogenous regressor
    ####################################################
    
    v <- rnorm(n)
    
    X <- pi_strength*rowMeans(Z) + v
    
    
    ####################################################
    # Structural error
    ####################################################
    
    epsilon <- rnorm(n)
    
    # creates endogeneity:
    # Cov(X,u) != 0
    u <- 0.7*v + epsilon
    
    
    ####################################################
    # Original nonlinear structural model
    ####################################################
    
    Y <- theta0*X +
      0.5*X^2 +
      u
    
    
    ####################################################
    # GMM estimation (same as baseline)
    ####################################################
    
    prep <- prepare_moments(
      Y,
      X,
      Z
    )
    
    
    theta_initial <- gmm_est(
      prep
    )
    
    
    Omega <- estimate_Omega(
      theta_initial,
      Y,
      X,
      Z
    )
    
    
    W <- solve(Omega)
    
    
    theta_hat <- gmm_est(
      prep,
      W
    )
    
    
    Omega <- estimate_Omega(
      theta_hat,
      Y,
      X,
      Z
    )
    
    
    W <- solve(Omega)
    
    
    copula <- prepare_copula(
      Omega
    )
    
    
    
    ####################################################
    # Wald
    ####################################################
    
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
    
    
    
    ####################################################
    # Bootstrap
    ####################################################
    
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
      
      
      theta_b <- gmm_est(
        prep_b
      )
      
      
      boot[b] <- theta_b
      
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
    
    
    
    ####################################################
    # CGFI
    ####################################################
    
    start <- Sys.time()
    
    
    cgfi <- numeric(B)
    
    
    for(b in 1:B){
      
      cgfi[b] <- cgfi_draw(
        theta_hat,
        prep,
        copula,
        W
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
        "\n"
      )
      
    }
    
    
  }
  
  
  cover <- c(
    mean(coverage[,1]),
    mean(coverage[,2]),
    mean(coverage[,3])
  )
  
  
  data.frame(
    
    Instrument_Strength=pi_strength,
    
    Method=c(
      "Wald",
      "Bootstrap",
      "CGFI"
    ),
    
    Coverage=cover,
    
    MCSE=sqrt(
      cover*(1-cover)/R
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
# Run weak-IV experiment
############################################################


weak_iv_results <- do.call(
  rbind,
  lapply(
    c(
      2,
      1,
      0.5,
      0.25,
      0.1
    ),
    function(pi_strength){
      
      simulate_weak_iv(
        pi_strength=pi_strength,
        q=100,
        n=500,
        R=100,
        B=200
      )
      
    }
  )
)


print(weak_iv_results)



############################################################
# Plot coverage
############################################################


ggplot(
  weak_iv_results,
  aes(
    x=Instrument_Strength,
    y=Coverage,
    color=Method,
    group=Method
  )
)+
  geom_point(
    size=3
  )+
  geom_line(
    linewidth=1
  )+
  geom_hline(
    yintercept=.95,
    linetype="dashed"
  )+
  ylim(.5,1)+
  theme_minimal()+
  labs(
    x="Instrument Strength ($\\pi$)",
    y="Coverage Probability"
  )