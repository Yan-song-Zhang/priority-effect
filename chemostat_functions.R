####nonlinear

simulate_chemostat_nonlinear <- function(params, state, times, mode = "substitu") {
  dt <- times[2] - times[1]
  
  D <- params$D
  R1_in <- params$R1_in
  R2_in <- params$R2_in
  c11 <- params$c11; c12 <- params$c12
  w11 <- params$w11; w12 <- params$w12
  K11 <- params$K11; K12 <- params$K12  # Half-saturation for N1
  c21 <- params$c21; c22 <- params$c22
  w21 <- params$w21; w22 <- params$w22
  K21 <- params$K21; K22 <- params$K22  # Half-saturation for N2
  sigma <- params$sigma
  ex_threshold <- 0.01
  
  # Initialize state
  n <- length(times)
  R1 <- numeric(n); R2 <- numeric(n)
  N1 <- numeric(n); N2 <- numeric(n)
  
  R1[1] <- state["R1"]; R2[1] <- state["R2"]
  N1[1] <- state["N1"]; N2[1] <- state["N2"]
  
  # Track extinction
  N1_extinct <- N1[1] < ex_threshold
  N2_extinct <- N2[1] < ex_threshold
  
  for (i in 2:n) {
    r1 <- R1[i-1]; r2 <- R2[i-1]
    n1 <- N1[i-1]; n2 <- N2[i-1]
    
    # Update extinction status
    if (n1 < ex_threshold) N1_extinct <- TRUE
    if (n2 < ex_threshold) N2_extinct <- TRUE
    
    # Resource dynamics WITH MONOD KINETICS
    dR1_dt <- D * (R1_in - r1)
    dR2_dt <- D * (R2_in - r2)
    
    if (!N1_extinct) {
      uptake_R1_N1 <- (c11 * r1) / (K11 + r1) * n1
      uptake_R2_N1 <- (c12 * r2) / (K12 + r2) * n1
      dR1_dt <- dR1_dt - uptake_R1_N1
      dR2_dt <- dR2_dt - uptake_R2_N1
    }
    if (!N2_extinct) {
      uptake_R1_N2 <- (c21 * r1) / (K21 + r1) * n2
      uptake_R2_N2 <- (c22 * r2) / (K22 + r2) * n2
      dR1_dt <- dR1_dt - uptake_R1_N2
      dR2_dt <- dR2_dt - uptake_R2_N2
    }
    
    # Population dynamics WITH MONOD KINETICS
    noise_N1 <- if (N1_extinct) 0 else sigma * n1 * rnorm(1, 0, sqrt(dt))
    noise_N2 <- if (N2_extinct) 0 else sigma * n2 * rnorm(1, 0, sqrt(dt))
    
    if (N1_extinct) {
      dN1_dt <- 0
    } else if (mode == "substitu") {
      # Substitutive: Weighted sum of resource contributions
      dN1_dt <- (w11 * (c11 * r1)/(K11 + r1) + 
                   w12 * (c12 * r2)/(K12 + r2) - D) * n1
    } else {
      # Essential: Growth limited by most scarce resource
      growth_R1 <- w11 * (c11 * r1)/(K11 + r1)
      growth_R2 <- w12 * (c12 * r2)/(K12 + r2)
      dN1_dt <- (min(growth_R1, growth_R2) - D) * n1
    }
    
    # Same for N2
    if (N2_extinct) {
      dN2_dt <- 0
    } else if (mode == "substitu") {
      dN2_dt <- (w21 * (c21 * r1)/(K21 + r1) + 
                   w22 * (c22 * r2)/(K22 + r2) - D) * n2
    } else {
      growth_R1_N2 <- w21 * (c21 * r1)/(K21 + r1)
      growth_R2_N2 <- w22 * (c22 * r2)/(K22 + r2)
      dN2_dt <- (min(growth_R1_N2, growth_R2_N2) - D) * n2
    }
    
    # Euler update with noise
    R1[i] <- max(r1 + dR1_dt * dt, 0)
    R2[i] <- max(r2 + dR2_dt * dt, 0)
    N1[i] <- if (N1_extinct) 0 else max(n1 + (dN1_dt+noise_N1) * dt , 0)
    N2[i] <- if (N2_extinct) 0 else max(n2 + (dN2_dt++noise_N2) * dt, 0)
  }
  
  data.frame(time = times, R1 = R1, R2 = R2, N1 = N1, N2 = N2)
}






#############
##linear####
############
simulate_chemostat <- function(params, state, times, mode = "substitu") {
  dt <- times[2] - times[1]
  
  # Extract parameters
  D <- params$D
  R1_in <- params$R1_in
  R2_in <- params$R2_in
  c11 <- params$c11; c12 <- params$c12
  w11 <- params$w11; w12 <- params$w12
  #m1 <- params$m1  # Mortality for species 1
  c21 <- params$c21; c22 <- params$c22
  w21 <- params$w21; w22 <- params$w22
  #m2 <- params$m2  # Mortality for species 2
  sigma <- params$sigma
  ex_threshold <- 0.001
  
  
  # Initialize state
  n <- length(times)
  R1 <- numeric(n); R2 <- numeric(n)
  N1 <- numeric(n); N2 <- numeric(n)
  
  R1[1] <- state["R1"]; R2[1] <- state["R2"]
  N1[1] <- state["N1"]; N2[1] <- state["N2"]
  
  # Track extinction
  N1_extinct <- N1[1] < ex_threshold
  N2_extinct <- N2[1] < ex_threshold
  
  for (i in 2:n) {
    # Current state
    r1 <- R1[i-1]; r2 <- R2[i-1]
    n1 <- N1[i-1]; n2 <- N2[i-1]
    
    # Update extinction status
    if (n1 < ex_threshold) N1_extinct <- TRUE
    if (n2 < ex_threshold) N2_extinct <- TRUE
    
    # Resource dynamics
    dR1_dt <- D * (R1_in - r1)
    dR2_dt <- D * (R2_in - r2)
    
    if (!N1_extinct) {
      dR1_dt <- dR1_dt - c11 * n1 * r1
      dR2_dt <- dR2_dt - c12 * n1 * r2
    }
    if (!N2_extinct) {
      dR1_dt <- dR1_dt - c21 * n2 * r1
      dR2_dt <- dR2_dt - c22 * n2 * r2
    }
    
    # Population dynamics
    noise_N1 <- if (N1_extinct) 0 else sigma * n1 * rnorm(1, mean = 0, sd = sqrt(dt))
    noise_N2 <- if (N2_extinct) 0 else sigma * n2 * rnorm(1, mean = 0, sd = sqrt(dt))
    
    
    if (N1_extinct) {
      dN1_dt <- 0
    } else if (mode == "substitu") {
      dN1_dt <- (c11 * r1 * w11 + c12 * r2 * w12 - D ) * n1 + noise_N1
    } else {
      dN1_dt <- (min(c11 * r1 * w11, c12 * r2 * w12) - D ) * n1 + noise_N1
    }
    
    if (N2_extinct) {
      dN2_dt <- 0
    } else if (mode == "substitu") {
      dN2_dt <- (c21 * r1 * w21 + c22 * r2 * w22 - D ) * n2 + noise_N2
    } else {
      dN2_dt <- (min(c21 * r1 * w21, c22 * r2 * w22) - D ) * n2 + noise_N2
    }
    
    # Euler update
    R1[i] <- max(r1 + dR1_dt * dt, 0)
    R2[i] <- max(r2 + dR2_dt * dt, 0)
    N1[i] <- if (N1_extinct) 0 else max(n1 + dN1_dt * dt, 0)
    N2[i] <- if (N2_extinct) 0 else max(n2 + dN2_dt * dt, 0)
  }
  
  data.frame(time = times, R1 = R1, R2 = R2, N1 = N1, N2 = N2)
}









#######
###N1&N2
#########

simulate_chemostat_manual_N1N2 <- function(params, state, times, mode = "substitu") {
  dt <- times[2] - times[1]  # timestep
  
  D <- params$D
  R1_in <- params$R1_in
  R2_in <- params$R2_in
  c11 <- params$c11
  c12 <- params$c12
  w11 <- params$w11
  w12 <- params$w12
  c21 <- params$c21
  c22 <- params$c22
  w21 <- params$w21
  w22 <- params$w22
  sigma <- 0.1
  R_1 <- numeric(length(times))
  R_2 <- numeric(length(times))
  N_1 <- numeric(length(times))
  N_2 <- numeric(length(times))
  
  
  R_1[1] <- state["R1"]
  R_2[1] <- state["R2"]
  N_1[1] <- state["N1"]
  N_2[1] <- state["N2"]
  
  N1_extinct <- FALSE
  N2_extinct <- FALSE
  
  
  for (i in 2:length(times)) {
    R1 = R_1[i-1] 
    R2 = R_2[i-1]
    N1 = N_1[i-1]
    N2 = N_2[i-1]
    
    
    if (N1 <= 0.001) N1_extinct <- TRUE
    if (N2 <= 0.001) N2_extinct <- TRUE
    
    # Resource dynamics
    
    
    if (!N1_extinct & !N2_extinct) {
      dR1_dt <- D * (R1_in - R1) - c11 * N1 * R1 - c21 * N2 * R1
      dR2_dt <- D * (R2_in - R2) - c12 * N1 * R2 - c22 * N2 * R2
      noise_N1 <- sigma * N1 * rnorm(1, mean = 0, sd = sqrt(dt))
      noise_N2 <- sigma * N2 * rnorm(1, mean = 0, sd = sqrt(dt))
      if (mode == "substitu") {
        dN1_dt <- (c11 * R1 * w11 + c12 * R2 * w12 - D) * N1 +  noise_N1
        dN2_dt <- (c21 * R1 * w21 + c22 * R2 * w22 - D) * N2 +  noise_N2
        
      } else if (mode == "essen") {
        dN1_dt <- (min(c11 * R1 * w11, c12 * R2 * w12) - D) * N1 + noise_N1
        dN2_dt <- (min(c21 * R1 * w21, c22 * R2 * w22) - D) * N2 + noise_N2
        
      } else {
        stop("Invalid mode!")
      }
      
    }
    
    
    if (N1_extinct & !N2_extinct) {
      dR1_dt <- D * (R1_in - R1)  - c21 * N2 * R1
      dR2_dt <- D * (R2_in - R2)  - c22 * N2 * R2
      noise_N2 <- sigma * N2 * rnorm(1, mean = 0, sd = sqrt(dt))
      if (mode == "substitu") {
        dN1_dt <- 0
        dN2_dt <- (c21 * R1 * w21 + c22 * R2 * w22 - D) * N2 +  noise_N2
        
      } else if (mode == "essen") {
        dN1_dt <- 0
        dN2_dt <- (min(c21 * R1 * w21, c22 * R2 * w22) - D) * N2 + noise_N2
        
      } else {
        stop("Invalid mode!")
      }
      
    }
    
    
    if (!N1_extinct & N2_extinct) {
      dR1_dt <- D * (R1_in - R1) - c11 * N1 * R1 
      dR2_dt <- D * (R2_in - R2) - c12 * N1 * R2 
      noise_N1 <- sigma * N1 * rnorm(1, mean = 0, sd = sqrt(dt))
      if (mode == "substitu") {
        dN1_dt <- (c11 * R1 * w11 + c12 * R2 * w12 - D) * N1 +  noise_N1
        dN2_dt <- 0
        
      } else if (mode == "essen") {
        dN1_dt <- (min(c11 * R1 * w11, c12 * R2 * w12) - D) * N1 + noise_N1
        dN2_dt <- 0
        
      } else {
        stop("Invalid mode!")
      }
      
    }
    
    
    if (N1_extinct & N2_extinct) {
      dR1_dt <- D * (R1_in - R1) 
      dR2_dt <- D * (R2_in - R2)
      
      if (mode == "substitu") {
        dN1_dt <- 0
        dN2_dt <- 0
        
      } else if (mode == "essen") {
        dN1_dt <- 0
        dN2_dt <- 0
        
      } else {
        stop("Invalid mode!")
      }
      
    }
    # Euler update
    
    
    
    R_1[i] <- max(R1 + dR1_dt * dt, 0)
    R_2[i] <- max(R2 + dR2_dt * dt, 0)
    N_1[i] <- if (N1_extinct) 0 else max(N1 + dN1_dt * dt, 0)
    N_2[i] <- if (N2_extinct) 0 else max(N2 + dN2_dt * dt, 0)
    
    
  }
  
  data.frame(time = times, R1 = R_1, R2 = R_2, N1 = N_1, N2 = N_2)
}
