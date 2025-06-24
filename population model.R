#logistic growth model
setwd("D:/zhang lab/experiment priority effect/Rcode")
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyverse)

load('2025.6.11.RData')

logistic_growth <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dn_dt <- r * (1 - n / K) * n
    return(list(c(dn_dt)))
  })
}

parameters <- c(r = 0.05,K = 100) 
state <- c(n = 20)  
times <- seq(0, 200, by = 10) 
out <- ode(y = state, times = times, func = logistic_growth, parms = parameters)
out_df <- as.data.frame(out)

# 绘制结果
plot(out_df$time, out_df$n, type = "b", col = "blue", xlab = "Time", ylab = "种群大小",
     main = "Logistic Growth Model")


#数值积分
r = 0.05
K = 100
N0 = 20
dt = 10
tmax = 200

time = seq(0, tmax, by = dt)
N <- numeric(length(time))
N[1] <- N0

for (i in 2: length(time)){
  dN = r*N[i-1]*(1-N[i-1]/K)
  N[i] = N[i-1]+dN*dt
}
plot(time, N, type = 'b', col = 'blue', lwd = 2, pch = 19,xlab="Time", ylab="Population Size",
     main="Discrete Logistic Growth (Ricker Model)" )


#contrast of the two methods
layout(matrix(c(1, 2), 2, 1))

plot(time, N, type = 'b', col = 'blue', lwd = 2, pch = 19,xlab="Time", ylab="Population Size",
     main="Discrete Logistic Growth (Ricker Model)" )
plot(time, out_df$n, type = "b", col = "blue",lwd = 2, pch = 19, xlab = "Time", ylab = "种群大小",
     main = "Logistic Growth Model(deSolve)")
layout(1)



##对数变化率 离散时间

r = 0.05
K = 100
N0 = 20
times =seq(1,200,1)

N = numeric(length(times)+1)
N[1]=N0
for (i in 2:(length(times)+1)){ 
  N[i] = N[i-1]*exp(r * (1 - N[i-1] / K))
}
times = c(c(0),times)
plot(times,N,type = 'l',main = "Discrete Time Logistic Growth Model", xlab = "Time", ylab = "Population Size" )


##离散时间 deSolve
parameters = c(r = 0.05,K = 100)
state = c(N = 20)
t = c(2,4,7,11,15,23,29,35,41,68,90,100,150,200)

discrete_logistic = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    dn.dt = r *(1-N/K)*N
    return(list(c(dn.dt)))
  }
       ) 
}
out = ode(y = state, times = t, func = discrete_logistic,parms = parameters)
out_df = as.data.frame(out)


# 绘制结果
plot(out_df$time, out_df$N, type = "b", col = "blue", xlab = "Time", ylab = "种群大小",
     main = "Logistic Growth Model")


# Monod equation

monod_model = function(t, state, parameters){
  N = state[1]
  R = state[2]
  
  r_max = parameters_m['r_max']
  Ks = parameters_m['Ks']
  c = parameters_m['c']
  
  dN_dt = r_max * N * R /(Ks + R)
  dR_dt = -c *  dn_dt
  if (R + dR_dt < 0) {
    dR_dt = 0  
  }
  
  return(list(c(dN_dt,dR_dt)))
}

state_m = c(N0 = 10, R0 = 5)
parameters_m = c(r_max = 0.05, Ks = 10, c = 0.1)
times_m = seq(0,80,by=1)

ode_m = ode(y= state_m, times = times_m,func = monod_model,parms = parameters_m ) 

matplot(ode_m[, "time"], ode_m[, -1], type = "l", lty = 1, col = c("blue", "red"), 
        xlab = "Time", ylab = "Value", main = "Monod Model (Population & Resource)")
legend("topleft", legend = c("Population (N)", "Resource (R)"), col = c("blue", "red"), lty = 1)



# Monod equation 1 species 1 resources chemostat

chemostat_model = function(t,state,parameters){
  with(as.list(c(state,parameters)), {
    
    r = r_max*R/(R+Ks)
    dR_dt = D*(R_in - R) - r*N/Y
    dN_dt = (r-D)*N
    
    list(c(dR_dt,dN_dt))
  }
       )
  
}
parameters_chemo = c(r_max = 0.6,Ks = 0.1, D = 0.5, Y = 0.4,  R_in = 10)

state_chemo = c(R = 10, N = 0.1)

times_chemo = seq(0,50,by = 0.1)

ode_chemo = ode (y = state_chemo, times = times_chemo, func = chemostat_model, parms = parameters_chemo)

out_df = as.data.frame(ode_chemo)

library(ggplot2)
ggplot(out_df, aes(x = time)) +
  geom_line(aes(y = R, color = "Resource"),linewidth = 1) +
  geom_line(aes(y = N, color = "Biomass"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  theme_minimal()



# Monod equation 1 species 2 resources chemostat substitutable

chemostat_model_2 = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    r1 = r1_max * R1 / (Ks1+R1)
    r2 = r2_max * R2 / (Ks2+R2)
    dR1_dt = D*(R1_in - R1) - r1*N/Y1
    dR2_dt = D*(R2_in - R2) - r2*N/Y2
    dN_dt = (r1+r2-D)*N
    list(c(dR1_dt,dR2_dt,dN_dt))
  })
}

parameters_chemo_2 = c(r1_max = 0.5,Ks1 = 0.1,r2_max = 0.8,Ks2 = 0.4, D = 0.1, Y1 = 0.4,Y2 = 0.6,  
                       R1_in = 10,R2_in = 20)

state_chemo_2 = c(R1 = 10,R2 = 20, N = 0.1)

times_chemo_2 = seq(0,50,by = 0.1)

ode_chemo_2 = ode (y = state_chemo_2, times = times_chemo_2, func = chemostat_model_2, parms = parameters_chemo_2)

out_df_2 = as.data.frame(ode_chemo_2)

library(ggplot2)
ggplot(out_df_2, aes(x = time)) +
  geom_line(aes(y = R1, color = "Resource_1"),linewidth = 1) +
  geom_line(aes(y = R2, color = "Resource_2"),linewidth = 1) +
  geom_line(aes(y = N, color = "Biomass"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation-substitute",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  scale_y_log10() +
  theme_minimal()


# Monod equation 1 species 2 resources chemostat essential

chemostat_model_3 = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    r1 = r1_max * R1 / (Ks1+R1)
    r2 = r2_max * R2 / (Ks2+R2)
    dR1_dt = D*(R1_in - R1) - r1*N/Y1
    dR2_dt = D*(R2_in - R2) - r2*N/Y2
    dN_dt = (min(r1,r2)-D)*N
    list(c(dR1_dt,dR2_dt,dN_dt))
  })
}


parameters_chemo_3 = c(r1_max = 0.5,Ks1 = 0.1,r2_max = 0.8,Ks2 = 0.4, D = 0.1, Y1 = 0.4,Y2 = 0.6,  
                       R1_in = 10,R2_in = 20)

state_chemo_3 = c(R1 = 10,R2 = 20, N = 0.1)

times_chemo_3 = seq(0,50,by = 0.1)

ode_chemo_3 = ode (y = state_chemo_3, times = times_chemo_3, func = chemostat_model_3, parms = parameters_chemo_3)

out_df_3 = as.data.frame(ode_chemo_3)

library(ggplot2)
ggplot(out_df_3, aes(x = time)) +
  geom_line(aes(y = R1, color = "Resource_1"),linewidth = 1) +
  geom_line(aes(y = R2, color = "Resource_2"),linewidth = 1) +
  geom_line(aes(y = N, color = "Biomass"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation-essential",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  theme_minimal()


###########
###########
###########
## Monod equation 1 species 2 resources chemostat 1 function
chemostat_model_1 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    r1 = r1_max * R1# / (Ks1+R1)
    r2 = r2_max * R2#/ (Ks2+R2)
    dR1_dt = D*(R1_in - R1) - r1*N/Y1
    dR2_dt = D*(R2_in - R2) - r2*N/Y2
   if (mode == "substitu") {
      dN_dt = (r1 + r2 - D) * N  
    } else if (mode == "essen") {
     dN_dt = (min(r1, r2) - D) * N  
   } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
  }
    list(c(dR1_dt,dR2_dt,dN_dt))
  })
}

parameters_chemo = c(r1_max = 0.5,Ks1 = 0.1,r2_max = 0.1,Ks2 = 0.4, D = 0.1, Y1 = 0.4,Y2 = 0.4,
                     R1_in = 10,R2_in = 10)

state_chemo = c(R1 = 10,R2 = 10, N = 0.1)

times_chemo = seq(0,50,by = 0.1)

ode_chemo = ode (y = state_chemo, times = times_chemo, func = chemostat_model, parms = parameters_chemo,
                 mode='essen')

out_df_chemo = as.data.frame(ode_chemo)

ggplot(out_df_chemo, aes(x = time)) +
  geom_line(aes(y = R1, color = "Resource_1"),linewidth = 1) +
  geom_line(aes(y = R2, color = "Resource_2"),linewidth = 1) +
  geom_line(aes(y = N, color = "Biomass"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation-substitu",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  scale_y_log10() +
  theme_minimal()




#2 species 1 resource

chemostat_model_2s1r = function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    r1 = r1_max * R #/ (Ks1+R)
    r2 = r2_max * R #/ (Ks2+R)
    dR_dt = D*(R_in - R) - r1*N1/Y1- r2*N2/Y2
    dN1_dt = (r1  - D) * N1  
    dN2_dt = (r2  - D) * N2 
    list(c(dR_dt,dN1_dt,dN2_dt))
  })
}

parameters_chemo = c(r1_max = 0.1,Ks1 = 0.1,r2_max = 0.15,Ks2 = 0.4, D = 0.1, Y1 = 0.4,Y2 = 0.6,
                     R_in = 10)

state_chemo = c(R = 10, N1 = 1,N2 = 0.1)

times_chemo = seq(0,50,by = 0.1)

ode_chemo = ode (y = state_chemo, times = times_chemo, func = chemostat_model_2s1r, parms = parameters_chemo)

out_df_chemo = as.data.frame(ode_chemo)

ggplot(out_df_chemo, aes(x = time)) +
  geom_line(aes(y = N1, color = "Species_1"),linewidth = 1) +
  geom_line(aes(y = N2, color = "Species_2"),linewidth = 1) +
  geom_line(aes(y = R, color = "Resource"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation 2 species 1 resource",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  scale_y_log10() +
  theme_minimal()

##

#2 species 2 resource

chemostat_model_2s2r = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    r11 = r11_max * R1 #/ (Ks1+R)
    r12 = r12_max * R2 #/ (Ks2+R).
    r21 = r21_max * R1 #/ (Ks1+R)
    r22 = r22_max * R2 #/ (Ks2+R)
    dR1_dt = D*(R1_in - R1) - r11*N1/Y11- r21*N2/Y21
    dR2_dt = D*(R2_in - R2) - r12*N1/Y12- r22*N2/Y22
    if (mode == "substitu") {
      dN1_dt = (r11 + r12 - D) * N1 
      dN2_dt = (r21 + r22 - D) * N2 

    } else if (mode == "essen") {
      dN1_dt = (min(r11, r12) - D) * N1 
      dN2_dt = (min(r21, r22) - D) * N2 
      } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
  
    list(c(dR1_dt,dR2_dt,dN1_dt,dN2_dt))
  })
}

param_grid <- expand.grid(
  r11_max = seq(0.1, 0.9, by = 0.1),
  r12_max = seq(0.1, 0.9, by = 0.1),
  r21_max = seq(0.1, 0.9, by = 0.1),
  r22_max = seq(0.1, 0.9, by = 0.1)
)

##param_grid = param_grid [((param_grid[,4]>param_grid[,2])&(param_grid[,2]>param_grid[,1])&(param_grid[,1]>param_grid[,3])),]

initial_conditions <- list(
  list(name = "0", state_chemo = c(R1 = 10, R2 = 10, N1 = 0.1, N2 = 1)),
  list(name = "1", state_chemo = c(R1 = 10, R2 = 10, N1 = 1, N2 = 0.1))
)


fixed_params = c(D = 0.1, Y11 = 0.8,Y12 = 0.8,Y21 = 0.8,Y22 = 0.8,R1_in = 10,R2_in = 10)

times_chemo = seq(0,50,by = 0.1)

results <- lapply(1:nrow(param_grid), function(i) {
  params <- c( param_grid[i, ],fixed_params)
  
  lapply(initial_conditions, function(init) {
    sim <- ode(y = init$state_chemo, times = times_chemo, func = chemostat_model_2s2r, parms = params, mode = "essen")
    sim_df <- as.data.frame(sim)
    sim_df$param_set <- i  # Track parameter set
    sim_df$initial_condition <- init$name  # Track initial condition
    sim_df
    
  })
})

results_df <- bind_rows(unlist(results, recursive = FALSE))


for (i in 1:126) {
  for(j in 0:1){
    plot.df <- results_df[(results_df$param_set==i) & (results_df$initial_condition==j),]
    
    plot <- ggplot(plot.df, aes(x = time)) +
      geom_line(aes(y = N1, color = "Species_1"), linewidth = 1) +
      geom_line(aes(y = N2, color = "Species_2"), linewidth = 1) +
      #geom_line(aes(y = R1, color = "Resource_1"), linewidth = 1) +
      #geom_line(aes(y = R2, color = "Resource_2"), linewidth = 1) +
      labs(title = paste("Chemostat Simulation:", i,j),
           x = "Time",
           y = "Concentration",
           color = "Legend") +
      scale_y_log10() +
      theme_minimal()
    setwd("D:/zhang lab/simulate2")
    ggsave(filename = paste0("chemostat_sim_", i, "_", j, ".png"), plot = plot, width = 8, height = 6)
  }
}

##
##
##
##
###substitute resources in different Y
set.seed(123)
n_samples <- 1000
param_df <- data.frame(
  r11_max = numeric(n_samples),
  r12_max = numeric(n_samples),
  r21_max = numeric(n_samples),
  r22_max = numeric(n_samples),
  Y11 = numeric(n_samples),
  Y12 = numeric(n_samples),
  Y21 = numeric(n_samples),
  Y22 = numeric(n_samples)
)

for (i in 1:n_samples) {
  sample_vals <- sample(1:100, 8, replace = TRUE) / 100  # Normalize to [0,1]
  param_df[i, ] <- sample_vals
}


initial_conditions <- list(
  list(name = "0", state_chemo = c(R1 = 10, R2 = 10, N1 = 0.1, N2 = 1)),
  list(name = "1", state_chemo = c(R1 = 10, R2 = 10, N1 = 1, N2 = 0.1))
)

fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,100,by = 0.1)

results.Y <- lapply(1:nrow(param_df), function(i) {
  params <- c(param_df[i, ],fixed_params)
  
  lapply(initial_conditions, function(init) {
    sim <- ode(y = init$state_chemo, times = times_chemo, func = chemostat_model_2s2r, parms = params, mode = "substitu")
    sim_df <- as.data.frame(sim)
    sim_df$param_set <- i  # Track parameter set
    sim_df$initial_condition <- init$name  # Track initial condition
    sim_df
    
  })
})

results_df.Y <- bind_rows(unlist(results.Y, recursive = FALSE))

priority.results = data.frame()
priority.results[1,1] =NA

for (i in 1:1000) {
  
   df.Y.0 <- results_df.Y[(results_df.Y$param_set==i) & (results_df.Y$initial_condition==0),]
   df.Y.1 <- results_df.Y[(results_df.Y$param_set==i) & (results_df.Y$initial_condition==1),]
   if  (((df.Y.0[980,4]-df.Y.0[970,4])*(df.Y.1[980,4]-df.Y.1[970,4]))<0){
     priority.results[i,1]=i
   }
   if  (((df.Y.0[980,5]-df.Y.0[970,5])*(df.Y.1[980,5]-df.Y.1[970,5]))<0){
     priority.results[i,2]=i
   }
   if  (((df.Y.0[980,5]-df.Y.0[970,5])*(df.Y.1[980,5]-df.Y.1[970,5]))==0){
     priority.results[i,3]=i
   }
   if  (((df.Y.0[980,5]-df.Y.0[970,5])*(df.Y.1[980,5]-df.Y.1[970,5]))==0){
     priority.results[i,4]=i
   }
}

priority.results <- priority.results[rowSums(is.na(priority.results)) != ncol(priority.results), ]
priority.results.final <- na.omit(priority.results)



###substitute resources with nonlinear growth rate

chemostat_model_2s2r_nonlinear = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    r11 = r11_max * R1 / (Ks11+R1)
    r12 = r12_max * R2 / (Ks12+R2)
    r21 = r21_max * R1 / (Ks21+R1)
    r22 = r22_max * R2 / (Ks22+R2)
    dR1_dt = D*(R1_in - R1) - r11*N1/Y11- r21*N2/Y21
    dR2_dt = D*(R2_in - R2) - r12*N1/Y12- r22*N2/Y22
    if (mode == "substitu") {
      dN1_dt = (r11 + r12 - D) * N1 
      dN2_dt = (r21 + r22 - D) * N2 
      
    } else if (mode == "essen") {
      dN1_dt = (min(r11, r12) - D) * N1 
      dN2_dt = (min(r21, r22) - D) * N2 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt,dN2_dt))
  })
}
set.seed(123)
n_samples <- 1000
param_df_nonlinear <- data.frame(
  r11_max = numeric(n_samples),
  r12_max = numeric(n_samples),
  r21_max = numeric(n_samples),
  r22_max = numeric(n_samples),
  Y11 = numeric(n_samples),
  Y12 = numeric(n_samples),
  Y21 = numeric(n_samples),
  Y22 = numeric(n_samples),
  Ks11 = numeric(n_samples),
  Ks12 = numeric(n_samples),
  Ks21 = numeric(n_samples),
  Ks22 = numeric(n_samples)
)

for (i in 1:n_samples) {
  sample_vals <- sample(1:100, 12, replace = TRUE) / 100  # Normalize to [0,1]
  param_df_nonlinear[i, ] <- sample_vals
}


initial_conditions <- list(
  list(name = "0", state_chemo = c(R1 = 10, R2 = 10, N1 = 0.1, N2 = 1)),
  list(name = "1", state_chemo = c(R1 = 10, R2 = 10, N1 = 1, N2 = 0.1))
)

fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,100,by = 0.1)

results.nonlinear <- lapply(1:nrow(param_df_nonlinear), function(i) {
  params <- c(param_df_nonlinear[i, ],fixed_params)
  
  lapply(initial_conditions, function(init) {
    sim <- ode(y = init$state_chemo, times = times_chemo, func = chemostat_model_2s2r_nonlinear, parms = params, mode = "substitu")
    sim_df <- as.data.frame(sim)
    sim_df$param_set <- i  # Track parameter set
    sim_df$initial_condition <- init$name  # Track initial condition
    sim_df
    
  })
})

results_df.nonlinear <- bind_rows(unlist(results.nonlinear, recursive = FALSE))

priority.results.nonlinear = data.frame()
priority.results.nonlinear[1,1] =NA

for (i in 1:1000) {
  
  df.Y.0 <- results_df.nonlinear[(results_df.nonlinear$param_set==i) & (results_df.nonlinear$initial_condition==0),]
  df.Y.1 <- results_df.nonlinear[(results_df.nonlinear$param_set==i) & (results_df.nonlinear$initial_condition==1),]
  if  (((df.Y.0[980,4]-df.Y.0[970,4])*(df.Y.1[980,4]-df.Y.1[970,4]))<0){
    priority.results.nonlinear[i,1]=i
  }
  if  (((df.Y.0[980,5]-df.Y.0[970,5])*(df.Y.1[980,5]-df.Y.1[970,5]))<0){
    priority.results.nonlinear[i,2]=i
  }
  if  (((df.Y.0[980,5]-df.Y.0[970,5])*(df.Y.1[980,5]-df.Y.1[970,5]))==0){
    priority.results.nonlinear[i,3]=i
  }
  if  (((df.Y.0[980,5]-df.Y.0[970,5])*(df.Y.1[980,5]-df.Y.1[970,5]))==0){
    priority.results.nonlinear[i,4]=i
  }
}

priority.results.nonlinear <- priority.results.nonlinear[rowSums(is.na(priority.results.nonlinear)) != ncol(priority.results.nonlinear), ]
priority.results.nonlinear.final <- na.omit(priority.results.nonlinear)



##invasion growth rate
chemostat_model_invasion_growth_rate_N1 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    r11 = r11_max * R1# / (Ks11+R1)
    r12 = r12_max * R2# / (Ks12+R2)
    dR1_dt = D*(R1_in - R1) - r11*N1/Y11
    dR2_dt = D*(R2_in - R2) - r12*N1/Y12
    if (mode == "substitu") {
      dN1_dt = (r11 + r12 - D) * N1 

    } else if (mode == "essen") {
      dN1_dt = (min(r11, r12) - D) * N1 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt))
  })
}
chemostat_model_invasion_growth_rate_N2 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    r21 = r21_max * R1# / (Ks11+R1)
    r22 = r22_max * R2# / (Ks12+R2)
    dR1_dt = D*(R1_in - R1) -  r21*N2/Y21
    dR2_dt = D*(R2_in - R2) -  r22*N2/Y22
    if (mode == "substitu") {
      dN2_dt = (r21 + r22 - D) * N2 
      
    } else if (mode == "essen") {
      dN2_dt = (min(r21, r22) - D) * N2 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN2_dt))
  })
}


set.seed(121)
set.seed(456)
n_samples <- 1000
param_df_invasion <- data.frame(
  r11_max = runif(n_samples, 0.1, 1),
  r12_max = runif(n_samples, 0.1, 1),
  Y11 = runif(n_samples, 0.1, 1),
  Y12 = runif(n_samples, 0.1, 1),
  r21_max = runif(n_samples, 0.1, 1),
  r22_max = runif(n_samples, 0.1, 1),
  Y21 = runif(n_samples, 0.1, 1),
  Y22 = runif(n_samples, 0.1, 1)
)




set.seed(123)
n_samples <- 1000
param_df_invasion <- data.frame(
  r11_max = numeric(n_samples),
  r12_max = numeric(n_samples),
  Y11 = numeric(n_samples),
  Y12 = numeric(n_samples),
  r21_max = numeric(n_samples),
  r22_max = numeric(n_samples),
  Y21 = numeric(n_samples),
  Y22 = numeric(n_samples)
)

for (i in 1:n_samples) {
  sample_vals <- sample(1:100, 8, replace = TRUE) / 100  # Normalize to [0,1]
  param_df_invasion[i, ] <- sample_vals
}



#N1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N1 = c(R1 = 10, R2 = 10, N1 = 0.1)

results.N1 <- lapply(1:nrow(param_df_invasion[,1:4]), function(i) {
    params <- c(param_df_invasion[i,1:4],fixed_params)
    sim <- ode(y = state_chemo_N1 , times = times_chemo, func = chemostat_model_invasion_growth_rate_N1, parms = params, mode = "substitu")
    sim_df <- as.data.frame(sim)
    sim_df$param_set <- i  # Track parameter set
    sim_df
    
})

results_df.N1 <- bind_rows(results.N1)

results_df.N1.equili = results_df.N1[results_df.N1$time==200.0,]


priority.results.invasion = data.frame()
for (i in 1:1000){
  params.2 <- param_df_invasion[i,5:8]
  params.1 <- results_df.N1.equili[i,]
  D=0.1
  
  if (((params.2$r21_max)*(params.1$R1)+(params.2$r22_max)*(params.1$R2)-D)<0){
    priority.results.invasion[i,1] = i 
  }
  else  {
    priority.results.invasion[i,1] = NA 
  }
  
}


#N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N2 = c(R1 = 10, R2 = 10, N2 = 0.1)

results.N2 <- lapply(1:nrow(param_df_invasion[,5:8]), function(i) {
  params <- c(param_df_invasion[i,5:8],fixed_params)
  sim <- ode(y = state_chemo_N2 , times = times_chemo, func = chemostat_model_invasion_growth_rate_N2, parms = params, mode = "substitu")
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N2 <- bind_rows(results.N2)

results_df.N2.equili = results_df.N2[results_df.N2$time==200.0,]


for (i in 1:1000){
  params.2 <- param_df_invasion[i,1:4]
  params.1 <- results_df.N2.equili[i,]
  D=0.1
  
  if (((params.2$r11_max)*(params.1$R1)+(params.2$r12_max)*(params.1$R2)-D)<0){
    priority.results.invasion[i,2] = i 
  }
  else  {
    priority.results.invasion[i,2] = NA 
  }
  
}


library(tidyverse)

#######
#
priority.results.invasion %>% filter(!is.na(V1), !is.na(V2)) %>% dim

#priority.results.invasion = priority.results.invasion[rowSums(is.na(priority.results.invasion))==2,]

priority.results.invasion = na.omit( priority.results.invasion)





## ## ## ## 
## ## ## ## 
## invasion growth rate in nonlinear growth
## ## ## ## 
## ## ## ## 

chemostat_model_invasion_growth_rate_nonlinear_N1 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    X11 = c11*R1/(s11+R1)
    X12 = c12*R2/(s12+R2)
    dR1_dt = D*(R1_in - R1) - X11* N1 
    dR2_dt = D*(R2_in - R2) - X12* N1 
    if (mode == "substitu") {
      dN1_dt = (w11*X11/(q11+X11)+  w12*X12/(q12+X12) - D) * N1 
      
    } else if (mode == "essen") {
      dN1_dt = (min((w11*X11/(q11+X11)), (w12*X12/(q12+X12))) - D) * N1 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt))
  })
}
chemostat_model_invasion_growth_rate_nonlinear_N2 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    X21 = c21*R1/(s21+R1)
    X22 = c22*R2/(s22+R2)
    dR1_dt = D*(R1_in - R1) - X21* N2
    dR2_dt = D*(R2_in - R2) - X22* N2 
    if (mode == "substitu") {
      dN2_dt = (w21*X21/(q21+X21)+  w22*X22/(q22+X22) - D) * N2
      
    } else if (mode == "essen") {
      dN2_dt = (min((w21*X21/(q21+X21)), (w22*X22/(q22+X22))) - D) * N2
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN2_dt))
  })
}


set.seed(123)
set.seed(124)
set.seed(125)
set.seed(126)
set.seed(127)

n_samples <- 1000
param_df_invasion_nonlinear <- data.frame(
  c11 = runif(n_samples, 0, 1),
  s11 = runif(n_samples, 0, 1),
  c12 = runif(n_samples, 0, 1),
  s12 = runif(n_samples, 0, 1),
  w11 = runif(n_samples, 0, 1),
  q11 = runif(n_samples, 0, 1),
  w12 = runif(n_samples, 0, 1),
  q12 = runif(n_samples, 0, 1),
  c21 = runif(n_samples, 0, 1),
  s21 = runif(n_samples, 0, 1),
  c22 = runif(n_samples, 0, 1),
  s22 = runif(n_samples, 0, 1),
  w21 = runif(n_samples, 0, 1),
  q21 = runif(n_samples, 0, 1),
  w22 = runif(n_samples, 0, 1),
  q22 = runif(n_samples, 0, 1)
)


#N1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N1 = c(R1 = 10, R2 = 10, N1 = 0.1)

results.N1 <- lapply(1:nrow(param_df_invasion_nonlinear[,1:8]), function(i) {
  params <- c(param_df_invasion_nonlinear[i,1:8],fixed_params)
  sim <- ode(y = state_chemo_N1 , times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_N1, parms = params, mode = "essen")  ##essen##substitu
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N1 <- bind_rows(results.N1)

results_df.N1.equili = results_df.N1 %>% filter(time==200.0,N1>0.1)

#results_df.N1.equili = results_df.N1[results_df.N1$time==200.0,]
priority.results.invasion_nonlinear = data.frame()

for (i in 1:nrow(results_df.N1.equili)){
  j = results_df.N1.equili[i,5]
  params.2 <- param_df_invasion_nonlinear[j,9:16]
  params.1 <- results_df.N1.equili %>% filter(param_set==j)
  D=0.1
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
    X21 = c21*R1/(s21+R1)
    X22 = c22*R2/(s22+R2)
    
#  if ((w21*X21/(q21+X21)+ w22*X22/(q22+X22) - D)<0){       ##substitute    
  if ((min((w21*X21/(q21+X21)), (w22*X22/(q22+X22))) - D)<0){       ##essential    
    
    
    priority.results.invasion_nonlinear[j,1] = j 
  }
  else  {
    priority.results.invasion_nonlinear[j,1] = NA 
  }
    priority.results.invasion_nonlinear
  })
}




#X21 = (params.2$c21)*(params.1$R1)/((params.2$s21)+(params.1$R1))
#X22 = (params.2$c22)*(params.1$R2)/((params.2$s22)+(params.1$R2))
#(params.2$w21)*X21/((params.2$q21)+X21)+ (params.2$w22)*X22/((params.2$q22)+X22) - D

#N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N2 = c(R1 = 10, R2 = 10, N2 = 0.1)

results.N2 <- lapply(1:nrow(param_df_invasion_nonlinear[,9:16]), function(i) {
  params <- c(param_df_invasion_nonlinear[i,9:16],fixed_params)
  sim <- ode(y = state_chemo_N2 , times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_N2, parms = params, mode = "essen") ##essen##substitu
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N2 <- bind_rows(results.N2)
results_df.N2.equili = results_df.N2 %>% filter(time==200.0,N2>0.1)
#results_df.N2.equili = results_df.N2[results_df.N2$time==200.0,]

for (i in 1:nrow(results_df.N2.equili)){
    j = results_df.N2.equili[i,5]
    params.2 <- param_df_invasion_nonlinear[j,1:8]
    params.1 <- results_df.N2.equili %>% filter(param_set==j)
    D=0.1 
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
    X11 = c11*R1/(s11+R1)
    X12 = c12*R2/(s12+R2)

      
#    if ((w11*X11/(q11+X11)+  w12*X12/(q12+X12) - D)<0){            ##substitute
    if (min((w11*X11/(q11+X11)), (w12*X12/(q12+X12))) - D<0){       ##essential    
   
      priority.results.invasion_nonlinear[j,2] = j
    }
    else  {
      priority.results.invasion_nonlinear[j,2] = NA 
    }
    priority.results.invasion_nonlinear
  })
}



priority.results.invasion_nonlinear.final = na.omit( priority.results.invasion_nonlinear)

priority.results.invasion_nonlinear.match = merge(results_df.N2.equili,results_df.N1.equili,by="param_set",all = F)
priority.results.invasion_nonlinear.match
####### coexistence
priority.results.invasion_nonlinear %>% filter(!is.na(V1), !is.na(V2)) %>% dim




#priority.results.invasion = priority.results.invasion[rowSums(is.na(priority.results.invasion))==2,]
data.nonlinear = data.frame(
       set = rep(1:5, each = 3),
       R2_R1 = rep(c(1, 2, 3), times = 5),
       value = rnorm(15, mean = rep(c(1, 2, 3), times = 5))  )
data.nonlinear$value = c(22,15,17,
                         11,10,9,
                         12,8,9,
                         15,12,16,
                         16,22,20)

anova_model <- aov(value ~ as.factor(R2_R1), data = data.nonlinear)
summary(anova_model)
ggplot(data.nonlinear, aes(x = R2_R1, y = value, color = as.factor(set))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


##
##
##
##
## test codes nonlinear
chemostat_model_invasion_growth_rate_nonlinear_N1N2 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    X11 = c11*R1/(s11+R1)
    X12 = c12*R2/(s12+R2)
    X21 = c21*R1/(s21+R1)
    X22 = c22*R2/(s22+R2)
    
    dR1_dt = D*(R1_in - R1) - X21* N2 - X11* N1
    dR2_dt = D*(R2_in - R2) - X22* N2 - X12* N1
    if (mode == "substitu") {
      dN2_dt = (w21*X21/(q21+X21)+  w22*X22/(q22+X22) - D) * N2
      dN1_dt = (w11*X11/(q11+X11)+  w12*X12/(q12+X12) - D) * N1
      
    } else if (mode == "essen") {
      dN2_dt = (min((w21*X21/(q21+X21)), (w22*X22/(q22+X22))) - D) * N2
      dN1_dt = (min((w11*X11/(q11+X11)), (w12*X12/(q12+X12))) - D) * N1
      
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt,dN2_dt))
  })
}

parameters_chemo = c(c11 = 0.7445759  ,s11 = 0.2133924   ,w11 = 0.1008311   ,q11 = 0.228527   , 
                     c12 = 0.4972046   ,s12 = 0.4546738   ,w12 = 0.6009141   ,q12 = 0.5336715  , 
                     c21 = 0.333766   ,s21 = 0.4559487   ,w21 = 0.9793974   ,q21 = 0.2761582  , 
                     c22 = 0.8024584  ,s22 = 0.7254989   ,w22 = 0.2637167   ,q22 = 0.5339131, 
                     D = 0.1,R1_in = 10,R2_in = 10)


state_chemo = c(R1 = 10,R2=10, N1 = 1,N2 = .1)
#state_chemo = c(R1 = 10,R2=10, N1 = .1,N2 = 1)

times_chemo = seq(0,500,by = 0.1)

ode_chemo = ode (y = state_chemo, times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_N1N2, parms = parameters_chemo,mode='substitu')

out_df_chemo = as.data.frame(ode_chemo)


p = ggplot(out_df_chemo, aes(x = time)) +
  geom_line(aes(y = N1, color = "Species_1"),linewidth = 1) +
  geom_line(aes(y = N2, color = "Species_2"),linewidth = 1) +
  #geom_line(aes(y = R1, color = "Resource_1"),linewidth = 1) +
  #geom_line(aes(y = R2, color = "Resource_2"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation 2 species 2 resource",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  scale_y_log10() +
  theme_minimal() 

setwd("D:/zhang lab/try")
ggsave(filename = paste0("chemostat_sim_1",  ".png"), plot = p, width = 8, height = 6)

#ggsave(filename = paste0("chemostat_sim_0",  ".png"), plot = p, width = 8, height = 6)



#####
#####
#####
#####
#####
for (i in 1:1000){
  params.2 <- param_df_invasion_nonlinear[i,9:16]
  params.1 <- results_df.N1.equili[i,]
  D=0.1
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
    X21 = c21*R1/(s21+R1)
    X22 = c22*R2/(s22+R2)
    
    if ((w21*X21/(q21+X21)+ w22*X22/(q22+X22) - D)<0){
      priority.results.invasion_nonlinear[i,1] =i 
    }
    else  {
      priority.results.invasion_nonlinear[i,1] = NA 
    }
    priority.results.invasion_nonlinear
  })
}


for (i in 1:1000){
params.2 <- param_df_invasion_nonlinear[i,1:8]
params.1 <- results_df.N2.equili[i,]
D=0.1 

priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
  X11 = c11*R1/(s11+R1)
  X12 = c12*R2/(s12+R2)
  
  
  if ((w11*X11/(q11+X11)+  w12*X12/(q12+X12) - D)<0){
    priority.results.invasion_nonlinear[i,2] = i
  }
  else  {
    priority.results.invasion_nonlinear[i,2] = NA 
  }
  priority.results.invasion_nonlinear
})
}







## test codes


parameters_chemo = c(r11_max = 0.89,r12_max = 0.33,r21_max = 0.74,r22_max = 0.97, D = 0.1, 
                     Y11 = 0.66,Y12 = 0.04,Y21 = 0.05,Y22 = 0.25,R1_in = 10,R2_in = 10)


state_chemo = c(R1 = 10,R2=10, N1 = 3,N2 = .1)
#state_chemo = c(R1 = 10,R2=10, N1 = .1,N2 = 3)

times_chemo = seq(0,500,by = 0.1)

ode_chemo = ode (y = state_chemo, times = times_chemo, func = chemostat_model_2s2r, parms = parameters_chemo,mode='substitu')

out_df_chemo = as.data.frame(ode_chemo)


p = ggplot(out_df_chemo, aes(x = time)) +
  geom_line(aes(y = N1, color = "Species_1"),linewidth = 1) +
 # geom_line(aes(y = N2, color = "Species_2"),linewidth = 1) +
  geom_line(aes(y = R1, color = "Resource_1"),linewidth = 1) +
  geom_line(aes(y = R2, color = "Resource_2"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation 2 species 2 resource",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  scale_y_log10() +
  theme_minimal() 

setwd("D:/zhang lab/try")
ggsave(filename = paste0("chemostat_sim_1",  ".png"), plot = p, width = 8, height = 6)

#ggsave(filename = paste0("chemostat_sim_0",  ".png"), plot = p, width = 8, height = 6)




## test codes of experiment 

chemostat_model_invasion_growth_rate_nonlinear_experiment = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{

    dR1_dt = D*(R1_in - R1) - c21*N2*R1/(s21+R1) - c11*N1*R1/(s11+R1)
    dR2_dt = D*(R2_in - R2) - c22*N2*R2/(s22+R2) - c12*N1*R2/(s12+R2)
    if (mode == "substitu") {
      dN2_dt = (u21*R1/(k21+R1)+  u22*R2/(k22+R2) - D) * N2
      dN1_dt = (u11*R1/(k11+R1)+  u12*R2/(k12+R2) - D) * N1
      
    } else if (mode == "essen") {
      dN2_dt = (min((u21*R1/(k21+R1)), (u22*R2/(k22+R2))) - D) * N2
      dN1_dt = (min((u11*R1/(k11+R1)), (u12*R2/(k12+R2))) - D) * N1
      
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt,dN2_dt))
  })
}

parameters_chemo = c(u11 = 1.768724   ,k11 = 0.2897989      ,c11 = 11.32309     ,s11 = 0.2346575      , 
                     u12 = 1.815949     ,k12 = 0.5255631      ,c12 = 12.26484     ,s12 = 1.552779     , 
                     u21 = 1.455718     ,k21 = 0.08491354      ,c21 = 13.44069     ,s21 = 7.942029    , 
                     u22 = 1.634622   ,k22 = 1.0439      ,c22 = 2.907819      ,s22 = 0.8035113 , 
                     D = 0.1,R1_in = 10,R2_in = 10)


state_chemo = c(R1 = 10,R2=10, N1 = 6,N2 = .1)
#state_chemo = c(R1 = 10,R2=10, N1 = .1,N2 = 6)

times_chemo = seq(0,500,by = 0.1)

ode_chemo = ode (y = state_chemo, times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_experiment, parms = parameters_chemo,mode='substitu')

out_df_chemo = as.data.frame(ode_chemo)


p = ggplot(out_df_chemo, aes(x = time)) +
  geom_line(aes(y = N1, color = "Species_1"),linewidth = 1) +
  geom_line(aes(y = N2, color = "Species_2"),linewidth = 1) +
  #geom_line(aes(y = R1, color = "Resource_1"),linewidth = 1) +
  #geom_line(aes(y = R2, color = "Resource_2"),linewidth = 1) +
  labs(title = "Chemostat Model Simulation 2 species 2 resource",
       x = "Time",
       y = "Concentration",
       color = "Legend") +
  scale_y_log10() +
  theme_minimal() 
p
setwd("D:/zhang lab/try")
ggsave(filename = paste0("chemostat_sim_1",  ".png"), plot = p, width = 8, height = 6)

#ggsave(filename = paste0("chemostat_sim_0",  ".png"), plot = p, width = 8, height = 6)


## accuracy of experiment 
## accuracy of experiment 
## accuracy of experiment 
## accuracy of experiment 
load('paras_2yansong.RData')
para = data.frame()
para = t_res %>% filter(sp %in% c('c','13','16','38','65')) %>% select(sp, u_nh,k_nh,c_nh,s_nh,u_no,k_no,c_no,s_no)
para$sp = c('2','3','4','5','1')

species_pairs <- combn(para$sp, 2) 

paired_para <- do.call(rbind, apply(species_pairs, 2, function(pair) {
  species1 <- para[para$sp == pair[1], ]
  species2 <- para[para$sp == pair[2], ]
  combined <- cbind(species1, species2)
  combined$lst_sp <- paste(pair[1], pair[2], sep = "_")
  
  return(combined)
}))
paired_para = paired_para[,!colnames(paired_para)=='sp']
colnames(paired_para) = c('u11','k11','c11','s11','u12','k12','c12','s12','u21','k21','c21','s21','u22','k22','c22','s22','sp')



chemostat_model_invasion_growth_rate_nonlinear_experiment_N1 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{

    dR1_dt = D*(R1_in - R1) - c11* N1 *R1/(R1+s11) 
    dR2_dt = D*(R2_in - R2) - c12* N1 *R2/(R2+s12)
    if (mode == "substitu") {
      dN1_dt = (u11*R1/(k11+R1)+  u12*R2/(k12+R2) - D) * N1 
      
    } else if (mode == "essen") {
      dN1_dt = (min((u11*R1/(k11+R1)), (u12*R2/(k12+R2))) - D) * N1 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt))
  })
}
chemostat_model_invasion_growth_rate_nonlinear_experiment_N2 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    dR1_dt = D*(R1_in - R1) - c21*N2*R1/(R1+s21)
    dR2_dt = D*(R2_in - R2) - c22*N2*R2/(R2+s22) 
    if (mode == "substitu") {
      dN2_dt = (u21*R1/(k21+R1)+  u22*R2/(k22+R2) - D) * N2
      
    } else if (mode == "essen") {
      dN2_dt = (min((u21*R1/(k21+R1)), (u22*R2/(k22+R2))) - D) * N2
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN2_dt))
  })
}


#N1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N1 = c(R1 = 10, R2 = 10, N1 = 0.1)

results.N1 <- lapply(1:nrow(paired_para[,1:8]), function(i) {
  params <- c(paired_para[i,1:8],fixed_params)
  sim <- ode(y = state_chemo_N1 , times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_experiment_N1, parms = params, mode = "substitu")  ##essen##substitu
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N1 <- bind_rows(results.N1)

results_df.N1.equili = results_df.N1 %>% filter(time==200.0,N1>0.1)

#results_df.N1.equili = results_df.N1[results_df.N1$time==200.0,]
priority.results.invasion_nonlinear = data.frame()

for (i in 1:nrow(results_df.N1.equili)){
  j = results_df.N1.equili[i,5]
  params.2 <- paired_para[j,9:16]
  params.1 <- results_df.N1.equili %>% filter(param_set==j)
  D=0.1
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{

    
      if ((u21*R1/(k21+R1)+  u22*R2/(k22+R2) - D)<0){       ##substitute    
   # if ((min((u21*R1/(k21+R1)), (u22*R2/(k22+R2))) - D)<0){       ##essential    
      
      
      priority.results.invasion_nonlinear[j,1] = j 
    }
    else  {
      priority.results.invasion_nonlinear[j,1] = NA 
    }
    priority.results.invasion_nonlinear
  })
}




#X21 = (params.2$c21)*(params.1$R1)/((params.2$s21)+(params.1$R1))
#X22 = (params.2$c22)*(params.1$R2)/((params.2$s22)+(params.1$R2))
#(params.2$w21)*X21/((params.2$q21)+X21)+ (params.2$w22)*X22/((params.2$q22)+X22) - D

#N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N2 = c(R1 = 10, R2 = 10, N2 = 0.1)

results.N2 <- lapply(1:nrow(paired_para[,9:16]), function(i) {
  params <- c(paired_para[i,9:16],fixed_params)
  sim <- ode(y = state_chemo_N2 , times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_experiment_N2, parms = params, mode = "substitu") ##essen##substitu
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N2 <- bind_rows(results.N2)
results_df.N2.equili = results_df.N2 %>% filter(time==200.0,N2>0.1)
#results_df.N2.equili = results_df.N2[results_df.N2$time==200.0,]

for (i in 1:nrow(results_df.N2.equili)){
  j = results_df.N2.equili[i,5]
  params.2 <- paired_para[j,1:8]
  params.1 <- results_df.N2.equili %>% filter(param_set==j)
  D=0.1 
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
    
        if ((u11*R1/(k11+R1)+  u12*R2/(k12+R2) - D)<0){            ##substitute
   # if ((min((u11*R1/(k11+R1)), (u12*R2/(k12+R2))) - D)<0){       ##essential    
      
      priority.results.invasion_nonlinear[j,2] = j
    }
    else  {
      priority.results.invasion_nonlinear[j,2] = NA 
    }
    priority.results.invasion_nonlinear
  })
}



priority.results.invasion_nonlinear.final = na.omit( priority.results.invasion_nonlinear)

priority.results.invasion_nonlinear.match = merge(results_df.N2.equili,results_df.N1.equili,by="param_set",all = F)
priority.results.invasion_nonlinear.match
####### coexistence
priority.results.invasion_nonlinear %>% filter(!is.na(V1), !is.na(V2)) %>% dim



## ## ## ## 
## ## ## ## 
## simulation in discrete time
## ## ## ## 
## ## ## ## 
set.seed(111)
set.seed(222)
set.seed(333)


R1_in <- 50
R2_in <- 50
n_samples <- 1000
param_df_discrete_nonlinear <- data.frame(
  u11 = runif(n_samples, 0, 1),
  d11 = runif(n_samples, 0, 0.1),
  u12 = runif(n_samples, 0, 1),
  d12 = runif(n_samples, 0, 0.1),
  c11 = runif(n_samples, 0, 1),
  k11 = runif(n_samples, 0, 1),
  s11 = runif(n_samples, 0, 1),
  c12 = runif(n_samples, 0, 1),
  k12 = runif(n_samples, 0, 1),
  s12 = runif(n_samples, 0, 1),
  u21 = runif(n_samples, 0, 1),
  d21 = runif(n_samples, 0, 0.1),
  u22 = runif(n_samples, 0, 1),
  d22 = runif(n_samples, 0, 0.1),
  c21 = runif(n_samples, 0, 1),
  k21 = runif(n_samples, 0, 1),
  s21 = runif(n_samples, 0, 1),
  c22 = runif(n_samples, 0, 1),
  k22 = runif(n_samples, 0, 1),
  s22 = runif(n_samples, 0, 1)
)
param_df_discrete_nonlinear$c11 = 1-param_df_discrete_nonlinear$u11
param_df_discrete_nonlinear$c12 = 1-param_df_discrete_nonlinear$u12
param_df_discrete_nonlinear$c21 = 1-param_df_discrete_nonlinear$u21
param_df_discrete_nonlinear$c22 = 1-param_df_discrete_nonlinear$u22







para = data.frame()
para = t_res %>% filter(sp %in% c('c','13','16','38','65')) %>% select(sp, u_nh,k_nh,d_nh,c_nh,s_nh,u_no,k_no,d_no,c_no,s_no,)
para$sp = c('2','3','4','5','1')

species_pairs <- combn(para$sp, 2) 

paired_para <- do.call(rbind, apply(species_pairs, 2, function(pair) {
  species1 <- para[para$sp == pair[1], ]
  species2 <- para[para$sp == pair[2], ]
  combined <- cbind(species1, species2)
  combined$lst_sp <- paste(pair[1], pair[2], sep = "_")
  
  return(combined)
}))
paired_para = paired_para[,!colnames(paired_para)=='sp']
colnames(paired_para) = c('u11','k11','d11','c11','s11','u12','k12',"d12",'c12','s12','u21','k21',"d21",'c21','s21','u22','k22','d22','c22','s22','sp')



simulate_discrete_chemostat <- function(params,N1_0,N2_0,R1_0,R2_0, timesteps = 50, D = 0.8) {
  with(params, {
    
  R1 <- R2 <- N1 <- N2 <- numeric(timesteps)
  R1[1] <- R1_0   
  R2[1] <- R2_0 
  N1[1] <- N1_0
  N2[1] <- N2_0

  for (t in 1:(timesteps - 1)) {
   n1 = numeric(4)
   n2 = numeric(4)
   r1 = numeric(4)
   r2 = numeric(4)
    n1[1] = N1[t]
    n2[1] = N2[t]
    r1[1] = R1[t]
    r2[1] = R2[t]
    for (i in 1:3){
      n1[i+1] = exp((u11*r1[i]/(r1[i]+k11) + u12*r2[i]/(r2[i]+k12) - d12 - d11)) * n1[i]
      n2[i+1] = exp((u21*r1[i]/(r1[i]+k21) + u22*r2[i]/(r2[i]+k22) - d22 - d21)) * n2[i]
      r1[i+1] = exp( (-c11* n1[i]/(r1[i]+s11)-c21* n2[i]/(r1[i]+s21)) )* r1[i]
      r2[i+1] = exp( (-c12* n1[i]/(r2[i]+s12)-c22* n2[i]/(r2[i]+s22)) )* r2[i]
    }
    n1_h = exp(0.5*(u11*r1[4]/(r1[4]+k11) + u12*r2[4]/(r2[4]+k12) - d12- d11) )* n1[4]
    n2_h = exp(0.5*(u21*r1[4]/(r1[4]+k21) + u22*r2[4]/(r2[4]+k22) - d22- d21) )* n2[4]
    
    r1_h = exp(-0.5* (c11*n1[4]/(r1[4]+s11)+  c21*r1[i]/(r1[i]+s21)) ) * r1[4]
    r2_h = exp(-0.5* (c12*n2[4]/(r2[4]+s12)+  c22*r1[i]/(r1[i]+s22)) ) * r2[4]
    
    N1[t+1] <- (1 - D) * n1_h
    N2[t+1] <- (1 - D) * n2_h
    R1[t+1] <- (1 - D) * r1_h + D * R1_in
    R2[t+1] <- (1 - D) * r2_h + D * R2_in
    
  }
  data.frame(
    timestep = 1:timesteps,
    N_1 = N1,
    N_2 = N2,
    R_1 = R1,
    R_2 = R2
  )  
  
  
  })
}
#1
results_list.1 <- vector("list", n_samples)

for (i in 1:n_samples) {
  params <- param_df_discrete_nonlinear[i, ]
  #params <- paired_para[i, ]
  
  results_list.1[[i]] <- simulate_discrete_chemostat(params,N1_0=0.1,N2_0=100,R1_0=50,R2_0=50)
  results_list.1[[i]]$SimID <- i
}

results_discrete_1 <- bind_rows(results_list.1)
#2
results_list.2 <- vector("list", n_samples)

for (i in 1:n_samples) {
  params <- param_df_discrete_nonlinear[i, ]
  #params <- paired_para[i, ]
  
  results_list.2[[i]] <- simulate_discrete_chemostat(params,N1_0=100,N2_0=0.1,R1_0=50,R2_0=50)
  results_list.2[[i]]$SimID <- i
}

results_discrete_2 <- bind_rows(results_list.2)

r_1 = results_discrete_1 %>% filter(timestep==50) %>%  mutate(Dominant = case_when(
 (N_1 > 1000*N_2)&(N_1 > 0.1) ~ 1,
 (N_1 < N_2/1000)&(N_2 > 0.1) ~ 2,
 (N_1 < 0.1)&(N_2 < 0.1) ~ -100,
  TRUE ~ 100))
  

r_2 = results_discrete_2 %>% filter(timestep==50) %>%  mutate(Dominant = case_when(
 (N_1 > 1000*N_2)&(N_1 > 0.1) ~ 1,
 (N_1 < N_2/1000)&(N_2 > 0.1) ~ 2,
 (N_1 < 0.1)&(N_2 < 0.1) ~ -100,
  TRUE ~ 100))
sum((r_1$Dominant - r_2$Dominant)==1)


r_combined = r_2 %>% mutate(Dominant_r_1 = r_1$Dominant) %>% filter((Dominant==-100)|(Dominant_r_1==-100))
r_combined %>% dim






# Plot
ggplot(df, aes(x = time)) +
  geom_line(aes(y = Resource, color = "Resource"), linewidth = 1.2) +
  geom_line(aes(y = Biomass, color = "Biomass"), linewidth = 1.2) +
  labs(title = "Discrete-Time Chemostat Simulation",
       x = "Time (steps)",
       y = "Concentration",
       color = "Variable") +
  theme_minimal()

 

ggplot(t_res, aes(x = u_nh, y = c_nh)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


library(corrplot)

# Plot correlation matrix
numeric_df <- t_res[sapply(t_res, is.numeric)]


cor_matrix <- cor(numeric_df[,1:10], use = "complete.obs")
print(cor_matrix)
which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8)








chemostat_model_invasion_growth_rate_nonlinear_experiment_discrete_N1 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    
    dR1_dt = D*(R1_in - R1) - c11* N1 *R1/(R1+s11) 
    dR2_dt = D*(R2_in - R2) - c12* N1 *R2/(R2+s12)
    if (mode == "substitu") {
      dN1_dt = (u11*R1/(k11+R1)+  u12*R2/(k12+R2) - D) * N1 
      
    } else if (mode == "essen") {
      dN1_dt = (min((u11*R1/(k11+R1)), (u12*R2/(k12+R2))) - D) * N1 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt))
  })
}
chemostat_model_invasion_growth_rate_nonlinear_experiment_discrete_N2 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    dR1_dt = D*(R1_in - R1) - c21*N2*R1/(R1+s21)
    dR2_dt = D*(R2_in - R2) - c22*N2*R2/(R2+s22) 
    if (mode == "substitu") {
      dN2_dt = (u21*R1/(k21+R1)+  u22*R2/(k22+R2) - D) * N2
      
    } else if (mode == "essen") {
      dN2_dt = (min((u21*R1/(k21+R1)), (u22*R2/(k22+R2))) - D) * N2
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN2_dt))
  })
}


#N1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N1 = c(R1 = 10, R2 = 10, N1 = 0.1)

results.N1 <- lapply(1:nrow(paired_para[,1:8]), function(i) {
  params <- c(paired_para[i,1:8],fixed_params)
  sim <- ode(y = state_chemo_N1 , times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_experiment_N1, parms = params, mode = "substitu")  ##essen##substitu
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N1 <- bind_rows(results.N1)

results_df.N1.equili = results_df.N1 %>% filter(time==200.0,N1>0.1)

#results_df.N1.equili = results_df.N1[results_df.N1$time==200.0,]
priority.results.invasion_nonlinear = data.frame()

for (i in 1:nrow(results_df.N1.equili)){
  j = results_df.N1.equili[i,5]
  params.2 <- paired_para[j,9:16]
  params.1 <- results_df.N1.equili %>% filter(param_set==j)
  D=0.1
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
    
    
    if ((u21*R1/(k21+R1)+  u22*R2/(k22+R2) - D)<0){       ##substitute    
      # if ((min((u21*R1/(k21+R1)), (u22*R2/(k22+R2))) - D)<0){       ##essential    
      
      
      priority.results.invasion_nonlinear[j,1] = j 
    }
    else  {
      priority.results.invasion_nonlinear[j,1] = NA 
    }
    priority.results.invasion_nonlinear
  })
}




#X21 = (params.2$c21)*(params.1$R1)/((params.2$s21)+(params.1$R1))
#X22 = (params.2$c22)*(params.1$R2)/((params.2$s22)+(params.1$R2))
#(params.2$w21)*X21/((params.2$q21)+X21)+ (params.2$w22)*X22/((params.2$q22)+X22) - D

#N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N2 = c(R1 = 10, R2 = 10, N2 = 0.1)

results.N2 <- lapply(1:nrow(paired_para[,9:16]), function(i) {
  params <- c(paired_para[i,9:16],fixed_params)
  sim <- ode(y = state_chemo_N2 , times = times_chemo, func = chemostat_model_invasion_growth_rate_nonlinear_experiment_N2, parms = params, mode = "substitu") ##essen##substitu
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N2 <- bind_rows(results.N2)
results_df.N2.equili = results_df.N2 %>% filter(time==200.0,N2>0.1)
#results_df.N2.equili = results_df.N2[results_df.N2$time==200.0,]

for (i in 1:nrow(results_df.N2.equili)){
  j = results_df.N2.equili[i,5]
  params.2 <- paired_para[j,1:8]
  params.1 <- results_df.N2.equili %>% filter(param_set==j)
  D=0.1 
  
  priority.results.invasion_nonlinear = with(as.data.frame(c(params.1,params.2)),{
    
    if ((u11*R1/(k11+R1)+  u12*R2/(k12+R2) - D)<0){            ##substitute
      # if ((min((u11*R1/(k11+R1)), (u12*R2/(k12+R2))) - D)<0){       ##essential    
      
      priority.results.invasion_nonlinear[j,2] = j
    }
    else  {
      priority.results.invasion_nonlinear[j,2] = NA 
    }
    priority.results.invasion_nonlinear
  })
}



priority.results.invasion_nonlinear.final = na.omit( priority.results.invasion_nonlinear)

priority.results.invasion_nonlinear.match = merge(results_df.N2.equili,results_df.N1.equili,by="param_set",all = F)
priority.results.invasion_nonlinear.match
####### coexistence
priority.results.invasion_nonlinear %>% filter(!is.na(V1), !is.na(V2)) %>% dim



## ## ## ## 
## ## ## ## 
## simulation in discrete time  reciprocal invasion experiments
## ## ## ## 
## ## ## ## 
set.seed(111)
set.seed(222)

R1_in <- 50
R2_in <- 50
n_samples <- 1000
param_df_discrete_nonlinear <- data.frame(
  u11 = runif(n_samples, 0, 1),
  d11 = runif(n_samples, 0, 0.1),
  u12 = runif(n_samples, 0, 1),
  d12 = runif(n_samples, 0, 0.1),
  c11 = runif(n_samples, 0, 1),
  k11 = runif(n_samples, 0, 1),
  s11 = runif(n_samples, 0, 1),
  c12 = runif(n_samples, 0, 1),
  k12 = runif(n_samples, 0, 1),
  s12 = runif(n_samples, 0, 1),
  u21 = runif(n_samples, 0, 1),
  d21 = runif(n_samples, 0, 0.1),
  u22 = runif(n_samples, 0, 1),
  d22 = runif(n_samples, 0, 0.1),
  c21 = runif(n_samples, 0, 1),
  k21 = runif(n_samples, 0, 1),
  s21 = runif(n_samples, 0, 1),
  c22 = runif(n_samples, 0, 1),
  k22 = runif(n_samples, 0, 1),
  s22 = runif(n_samples, 0, 1)
)

simulate_discrete_chemostat_invasion_N1 <- function(params,N1_0,R1_0,R2_0, timesteps = 50, D = 0.8) {
  with(params, {
    
    R1 <- R2 <- N1  <- numeric(timesteps)
    R1[1] <- R1_0   
    R2[1] <- R2_0 
    N1[1] <- N1_0
    #N2[1] <- N2_0

    for (t in 1:(timesteps - 1)) {
      n1 = numeric(4)
      #n2 = numeric(4)
      r1 = numeric(4)
      r2 = numeric(4)
      n1[1] = N1[t]
      #n2[1] = N2[t]
      r1[1] = R1[t]
      r2[1] = R2[t]
      for (i in 1:3){
        n1[i+1] = exp((u11*r1[i]/(r1[i]+k11) + u12*r2[i]/(r2[i]+k12) - d12 - d11)) * n1[i]
        #n2[i+1] = exp((u21*r1[i]/(r1[i]+k21) + u22*r2[i]/(r2[i]+k22) - d22 - d21)) * n2[i]
        r1[i+1] = exp( (-c11* n1[i]/(r1[i]+s11)) )* r1[i]
        r2[i+1] = exp( (-c12* n1[i]/(r2[i]+s12)) )* r2[i]
      }
      n1_h = exp(0.5*(u11*r1[4]/(r1[4]+k11) + u12*r2[4]/(r2[4]+k12) - d12- d11) )* n1[4]
      #n2_h = exp(0.5*(u21*r1[4]/(r1[4]+k21) + u22*r2[4]/(r2[4]+k22) - d22- d21) )* n2[4]
      
      r1_h = exp(-0.5* (c11*n1[4]/(r1[4]+s11)) ) * r1[4]
      r2_h = exp(-0.5* (c12*n2[4]/(r2[4]+s12)) ) * r2[4]
      
      N1[t+1] <- (1 - D) * n1_h
      #N2[t+1] <- (1 - D) * n2_h
      R1[t+1] <- (1 - D) * r1_h + D * R1_in
      R2[t+1] <- (1 - D) * r2_h + D * R2_in
      
    }
    data.frame(
      timestep = 1:timesteps,
      N_1 = N1,
      R_1 = R1,
      R_2 = R2
    )  
    
    
  })
}
#1
results_list.1 <- vector("list", n_samples)

for (i in 1:n_samples) {
  params <- param_df_discrete_nonlinear[i, ]
  results_list.1[[i]] <- simulate_discrete_chemostat(params,N1_0=0.1,N2_0=10,R1_0=50,R2_0=50)
  results_list.1[[i]]$SimID <- i
}

results_discrete_1 <- bind_rows(results_list.1)
#2
results_list.2 <- vector("list", n_samples)

for (i in 1:n_samples) {
  params <- param_df_discrete_nonlinear[i, ]
  results_list.2[[i]] <- simulate_discrete_chemostat(params,N1_0=10,N2_0=0.1,R1_0=50,R2_0=50)
  results_list.2[[i]]$SimID <- i
}

results_discrete_2 <- bind_rows(results_list.2)

r_1 = results_discrete_1 %>% filter(timestep==50) %>%  mutate(Dominant = case_when(
  (N_1 > 1000*N_2)&(N_1 > 0.1) ~ 1,
  (N_1 < N_2/1000)&(N_2 > 0.1) ~ 2,
  (N_1 < 0.1)&(N_2 < 0.1) ~ -100,
  TRUE ~ 100))


r_2 = results_discrete_2 %>% filter(timestep==50) %>%  mutate(Dominant = case_when(
  (N_1 > 1000*N_2)&(N_1 > 0.1) ~ 1,
  (N_1 < N_2/1000)&(N_2 > 0.1) ~ 2,
  (N_1 < 0.1)&(N_2 < 0.1) ~ -100,
  TRUE ~ 100))
sum((r_1$Dominant - r_2$Dominant)==1)


r_combined = r_2 %>% mutate(Dominant_r_1 = r_1$Dominant) %>% filter((Dominant==-100)|(Dominant_r_1==-100))
r_combined %>% dim


## ## ## ## 
## ## ## ## 
## simulation in discrete time linear process
## ## ## ## 
## ## ## ## 
## ## ## ## 
## ## ## ## 
para = data.frame()
para = t_res %>% filter(sp %in% c('c','13','16','38','65')) %>% select(sp, u_nh,k_nh,d_nh,c_nh,s_nh,u_no,k_no,d_no,c_no,s_no,)
para$sp = c('2','3','4','5','1')

species_pairs <- combn(para$sp, 2) 

paired_para <- do.call(rbind, apply(species_pairs, 2, function(pair) {
  species1 <- para[para$sp == pair[1], ]
  species2 <- para[para$sp == pair[2], ]
  combined <- cbind(species1, species2)
  combined$lst_sp <- paste(pair[1], pair[2], sep = "_")
  
  return(combined)
}))
paired_para = paired_para[,!colnames(paired_para)=='sp']
colnames(paired_para) = c('u11','k11','d11','c11','s11','u12','k12',"d12",'c12','s12','u21','k21',"d21",'c21','s21','u22','k22','d22','c22','s22','sp')



set.seed(122)


R1_in <- 50
R2_in <- 50
n_samples <- 1000
param_df_discrete_linear <- data.frame(
  c11 = runif(n_samples, 0.1, 1),
  m11 = runif(n_samples, 0.1, 0.5),
  w11 = runif(n_samples, 0.01, 0.1),
  c12 = runif(n_samples, 0.1, 1),
  m12 = runif(n_samples, 0.1, 0.5),
  w12 = runif(n_samples, 0.01, 0.1),
  c21 = runif(n_samples, 0.1, 1),
  m21 = runif(n_samples, 0.1, 0.5),
  w21 = runif(n_samples, 0.01, 0.1),
  c22 = runif(n_samples, 0.1, 1),
  m22 = runif(n_samples, 0.1, 0.5),
  w22 = runif(n_samples, 0.01, 0.1)
)
param_df_discrete_linear$c11 = 1-param_df_discrete_linear$w11
param_df_discrete_linear$c12 = 1-param_df_discrete_linear$w12
param_df_discrete_linear$c21 = 1-param_df_discrete_linear$w21
param_df_discrete_linear$c22 = 1-param_df_discrete_linear$w22






simulate_discrete_chemostat_linear <- function(params,N1_0,N2_0,R1_0,R2_0, timesteps = 50, D = 0.8) {
  with(params, {
    
    R1 <- R2 <- N1 <- N2 <- numeric(timesteps)
    R1[1] <- R1_0   
    R2[1] <- R2_0 
    N1[1] <- N1_0
    N2[1] <- N2_0
    
    for (t in 1:(timesteps - 1)) {
      n1 = numeric(4)
      n2 = numeric(4)
      r1 = numeric(4)
      r2 = numeric(4)
      n1[1] = N1[t]
      n2[1] = N2[t]
      r1[1] = R1[t]
      r2[1] = R2[t]
      for (i in 1:3){

        r1[i+1] = exp( (-c11* n1[i]-c21*n2[i]) )* r1[i]
        r2[i+1] = exp( (-c12* n1[i]-c22*n2[i]) )* r2[i]
        n1[i+1] = exp((c11*w11*r1[i] + c12*w12*r2[i] - m12 - m11)) * n1[i]
        n2[i+1] = exp((c21*w21*r1[i] + c22*w22*r2[i] - m22 - m21)) * n2[i]
        
      }

      
      r1_h = exp(-0.5* (c11* n1[4]+c21*n2[4]) ) * r1[4]
      r2_h = exp(-0.5* (c12* n1[4]+c22*n2[4]) ) * r2[4]
      n1_h = exp(0.5*(c11*w11*r1[4] + c12*w12*r2[4] - m12 - m11) )* n1[4]
      n2_h = exp(0.5*(c21*w21*r1[4] + c22*w22*r2[4] - m22 - m21) )* n2[4]
      
      N1[t+1] <- (1 - D) * n1_h
      N2[t+1] <- (1 - D) * n2_h
      R1[t+1] <- (1 - D) * r1_h + D * R1_in
      R2[t+1] <- (1 - D) * r2_h + D * R2_in
      
    }
    data.frame(
      timestep = 1:timesteps,
      N_1 = N1,
      N_2 = N2,
      R_1 = R1,
      R_2 = R2
    )  
    
    
  })
}
#1
results_list.1 <- vector("list", n_samples)

for (i in 1:n_samples) {
  params <- param_df_discrete_linear[i, ]
  #params <- paired_para[i, ]
  
  results_list.1[[i]] <- simulate_discrete_chemostat_linear(params,N1_0=0.001,N2_0=0.1,R1_0=50,R2_0=50)
  results_list.1[[i]]$SimID <- i
}

results_discrete_1 <- bind_rows(results_list.1)


#2
results_list.2 <- vector("list", n_samples)

for (i in 1:n_samples) {
  params <- param_df_discrete_linear[i, ]
  #params <- paired_para[i, ]
  
  results_list.2[[i]] <- simulate_discrete_chemostat_linear(params,N1_0=0.1,N2_0=0.001,R1_0=50,R2_0=50)
  results_list.2[[i]]$SimID <- i
}

results_discrete_2 <- bind_rows(results_list.2)

r_1.30 = results_discrete_1 %>% filter(timestep==50) 
r_1.28 = results_discrete_1 %>% filter(timestep==48) 


r_2.30 = results_discrete_2 %>% filter(timestep==50) 
r_2.28 = results_discrete_2 %>% filter(timestep==48) 


priority.results.discret = data.frame()
priority.results.discret[1,1] =NA

for (i in 1:1000){
  if (((r_1.30[i,3]-r_1.28[i,3])>0) & ((r_1.30[i,2]-r_1.28[i,2])<0)){
     priority.results.discret[i,1]=1
  }
  
  if (((r_2.30[i,2]-r_2.28[i,2])>0) & ((r_2.30[i,3]-r_2.28[i,3])<0)){
    priority.results.discret[i,2]=1
  }
  
}
priority.results.discret %>% filter((V1==1)&(V2==1)) %>% dim



###try code
ID_pri <- d %>% filter(priority == 1) %>% pull(SimID)
length(ID_pri)/n_samples
d_raw %>% filter(SimID %in% ID_pri) %>% 
  filter(timestep == max(timestep)) %>% 
  arrange(SimID)



for (i in 1:1000){
  if (((r_2.30[i,2])>10) & ((r_2.30[i,3]>0.10))){
    priority.results.discret[i,1]=1
  }
  else
  {
     priority.results.discret[i,1]=NA
  }
   
  
}
priority.results.discret %>% filter((V1==1)) %>% dim


library(dplyr)
library(ggplot2)
library(patchwork)  # for side-by-side plots

# Assume both dataframes have a column `param_set` to identify the set,
# and columns: time, R1, R2, N1, etc.

# 1. Select 20 random param_set IDs
set.seed(1234)  # for reproducibility
set.seed(1111)
set.seed(2222)
set.seed(3333)

random_ids <- sample(1:1000, 5)

# 2. Filter both dataframes
df1_sub <- results_discrete_1 %>% filter(SimID %in% random_ids)
df2_sub <- results_discrete_2 %>% filter(SimID %in% random_ids)

df1_sub$Source <- "N2"
df2_sub$Source <- "N1"
combined_df <- dplyr::bind_rows(df1_sub, df2_sub)

long_df <- pivot_longer(
  combined_df,
  cols = c(N_1, N_2),
  names_to = "Population",
  values_to = "Abundance"
)

ggplot(long_df, aes(x = timestep, y = log(Abundance), color = Population)) +
  geom_line() +
  facet_grid(Source ~ SimID, scales = "free_y") +
  labs(
    title = "Dynamics of N.1 and N.2 for Random Simulations",
    x = "Timestep", y = "Abundance"
  ) +
  theme_minimal() +
  
  coord_cartesian(ylim = c(-10, 10)) +
  theme(strip.text.x = element_text(angle = 90))

ggsave("dynamics_plot.log.pdf")



check.2 = cbind(results_discrete_1,results_discrete_2)
check.1 = cbind(r_1.30,r_2.30)




r_combined = r_2 %>% mutate(Dominant_r_1 = r_1$Dominant) %>% filter((Dominant==-100)|(Dominant_r_1==-100))
r_combined %>% dim



###continuous linear




chemostat_model_invasion_growth_rate_N1 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    #N1 <- pmax(N1, 1e-6)#noisy
    dR1_dt = D*(R1_in - R1) - c11*N1*R1
    dR2_dt = D*(R2_in - R2) - c12*N1*R2
    if (mode == "substitu") {
      dN1_dt = (c11*R1*w11 + c12*R2*w12 - D) * N1 
     # dN1_dt = dN1_dt + rnorm(1, mean = 0, sd = noise_sd * N1) #noisy
    } else if (mode == "essen") {
      dN1_dt = (min(c11*R1*w11, c12*R2*w12) - D) * N1 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
    
    list(c(dR1_dt,dR2_dt,dN1_dt))
  })
}

chemostat_model_invasion_growth_rate_N2 = function(t,state,parameters,mode){
  with(as.list(c(state,parameters)),{
    
    dR1_dt = D*(R1_in - R1) -  c21*N2*R1
    dR2_dt = D*(R2_in - R2) -  c22*N2*R2
    if (mode == "substitu") {
      dN2_dt = (c21*R1*w21 + c22*R2*w22 - D) * N2 
    } else if (mode == "essen") {
      dN2_dt = (min(c21*R1*w21, c22*R2*w22) - D) * N2 
    } else {
      stop("Invalid mode! Choose 'substitu' or 'essen'.")
    }
 
    list(c(dR1_dt,dR2_dt,dN2_dt))
  })
}


set.seed(121)
n_samples <- 1000
param_df_invasion <- data.frame(
  c11 = runif(n_samples, 0.1, 1),
  w11 = runif(n_samples, 0.1, 1),
  m11 = runif(n_samples, 0, 0.1),
  c12 = runif(n_samples, 0.1, 1),
  w12 = runif(n_samples, 0.1, 1),
  m12 = runif(n_samples, 0, 0.1),
  c21 = runif(n_samples, 0.1, 1),
  w21 = runif(n_samples, 0.1, 1),
  m21 = runif(n_samples, 0, 0.1),
  c22 = runif(n_samples, 0.1, 1),
  w22 = runif(n_samples, 0.1, 1),
  m22 = runif(n_samples, 0, 0.1)
)
param_df_invasion$w11 = 1-param_df_invasion$c11
param_df_invasion$w12 = 1-param_df_invasion$c12
param_df_invasion$w21 = 1-param_df_invasion$c21
param_df_invasion$w22 = 1-param_df_invasion$c22



#N1
#noise_sd <- 0.01  # define noise strength


fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N1 = c(R1 = 10, R2 = 10, N1 = 0.1)

results.N1 <- lapply(1:nrow(param_df_invasion[,1:6]), function(i) {
  params <- c(param_df_invasion[i,1:6],fixed_params) # define noise strength
  sim <- ode(y = state_chemo_N1 , times = times_chemo, func = chemostat_model_invasion_growth_rate_N1, parms = params, mode = "substitu")
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N1 <- bind_rows(results.N1)

results_df.N1.equili = results_df.N1[results_df.N1$time==200.0,]


priority.results.invasion = data.frame()

for (i in 1:1000){
  params.2 <- param_df_invasion[i,7:12]
  params.1 <- results_df.N1.equili[i,]
   D=0.1
  
  if ((((params.2$c21)*(params.2$w21)*(params.1$R1)+(params.2$c22)*(params.2$w22)*(params.1$R2)-D)<0) ){
    priority.results.invasion[i,1] = i 
  }
  else  {
    priority.results.invasion[i,1] = NA 
  }
   
   
   if (params.1$N1>0.1) {
     priority.results.invasion[i,2] = i 
   }
   else  {
     priority.results.invasion[i,2] = NA 
   }
   
   
}


#N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,200,by = 0.1)

state_chemo_N2 = c(R1 = 10, R2 = 10, N2 = 0.1)

results.N2 <- lapply(1:nrow(param_df_invasion[,7:12]), function(i) {
  params <- c(param_df_invasion[i,7:12],fixed_params)
  sim <- ode(y = state_chemo_N2 , times = times_chemo, func = chemostat_model_invasion_growth_rate_N2, parms = params, mode = "substitu")
  sim_df <- as.data.frame(sim)
  sim_df$param_set <- i  # Track parameter set
  sim_df
  
})

results_df.N2 <- bind_rows(results.N2)

results_df.N2.equili = results_df.N2[results_df.N2$time==200.0,]


for (i in 1:1000){
  params.2 <- param_df_invasion[i,1:6]
  params.1 <- results_df.N2.equili[i,]

  
  if (((params.2$c11)*(params.2$w11)*(params.1$R1)+(params.2$c12)*(params.2$w12)*(params.1$R2)-D)<0){
    priority.results.invasion[i,3] = i 
  }
  else  {
    priority.results.invasion[i,3] = NA 
  }
  
  if (params.1$N2>0.1) {
    priority.results.invasion[i,4] = i 
  }
  else  {
    priority.results.invasion[i,4] = NA 
  }
  
}


library(tidyverse)

#######
#
priority.results.invasion %>% filter(!is.na(V1), !is.na(V2),!is.na(V3), !is.na(V4)) %>% dim



## ## ## ## 
## ## ## ## 
## simulation in discrete time linear process add noise and every day
## ## ## ## 
## ## ## ## 
## ## ## ## 
## ## ## ## 


set.seed(121)
set.seed(122)
set.seed(123)

n_samples <- 1000
param_df_invasion <- data.frame(
  c11 = runif(n_samples, 0.1, 1),
  w11 = runif(n_samples, 0.1, 1),
  m11 = runif(n_samples, 0, 0.1),
  c12 = runif(n_samples, 0.1, 1),
  w12 = runif(n_samples, 0.1, 1),
  m12 = runif(n_samples, 0, 0.1),
  c21 = runif(n_samples, 0.1, 1),
  w21 = runif(n_samples, 0.1, 1),
  m21 = runif(n_samples, 0, 0.1),
  c22 = runif(n_samples, 0.1, 1),
  w22 = runif(n_samples, 0.1, 1),
  m22 = runif(n_samples, 0, 0.1)
)
param_df_invasion$w11 = 1-param_df_invasion$c11
param_df_invasion$w12 = 1-param_df_invasion$c12
param_df_invasion$w21 = 1-param_df_invasion$c21
param_df_invasion$w22 = 1-param_df_invasion$c22



simulate_chemostat_manual_N1 <- function(params, state, times, mode = "substitu") {
  dt <- times[2] - times[1]  # timestep
  
  D <- params$D
  R1_in <- params$R1_in
  R2_in <- params$R2_in
  c11 <- params$c11
  c12 <- params$c12
  w11 <- params$w11
  w12 <- params$w12
  sigma <- 0
  R_1 <- numeric(length(times))
  R_2 <- numeric(length(times))
  N_1 <- numeric(length(times))
  
  R_1[1] <- state["R1"]
  R_2[1] <- state["R2"]
  N_1[1] <- state["N1"]
  
  for (i in 2:length(times)) {
    R1 = R_1[i-1] 
    R2 = R_2[i-1]
    N1 = N_1[i-1]
    # Resource dynamics
    
    dR1_dt <- D * (R1_in - R1) - c11 * N1 * R1
    dR2_dt <- D * (R2_in - R2) - c12 * N1 * R2
    noise <- sigma * N1 * rnorm(1, mean = 0, sd = sqrt(dt))
    # N2 growth
    if (mode == "substitu") {
      dN1_dt <- (c11 * R1 * w11 + c12 * R2 * w12 - D) * N1 + noise
    } else if (mode == "essen") {
      dN1_dt <- (min(c11 * R1 * w11, c12 * R2 * w12) - D) * N1 + noise
    } else {
      stop("Invalid mode!")
    }
    
    # Euler update
    R1_new <- R1 + dR1_dt * dt
    R2_new <- R2 + dR2_dt * dt
    N1_new <- N1 + dN1_dt * dt
    
    # Prevent negative values
    R_1[i] <- max(R1_new, 0)
    R_2[i] <- max(R2_new, 0)
    N_1[i] <- max(N1_new, 0)
    
  }
  
  data.frame(time = times, R1 = R_1, R2 = R_2, N1 = N_1)
}

simulate_chemostat_manual_N2 <- function(params, state, times, mode = "substitu") {
  dt <- times[2] - times[1]  # timestep
  
  D <- params$D
  R1_in <- params$R1_in
  R2_in <- params$R2_in
  c21 <- params$c21
  c22 <- params$c22
  w21 <- params$w21
  w22 <- params$w22
  sigma <- 0
  R_1 <- numeric(length(times))
  R_2 <- numeric(length(times))
  N_2 <- numeric(length(times))
  
  R_1[1] <- state["R1"]
  R_2[1] <- state["R2"]
  N_2[1] <- state["N2"]
  
  for (i in 2:length(times)) {
    R1 = R_1[i-1] 
    R2 = R_2[i-1]
    N2 = N_2[i-1]
      # Resource dynamics
     
      dR1_dt <- D * (R1_in - R1) - c21 * N2 * R1
      dR2_dt <- D * (R2_in - R2) - c22 * N2 * R2
      noise <- sigma * N2 * rnorm(1, mean = 0, sd = sqrt(dt))
      # N2 growth
      if (mode == "substitu") {
        dN2_dt <- (c21 * R1 * w21 + c22 * R2 * w22 - D) * N2 + noise
      } else if (mode == "essen") {
        dN2_dt <- (min(c21 * R1 * w21, c22 * R2 * w22) - D) * N2 + noise
      } else {
        stop("Invalid mode!")
      }
      
      # Euler update
      R1_new <- R1 + dR1_dt * dt
      R2_new <- R2 + dR2_dt * dt
      N2_new <- N2 + dN2_dt * dt
      
      # Prevent negative values
      R_1[i] <- max(R1_new, 0)
      R_2[i] <- max(R2_new, 0)
      N_2[i] <- max(N2_new, 0)
    
  }
  
  data.frame(time = times, R1 = R_1, R2 = R_2, N2 = N_2)
}


##N1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,100,by = 0.1)

state_chemo_N1 = c(R1 = 10, R2 = 10, N1 = 0.1)
# Loop over each row in param_df_invasion[, 1:6]
results.N1 <- lapply(1:nrow(param_df_invasion), function(i) {
  # Combine individual and fixed parameters
  params <- c(as.list(param_df_invasion[i, 1:6]), fixed_params)
  
  # Run the manual simulation
  sim_df <- simulate_chemostat_manual_N1(
    params = params,
    state = state_chemo_N1,
    times = times_chemo,
    mode = "substitu"  # or "essen"
  )
  
  sim_df$param_set <- i  # Track parameter set
  sim_df
})

# Combine all results into one data frame
all_results_df.N1 <- do.call(rbind, results.N1)

all_results_df.N1.equili = all_results_df.N1[all_results_df.N1$time==100.0,]




priority.results.invasion = data.frame()

##N2 invade N1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,5,by = 0.1)

results.N1.invasiom <- lapply(1:nrow(param_df_invasion), function(i) {
  # Combine individual and fixed parameters
  params <- c(as.list(param_df_invasion[i,]), fixed_params)
  state_chemo_N1 = c(R1 = all_results_df.N1.equili[i,2], R2 = all_results_df.N1.equili[i,3], N2 = 0.1,N1=all_results_df.N1.equili[i,4])
  # Run the manual simulation
  sim_df <- simulate_chemostat_manual_N1N2(
    params = params,
    state = state_chemo_N1,
    times = times_chemo,
    mode = "substitu"  # or "essen"
  )
  
  sim_df$param_set <- i  # Track parameter set
  sim_df
})

all_results.N1.invasiom <- do.call(rbind, results.N1.invasiom)

results.N1.invasiom.final = all_results.N1.invasiom[all_results.N1.invasiom$time==5.0,]


for (i in 1:1000){

  
  if (results.N1.invasiom.final[i,5]<0.1) {
    priority.results.invasion[i,1] = 1
  }
  else  {
    priority.results.invasion[i,1] = 0 
  }
  
  
  
}




##N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,100,by = 0.1)

state_chemo_N2 = c(R1 = 10, R2 = 10, N2 = 0.1)
# Loop over each row in param_df_invasion[, 7:12]
results.N2 <- lapply(1:nrow(param_df_invasion), function(i) {
  # Combine individual and fixed parameters
  params <- c(as.list(param_df_invasion[i, 7:12]), fixed_params)
  
  # Run the manual simulation
  sim_df <- simulate_chemostat_manual_N2(
    params = params,
    state = state_chemo_N2,
    times = times_chemo,
    mode = "substitu"  # or "essen"
  )
  
  sim_df$param_set <- i  # Track parameter set
  sim_df
})

# Combine all results into one data frame
all_results_df.N2 <- do.call(rbind, results.N2)

all_results_df.N2.equili = all_results_df.N2[all_results_df.N2$time==100.0,]




##N1 invade N2
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,5,by = 0.1)

results.N2.invasiom <- lapply(1:nrow(param_df_invasion), function(i) {
  # Combine individual and fixed parameters
  params <- c(as.list(param_df_invasion[i,]), fixed_params)
  state_chemo_N2 = c(R1 = all_results_df.N2.equili[i,2], R2 = all_results_df.N2.equili[i,3], N1 = 0.1,N2=all_results_df.N2.equili[i,4])
  # Run the manual simulation
  sim_df <- simulate_chemostat_manual_N1N2(
    params = params,
    state = state_chemo_N2,
    times = times_chemo,
    mode = "substitu"  # or "essen"
  )
  
  sim_df$param_set <- i  # Track parameter set
  sim_df
})

all_results.N2.invasiom <- do.call(rbind, results.N2.invasiom)

results.N2.invasiom.final = all_results.N2.invasiom[all_results.N2.invasiom$time==5.0,]


for (i in 1:1000){
  if (results.N2.invasiom.final[i,4]<0.1) {
    priority.results.invasion[i,2] = 1
  }
  else  {
    priority.results.invasion[i,2] = 0 
  }
}


priority.results.invasion %>% filter(V1==1, V2==1) %>% dim







########################################
####NEW codes of simulation with noise##
########################################

# 1. Define parameters and helper functions
library(tidyverse)
set.seed(123)

# Chemostat simulation function with mortality and noise
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

# 2. Generate parameter sets
n_sets <- 1000  # Number of parameter sets
param_df <- tibble(
  c11 = runif(n_sets, 0.1, 1),
  c12 = runif(n_sets, 0.1, 1),
  w11 = runif(n_sets, 0.1, 1),
  w12 = runif(n_sets, 0.1, 1),
  #m1 = runif(n_sets, 0.01, 0.1),  # Species 1 mortality
  c21 = runif(n_sets, 0.1, 1),
  c22 = runif(n_sets, 0.1, 1),
  w21 = runif(n_sets, 0.1, 1),
  w22 = runif(n_sets, 0.1, 1),
  #m2 = runif(n_sets, 0.01, 0.1),  # Species 2 mortality
  #sigma = runif(n_sets, 0.1, 0.5) # Noise magnitude
  sigma = 0.1
)

# Add fixed parameters
param_df <- param_df %>% mutate(
  D = 0.1,
  R1_in = 10,
  R2_in = 10,
  param_set = row_number(),

)

# 3. Mutual invasion simulation
results <- map_dfr(1:nrow(param_df), function(i) {
  params <- as.list(param_df[i, ])
  
  # Phase 1: Establish resident (N1)
  resident <- simulate_chemostat(
    params,
    state = c(R1 = 10, R2 = 10, N1 = 0.1, N2 = 0),
    times = seq(0, 500, 0.1),
    mode = "substitu"
  )
  
  # Get equilibrium state
  eq_state <- tail(resident, 1) %>% select(R1, R2, N1, N2) %>% unlist()
  
  # Phase 2: Invasion with noise (N2 invades)
  invader_state <- eq_state
  invader_state["N2"] <- 0.01  # Introduce invader
  
  invasion <- simulate_chemostat(
    params,
    state = invader_state,
    times = seq(0, 50, 0.1),
    mode = "substitu"
  )
  
  # Get final state
  final <- tail(invasion, 1)
  
  # Phase 3: Reverse invasion (N1 invades N2)
  # First establish N2 resident
  resident2 <- simulate_chemostat(
    params,
    state = c(R1 = 10, R2 = 10, N1 = 0, N2 = 0.1),
    times = seq(0, 500, 0.1),
    mode = "substitu"
  )
  
  eq_state2 <- tail(resident2, 1) %>% select(R1, R2, N1, N2) %>% unlist()
  
  # N1 invasion
  invader_state2 <- eq_state2
  invader_state2["N1"] <- 0.01
  
  invasion2 <- simulate_chemostat(
    params,
    state = invader_state2,
    times = seq(0, 50, 0.1),
    mode = "substitu"
  )
  
  final2 <- tail(invasion2, 1)
  
  # Return results
  tibble(
    param_set = i,
    # Invasion 1: N2 invading N1 resident
    N1_final_1 = final$N1,
    N2_final_1 = final$N2,
    invasion_success_1 = final$N2 > 0.01,  # Invasion threshold
    
    # Invasion 2: N1 invading N2 resident
    N1_final_2 = final2$N1,
    N2_final_2 = final2$N2,
    invasion_success_2 = final2$N1 > 0.01,
    
    # Priority effect classification
    priority_effect = !invasion_success_1 & !invasion_success_2
  )
})


# 4. Analyze results
priority_effect_summary <- results %>%
  group_by(priority_effect) %>%
  summarise(count = n(), percentage = n()/n_sets * 100)
priority_effect_summary


# 5. Visualize outcomes
# Plot priority effect cases
priority_cases <- results %>% filter(priority_effect)

if (nrow(priority_cases)) {
  ggplot(priority_cases, aes(x = N1_final_1, y = N2_final_1)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0.001, linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0.001, linetype = "dashed", color = "red") +
    labs(title = "Priority Effect Outcomes",
         subtitle = "Each species resists invasion by the other",
         x = "Resident (N1) Final Biomass", 
         y = "Invader (N2) Final Biomass") +
    theme_bw()
} else {
  message("No priority effects detected in this parameter set")
}

# 6. Full simulation example (for plotting)
example_params <- as.list(param_df[1, ])
example_sim <- simulate_chemostat(
  example_params,
  state = c(R1 = 10, R2 = 10, N1 = 0.1, N2 = 0.01),
  times = seq(0, 200, 0.5),
  mode = "substitu"
)

# Plot dynamics
example_sim_long <- example_sim %>%
  pivot_longer(cols = c(N1, N2), 
               names_to = "Species", 
               values_to = "Biomass")

ggplot(example_sim_long, aes(x = time, y = Biomass, color = Species)) +
  geom_line() +
  scale_y_log10() +
  labs(title = "Example Invasion Dynamics",
       subtitle = "Species N2 invading resident N1",
       x = "Time", y = "Biomass (log scale)") +
  theme_bw()




###
###
###
filtered_row_names <- priority.results.invasion %>%
  tibble::rownames_to_column("original_row") %>%  # Convert rownames to a column
  filter(V1==1, V2==1) %>%
  pull(original_row)

filtered_row_names
filtered_indices <- as.integer(filtered_row_names)

filtered_paras = param_df_invasion[filtered_indices,]

set.seed(123)
n_sets <- 1000  # Number of parameter sets
param_df_invasion <- tibble(
  c11 = runif(n_sets, 0.1, 1),
  c12 = runif(n_sets, 0.1, 1),
  w11 = runif(n_sets, 0.1, 1),
  w12 = runif(n_sets, 0.1, 1),
  #m1 = runif(n_sets, 0.01, 0.1),  # Species 1 mortality
  c21 = runif(n_sets, 0.1, 1),
  c22 = runif(n_sets, 0.1, 1),
  w21 = runif(n_sets, 0.1, 1),
  w22 = runif(n_sets, 0.1, 1),
  #m2 = runif(n_sets, 0.01, 0.1),  # Species 2 mortality
  #sigma = runif(n_sets, 0.1, 0.5) # Noise magnitude
  #sigma = 0
)

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

###1111
##N1 0.1 N2 10
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,500,by = 0.1)

state_chemo_111 = c(R1 = 10, R2 = 10,N1 = .1, N2 = 10)
# Loop over each row in param_df_invasion[, 7:12]
results.111 <- lapply(1:nrow(param_df_invasion), function(i) {
  # Combine individual and fixed parameters
  params <- c(as.list(param_df_invasion[i, ]), fixed_params)
  
  # Run the manual simulation
  sim_df <- simulate_chemostat_manual_N1N2(
    params = params,
    state = state_chemo_111,
    times = times_chemo,
    mode = "substitu"  # or "essen"
  )
  
  sim_df$param_set <- i  # Track parameter set
  sim_df
})

# Combine all results into one data frame
all_results_df.111 <- do.call(rbind, results.111)

all_results_df.111.equili = all_results_df.111[all_results_df.111$time==500.0,]

###222
##N1 10 N2 0.1
fixed_params = c(D = 0.1, R1_in = 10,R2_in = 10)

times_chemo = seq(0,500,by = 0.1)

state_chemo_222 = c(R1 = 10, R2 = 10,N1 = 10, N2 = .1)
# Loop over each row in param_df_invasion[, 7:12]
results.222 <- lapply(1:nrow(param_df_invasion), function(i) {
  # Combine individual and fixed parameters
  params <- c(as.list(param_df_invasion[i, ]), fixed_params)
  
  # Run the manual simulation
  sim_df <- simulate_chemostat_manual_N1N2(
    params = params,
    state = state_chemo_222,
    times = times_chemo,
    mode = "substitu"  # or "essen"
  )
  
  sim_df$param_set <- i  # Track parameter set
  sim_df
})

# Combine all results into one data frame
all_results_df.222 <- do.call(rbind, results.222)

all_results_df.222.equili = all_results_df.222[all_results_df.222$time==500.0,]


####summary


all_results.equili = cbind(all_results_df.111.equili[,2:6],all_results_df.222.equili[2:5])

names(all_results.equili) = c("R1_2","R2_2",'N1_2','N2_2',"para_set","R1_1","R2_1",'N1_1','N2_1')

df_all_results.equili <- all_results.equili %>%
  mutate(diff_N1 = abs(N1_2 - N1_1),
         diff_N2 = abs(N2_2 - N2_1))            %>%
  mutate(
    diff_N1_per = ifelse(diff_N1 != 0, abs(N1_2 - N1_1) / (N1_2 + N1_1), 0),
    diff_N2_per = ifelse(diff_N2 != 0, abs(N2_2 - N2_1) / (N2_2 + N2_1), 0)
  ) %>%
  mutate(
    results_1 = case_when(
      (N1_1 > 100*N2_1)  ~ 'N1 outcompete',
      (N2_1 > 100*N1_1)  ~ 'N2 outcompete',
      TRUE ~ 'coexsitence'),
    results_2 = case_when(
      (N1_2 > 100*N2_2)  ~ 'N1 outcompete',
      (N2_2 > 100*N1_2)  ~ 'N2 outcompete',
      TRUE ~ 'coexsitence'
)
  )

df_all_results.equili_new = df_all_results.equili[,c("para_set","N1_1","N1_2","diff_N1","diff_N1_per","N2_1","N2_2","diff_N2","diff_N2_per",'results_1','results_2')]

df_all_results.equili_new %>%
  group_by(results_1, results_2) %>%
  summarise(count = n(), .groups = 'drop')


filtered_row_names <- df_all_results.equili_new %>%
 filter(results_1=='N1 outcompete',results_2=='coexsitence',N2_2>5) %>%
 #  filter(results_2=='N2 outcompete',results_1=='coexsitence',N1_1>5) %>%
  pull(para_set)








library(dplyr)
library(ggplot2)
library(patchwork)  # for side-by-side plots

# Assume both dataframes have a column `param_set` to identify the set,
# and columns: time, R1, R2, N1, etc.

# 1. Select 20 random param_set IDs
set.seed(1234)  # for reproducibilit

random_ids <- sample(1:1000, 5)
filtered_row_names
random_ids = c( 577 ,715 ,809 ,811 ,815)

# 2. Filter both dataframes
df1_sub <- all_results_df.111 %>% filter(param_set %in% random_ids)
df2_sub <- all_results_df.222 %>% filter(param_set %in% random_ids)

df1_sub$Source <- "N2_first"
df2_sub$Source <- "N1_first"
combined_df <- dplyr::bind_rows(df1_sub, df2_sub)

long_df <- pivot_longer(
  combined_df,
  cols = c(N1, N2),
  names_to = "Population",
  values_to = "Abundance"
)


ggplot(long_df, aes(x = time, y = Abundance, color = Population)) +
  geom_line() +
  facet_grid(Source ~ param_set, scales = "free_y") +
  labs(
    title = "Dynamics of N.1 and N.2 for Random Simulations",
    x = "Timestep", y = "Abundance (log scale)"
  ) +
  theme_minimal() +
  scale_y_log10() +
  theme(strip.text.x = element_text(angle = 90))

ggsave("1.pdf")










########################################
####ESSENTIAL resources##
########################################
##
##
##
##
## linear
results.ess <- map_dfr(1:nrow(param_df), function(i) {
  params <- as.list(param_df[i, ])
  
  # Phase 1: Establish resident (N1)
  resident <- simulate_chemostat(
    params,
    state = c(R1 = 10, R2 = 10, N1 = 0.1, N2 = 0),
    times = seq(0, 100, 0.1),
    mode = "essen"
  )
  
  # Get equilibrium state
  eq_state <- tail(resident, 1) %>% select(R1, R2, N1, N2) %>% unlist()
  
  # Phase 2: Invasion with noise (N2 invades)
  invader_state <- eq_state
  invader_state["N2"] <- 0.01  # Introduce invader
  
  invasion <- simulate_chemostat(
    params,
    state = invader_state,
    times = seq(0, 50, 0.1),
    mode = "essen"
  )
  
  # Get final state
  final <- tail(invasion, 1)
  
  # Phase 3: Reverse invasion (N1 invades N2)
  # First establish N2 resident
  resident2 <- simulate_chemostat(
    params,
    state = c(R1 = 10, R2 = 10, N1 = 0, N2 = 0.1),
    times = seq(0, 100, 0.1),
    mode = "essen"
  )
  
  eq_state2 <- tail(resident2, 1) %>% select(R1, R2, N1, N2) %>% unlist()
  
  # N1 invasion
  invader_state2 <- eq_state2
  invader_state2["N1"] <- 0.01
  
  invasion2 <- simulate_chemostat(
    params,
    state = invader_state2,
    times = seq(0, 50, 0.1),
    mode = "essen"
  )
  
  final2 <- tail(invasion2, 1)
  
  # Return results
  tibble(
    param_set = i,
    # Invasion 1: N2 invading N1 resident
    N1_final_1 = final$N1,
    N2_final_1 = final$N2,
    invasion_success_1 = final$N2 > 0.01,  # Invasion threshold
    
    # Invasion 2: N1 invading N2 resident
    N1_final_2 = final2$N1,
    N2_final_2 = final2$N2,
    invasion_success_2 = final2$N1 > 0.01,
    
    # Priority effect classification
    priority_effect = !invasion_success_1 & !invasion_success_2
  )
})


# 4. Analyze results
priority_effect_summary.ess <- results.ess %>%
  group_by(priority_effect) %>%
  summarise(count = n(), percentage = n()/n_sets * 100)
priority_effect_summary.ess







###nonlinear 
source("chemostat_functions.R")
# Generate parameter sets WITH K VALUES
n_sets <- 1000
param_df_nonlinear <- tibble(
  # Species 1 parameters
  c11 = runif(n_sets, 0.1, 1),
  c12 = runif(n_sets, 0.1, 1),
  w11 = runif(n_sets, 0.1, 1),
  w12 = runif(n_sets, 0.1, 1),
  K11 = runif(n_sets, 0.1, 5),  
  K12 = runif(n_sets, 0.1, 5),  
  
  # Species 2 parameters
  c21 = runif(n_sets, 0.1, 1),
  c22 = runif(n_sets, 0.1, 1),
  w21 = runif(n_sets, 0.1, 1),
  w22 = runif(n_sets, 0.1, 1),
  K21 = runif(n_sets, 0.1, 5),  
  K22 = runif(n_sets, 0.1, 5),  
  
  sigma = 0.1
) %>% mutate(
  D = 0.1,
  R1_in = 50,
  R2_in = 50,
  param_set = row_number()
)

# Mutual invasion simulation (UNCHANGED)
results.ess <- map_dfr(1:nrow(param_df_nonlinear), function(i) {

    params <- as.list(param_df_nonlinear[i, ])
    
    # Phase 1: Establish resident (N1)
    resident <- simulate_chemostat_nonlinear(
      params,
      state = c(R1 = 50, R2 = 50, N1 = 0.1, N2 = 0),
      times = seq(0, 100, 0.1),
      mode = "essen"
    )
    
    # Get equilibrium state
    eq_state <- tail(resident, 1) %>% select(R1, R2, N1, N2) %>% unlist()
    
    # Phase 2: Invasion with noise (N2 invades)
    invader_state <- eq_state
    invader_state["N2"] <- 0.01  # Introduce invader
    
    invasion <- simulate_chemostat_nonlinear(
      params,
      state = invader_state,
      times = seq(0, 50, 0.1),
      mode = "essen"
    )
    
    # Get final state
    final <- tail(invasion, 1)
    
    # Phase 3: Reverse invasion (N1 invades N2)
    # First establish N2 resident
    resident2 <- simulate_chemostat_nonlinear(
      params,
      state = c(R1 = 50, R2 = 50, N1 = 0, N2 = 0.1),
      times = seq(0, 100, 0.1),
      mode = "essen"
    )
    
    eq_state2 <- tail(resident2, 1) %>% select(R1, R2, N1, N2) %>% unlist()
    
    # N1 invasion
    invader_state2 <- eq_state2
    invader_state2["N1"] <- 0.01
    
    invasion2 <- simulate_chemostat_nonlinear(
      params,
      state = invader_state2,
      times = seq(0, 50, 0.1),
      mode = "essen"
    )
    
    final2 <- tail(invasion2, 1)
    
    # Return results
    tibble(
      param_set = i,
      # Invasion 1: N2 invading N1 resident
      N1_final_1 = final$N1,
      N2_final_1 = final$N2,
      invasion_success_1 = final$N2 > 0.01,  # Invasion threshold
      
      # Invasion 2: N1 invading N2 resident
      N1_final_2 = final2$N1,
      N2_final_2 = final2$N2,
      invasion_success_2 = final2$N1 > 0.01,
      
      # Priority effect classification
      priority_effect = !invasion_success_1 & !invasion_success_2
    )
  })
  
  
# 4. Analyze results
priority_effect_summary <- results.ess %>%
    group_by(priority_effect) %>%
    summarise(count = n(), percentage = n()/n_sets * 100)
priority_effect_summary
  
  
###################################
#######only ONE resource#######
##############################
############################
n_sets <- 1000  # Number of parameter sets
param_df_1R <- tibble(
  c11 = runif(n_sets, 0.1, 1),
  c12 = runif(n_sets, 0.1, 1),
  w11 = runif(n_sets, 0.1, 1),
  w12 = runif(n_sets, 0.1, 1),
  #m1 = runif(n_sets, 0.01, 0.1),  # Species 1 mortality
  c21 = runif(n_sets, 0.1, 1),
  c22 = runif(n_sets, 0.1, 1),
  w21 = runif(n_sets, 0.1, 1),
  w22 = runif(n_sets, 0.1, 1),
  #m2 = runif(n_sets, 0.01, 0.1),  # Species 2 mortality
  sigma = 0.1 # Noise magnitude
  # sigma = 0
)

# Add fixed parameters
param_df_1R <- param_df_1R %>% mutate(
  D = 0.1,
  R1_in = 0,
  R2_in = 10,
  param_set = row_number(),
  # sigma = 0
)


results.1R <- map_dfr(1:nrow(param_df_1R), function(i) {
  params <- as.list(param_df_1R[i, ])
  
  # Phase 1: Establish resident (N1)
  resident <- simulate_chemostat(
    params,
    state = c(R1 = 0, R2 = 10, N1 = 0.1, N2 = 0),
    times = seq(0, 500, 0.1),
    mode = "substitu"
  )
  
  # Get equilibrium state
  eq_state <- tail(resident, 1) %>% select(R1, R2, N1, N2) %>% unlist()
  
  # Phase 2: Invasion with noise (N2 invades)
  invader_state <- eq_state
  invader_state["N2"] <- 0.01  # Introduce invader
  
  invasion <- simulate_chemostat(
    params,
    state = invader_state,
    times = seq(0, 50, 0.1),
    mode = "substitu"
  )
  
  # Get final state
  final <- tail(invasion, 1)
  
  # Phase 3: Reverse invasion (N1 invades N2)
  # First establish N2 resident
  resident2 <- simulate_chemostat(
    params,
    state = c(R1 = 0, R2 = 10, N1 = 0, N2 = 0.1),
    times = seq(0, 500, 0.1),
    mode = "substitu"
  )
  
  eq_state2 <- tail(resident2, 1) %>% select(R1, R2, N1, N2) %>% unlist()
  
  # N1 invasion
  invader_state2 <- eq_state2
  invader_state2["N1"] <- 0.01
  
  invasion2 <- simulate_chemostat(
    params,
    state = invader_state2,
    times = seq(0, 50, 0.1),
    mode = "substitu"
  )
  
  final2 <- tail(invasion2, 1)
  
  # Return results
  tibble(
    param_set = i,
    # Invasion 1: N2 invading N1 resident
    N1_final_1 = final$N1,
    N2_final_1 = final$N2,
    invasion_success_1 = final$N2 > 0.01,  # Invasion threshold
    
    # Invasion 2: N1 invading N2 resident
    N1_final_2 = final2$N1,
    N2_final_2 = final2$N2,
    invasion_success_2 = final2$N1 > 0.01,
    
    # Priority effect classification
    priority_effect = !invasion_success_1 & !invasion_success_2
  )
})


# 4. Analyze results
priority_effect_summary.1R <- results.1R %>%
  group_by(priority_effect) %>%
  summarise(count = n(), percentage = n()/n_sets * 100)
priority_effect_summary.1R





f <- function(x) {
  0.1*x / (x + 0.1)
}

# Create a sequence of x values (avoid x = -2 to prevent division by zero)
x_vals <- seq(-10, 10, by = 0.1)
x_vals <- x_vals[x_vals != -2]  # remove the discontinuity at x = -2

# Calculate corresponding y values
y_vals <- f(x_vals)

# Plot the function
plot(x_vals, y_vals, type = "l", lwd = 2, col = "blue",
     xlab = "x", ylab = "y",
     main = expression(y == frac(x, x+2)))

# Add a vertical dashed line to indicate the asymptote at x = -2
abline(v = -2, col = "red", lty = 2)




#####
####
####
######
#####
####
####
######
library(Sim.DiffProc)

# Define parameters
r <- 1        # growth rate
K <- 100000      # carrying capacity
sigma <- 0.8  # noise intensity
N0 <- 100      # initial population

# Define drift and diffusion terms
drift <- expression(r * N * (1 - N/K))
diffusion <- expression(sigma * N)

# Simulate the stochastic logistic model
set.seed(123)
sde.model <- snssde1d(drift = drift,
                      diffusion = diffusion,
                      x0 = N0,
                      M = 1,       # number of simulations
                      N = 1000,    # time steps
                      T = 20,      # total time
                      method = "euler")

plot(sde.model, main = "Stochastic Logistic Growth with Noise", ylab = "Population N")




### the experiment

##criteria
d = d_comp_cl %>% 
      filter(n_sp == 2) %>% 
      filter(medium != 'NA') %>%
      group_by(lst_sp, high, n_sp, medium, Well.Name, day) %>% 
      mutate(freq = n/sum(n)) 
d_last = d %>% filter(day == 49)

d_last_fre = d_last %>% group_by(comm,medium)%>%
       mutate(mean_fre = mean(freq, na.rm = TRUE))%>%  group_by(lst_sp, medium) %>%
       mutate(p = ifelse(all(mean_fre > 0.7), 1, 0)) %>% mutate(group = paste(lst_sp,'_',medium)) %>% filter(p==1)
unique(d_last_fre$group)



#logistic reregression


library(dplyr)


d_bin =  d_last %>% filter(lst_sp %in% c("1_2","2_3"))%>% mutate(group = paste0(lst_sp,medium))%>%
 mutate(p = ifelse(group %in% c("1_20","1_22","1_23","2_30"), 1, 0))%>%
  mutate(numr = ifelse(medium==1|medium==0, 1, 2))
 
d_bin$numr <- factor(d_bin$numr)  # make sure it's a factor
d_bin$medium <- factor(d_bin$medium)

model1 <- glm(p ~ numr, data = d_bin, family = binomial)
summary(model1)
model2 <- glm(p ~ medium, data = d_bin, family = binomial)
summary(model2)

ggplot(d_bin, aes(x = numr, y = p)) +
  geom_jitter(width = 0.1, height = 0.05, size = 2, alpha = 0.6) +
  labs(x = "Num of Resource", y = "Outcome (priority effect)", title = "The effect of Res.Num on the outcomes") +
  theme_minimal()

ggplot(d_bin, aes(x = medium, y = p)) +
  geom_jitter(width = 0.1, height = 0.05, size = 2, alpha = 0.6) +
  labs(x = "Medium", y = "Outcome (priority effect)", title = "The effect of medium on the outcomes") +
  theme_minimal()

ggsave(filename = paste0("The effect of Res.Num on the outcomes_2",  ".pdf"),  width = 8, height = 6)


table(d_bin$p)
summary(d_bin$numr)
