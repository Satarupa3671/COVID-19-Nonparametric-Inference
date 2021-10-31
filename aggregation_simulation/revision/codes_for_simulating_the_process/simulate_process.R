## aggregation simulation -- simulating the process

rm(list=ls())
library(tidyverse)
library(parallel)
load("./set1_wd_30.RData")
source("./func_101221.R")
source('./fitting_function_101620.R')

#### simulation study (eq (1)-(6))

## initial parameter:

Time <- 100

true.para <- list(gamma.true = 0.03,
                  rhoA.true = 0.06,
                  rhoH.t.true = c(rep(0.03,20),rep(0.04,40),rep(0.05,40)), 
                  alpha.true = 0.2,
                  delta.t.true = c(rep(0.007,50),rep(0.005,50)),
                  phi.t.true = 0.001+1:Time*0.00005)

covid_simulation <- function(dat.int, true.para,Time){
  
  A <- c(dat.int$At,rep(0,Time-1))
  S <- c(dat.int$St,rep(0,Time-1))
  R <- c(dat.int$Rt,rep(0,Time-1))
  Q <- c(dat.int$Qt,rep(0,Time-1))
  H <- c(dat.int$Ht,rep(0,Time-1))
  D <- c(dat.int$Dt,rep(0,Time-1))
  C <- c(dat.int$Ct,rep(0,Time-1))
  Test <- dat.int$Tt
  kappa <- dat.int$kappa
  Max.pop <- dat.int$Max.pop
  
  gamma <- true.para$gamma.true
  rhoA <- true.para$rhoA.true
  rhoH.t <- true.para$rhoH.t.true
  alpha <- true.para$alpha.true
  delta.t <- true.para$delta.t.true
  phi.t <- true.para$phi.t.true
  
  
  theta.t <- NULL
  
  for (t in 1:(Time-1)){
    
    theta.t[t] <- phi.t[t]*diff(Test)[t]/H[t] ## eq(7)
    
    ## 4 sub-poisson processes in eq(2) 
    ## delA can be negative
    ## (see notes 'remarks_of_simulation.pdf')
    
    P1 <- rpois(1, theta.t[t]*A[t])
    
    P2 <- rpois(1, gamma*A[t])
    
    P3 <- rpois(1, rhoA*A[t])
    
    P4 <- rpois(1, (S[t]/(S[t]+R[t]+A[t]))*alpha*kappa[t]^2*A[t])
    
    delA.t <- -(P1+P2+P3)+P4
    
    ## 2 sub-poisson processes in eq(3)
    
    P5 <- rpois(1, gamma*Q[t])
    
    P6 <- rpois(1, rhoA*Q[t])
    
    delQ.t <- P1-P5-P6
    
    ## 2 sub-poisson processes in eq(4)
    
    P7 <- rpois(1, rhoH.t[t]*H[t])
    
    P8 <- rpois(1, delta.t[t]*H[t])
    
    delH.t <- P2+P5-P7-P8
    
    ## eq(5) is deterministics by eq(4)
    
    delD.t <- P8
    
    ## eq(6) determined by eq(2)
    delC.t <- P1+P2
    
    ## eq(1)
    delS.t <- -delA.t
    
    ## eq(8)-(10) under assumption A5: rhoQ=rhoA
    delR.t <- P6+P7#P6+P7
    
    A[t+1] <- ifelse(A[t]+delA.t>=1,pmin(A[t]+delA.t, Max.pop), pmax(A[t]+delA.t, 1)) ## At not drop to 0
    S[t+1] <- ifelse(S[t]+delS.t>=0,pmin(S[t]+delS.t, Max.pop), pmax(S[t]+delS.t, 0))
    R[t+1] <- ifelse(R[t]+delR.t>=0,pmin(R[t]+delR.t, Max.pop), pmax(R[t]+delR.t, 0))
    Q[t+1] <- ifelse(Q[t]+delQ.t>=0,pmin(Q[t]+delQ.t, Max.pop), pmax(Q[t]+delQ.t, 0))
    H[t+1] <- ifelse(H[t]+delH.t>=0,pmin(H[t]+delH.t, Max.pop), pmax(H[t]+delH.t, 0))
    D[t+1] <- ifelse(D[t]+delD.t>=0,pmin(D[t]+delD.t, Max.pop), pmax(D[t]+delD.t, 0))
    C[t+1] <- ifelse(C[t]+delC.t>=0,pmin(C[t]+delC.t, Max.pop), pmax(C[t]+delC.t, 0))
    
    
  }
  
  dat.sim <- data.frame(At=A) %>%
    mutate(St = S,Rt = R,
           Qt = Q,Ht = H,
           Dt = D,Ct = C,
           Tt = Test[1:Time],
           kappa = kappa[1:Time])
  
  return(list(data=dat.sim,theta=theta.t))
  
}


## to gumbel

seed1 <- 490
seed2 <- 372

set.seed(10+seed1)
ten.initial <- data.frame(At = sample(100:300,10,replace = T),
                          St = sample((3*10^6):(5*10^6),10,replace = T),
                          Rt = sample(1:10,10,replace = T),
                          Qt = sample(50:100,10,replace = T),
                          Ht = sample(1:15,10,replace = T),
                          Dt = sample(0:7,10,replace = T),
                          Ct = sample(110:180,10,replace = T)) %>%
  mutate(Max.pop = At+St+Rt+Qt+Ht+Dt+Ct)


ten.Test <- sapply(1:10,function(i) dat.sim[[i]]$data$Tt)
ten.kappa <- sapply(1:10,function(i) dat.sim[[i]]$data$kappa)

dat.int.list <- lapply(1:10,function(i){
  
  return(list(At=ten.initial$At[i],St = ten.initial$St[i],
              Rt = ten.initial$Rt[i],Qt = ten.initial$Qt[i],
              Ht = ten.initial$Ht[i],Dt = ten.initial$Dt[i],
              Ct = ten.initial$Ct[i],Max.pop = ten.initial$Max.pop[i],
              Tt = ten.Test[,i],
              kappa = ten.kappa[,i])
  )
})

set.seed(200+seed2)
dat.sim.new <- lapply(1:10, function(i){
  
  return(covid_simulation(dat.int=dat.int.list[[i]],true.para,Time=Time))
})
