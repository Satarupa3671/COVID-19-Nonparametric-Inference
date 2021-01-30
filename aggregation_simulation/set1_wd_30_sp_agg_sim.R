rm(list=ls())
library(tidyverse)
library(parallel)
#### simulation study (eq (1)-(6))

## initial parameter:

Time <- 100

true.para <- list(gamma.true = 0.03,
                  rhoA.true = 0.06,#0.02 (relative smaller NIR but larger than 0.05)
                  rhoH.t.true = c(rep(0.03,20),rep(0.04,40),rep(0.05,40)), 
                  alpha.true = 0.2,
                  delta.t.true = c(rep(0.007,50),rep(0.005,50)),
                  phi.t.true = 0.001+1:Time*0.00005)

#set.seed(123)
ten.initial <- data.frame(At = sample(100:300,10,replace = T),
                          St = sample((3*10^6):(5*10^6),10,replace = T),
                          Rt = sample(1:10,10,replace = T),
                          Qt = sample(50:100,10,replace = T),
                          Ht = sample(1:15,10,replace = T),
                          Dt = sample(0:7,10,replace = T),
                          Ct = sample(110:180,10,replace = T)) %>%
  mutate(Max.pop = At+St+Rt+Qt+Ht+Dt+Ct)

ten.Test <- data.frame(k1 = cumsum(c(sample(1:300,20),
                                     sample(1000:3000,40),
                                     sample(4000:8000,30),
                                     sample(7500:12000,10))),
                       k2 = cumsum(c(sample(1:200,20),
                                     sample(400:1000,25),
                                     sample(2000:6000,30),
                                     sample(3000:6000,25))),
                       k3 = cumsum(c(sample(100:300,10),
                                     sample(1200:6000,40),
                                     sample(3000:5000,30),
                                     sample(5000:7000,20))),
                       k4 = cumsum(c(sample(50:300,20),
                                     sample(2000:4000,40),
                                     sample(2000:2500,30),
                                     sample(5000:8000,10))),
                       k5 = cumsum(c(sample(40:100,10),
                                     sample(300:1300,20),
                                     sample(2000:5000,40),
                                     sample(3000:4000,30))),
                       k6 = cumsum(c(sample(1:100,20),
                                     sample(500:3000,40),
                                     sample(2000:5000,10),
                                     sample(2500:4000,30))),
                       k7 = cumsum(c(sample(30:200,20),
                                     sample(5000:1500,30),
                                     sample(900:4000,40),
                                     sample(3000:5000,10))),
                       k8 = cumsum(c(sample(10:200,10),
                                     sample(500:1500,60),
                                     sample(2000:4000,10),
                                     sample(3000:5000,20))),
                       k9 = cumsum(c(sample(50:300,20),
                                     sample(1200:5000,60),
                                     sample(2500:4000,20))),
                       k10 = cumsum(c(sample(10:400,10),
                                      sample(1000:3500,40),
                                      sample(2000:5000,30),
                                      sample(3000:4000,20))))


ten.kappa <- data.frame(matrix(rep(lowess(c(rep(1.2,15),rep(0.8,10),rep(0.4,20),
                                            rep(0.6,15),rep(0.7,15),
                                            rep(0.9,15),rep(0.7,10)),
                                          f=1/16)$y,10),ncol=10))
#View(ten.kappa)


dat.int.list <- lapply(1:10,function(i){
  
  return(list(At=ten.initial$At[i],St = ten.initial$St[i],
              Rt = ten.initial$Rt[i],Qt = ten.initial$Qt[i],
              Ht = ten.initial$Ht[i],Dt = ten.initial$Dt[i],
              Ct = ten.initial$Ct[i],Max.pop = ten.initial$Max.pop[i],
              Tt = ten.Test[,i],
              kappa = ten.kappa[,i])
  )
})


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



source('./est_plot_func_101520.R')
source('./fitting_function_101620.R')



### obtain sampling distribuion for M=1000
ncore <- 2
sim.BS.N = 2
time <- proc.time()


est.BS.list <- mclapply(1:sim.BS.N, function(i){
  
  dat.sim <- lapply(1:10, function(i){
  
  return(covid_simulation(dat.int.list[[i]],true.para,Time=Time))
})
  
  ##Variables: all will become a Time * 10 matrix 
  delC = sapply(1:length(dat.sim),function(i){diff(dat.sim[[i]]$data$Ct)})
  delR = sapply(1:length(dat.sim),function(i){diff(dat.sim[[i]]$data$Rt)})
  delH = sapply(1:length(dat.sim),function(i){diff(dat.sim[[i]]$data$Ht)})
  delD = sapply(1:length(dat.sim),function(i){diff(dat.sim[[i]]$data$Dt)})
  Q = sapply(1:length(dat.sim),function(i){dat.sim[[i]]$data$Qt})
  H = sapply(1:length(dat.sim),function(i){dat.sim[[i]]$data$Ht})
  delTests = sapply(1:length(dat.sim),function(i){
    lowess((diff(dat.sim[[i]]$data$Tt)), f = 1/16)$y})
  TtoH_deno = sapply(1:length(dat.sim),function(i){
    lowess((dat.sim[[i]]$data$Ht[-dim(dat.sim[[i]]$data)[1]]), f = 1/16)$y})
  TtoH = delTests/TtoH_deno
  
  var = list(delC = delC, delR = delR, delH = delH, 
             delD = delD, Q = Q, H = H, TtoH = TtoH)
  
  
  
  TimeBlock = NULL
  numBlock = 9#8
  wd.length =  30#40
  n = Time
  shift = floor(((n-1) - wd.length)/(numBlock - 1))
  #
  for(j in 1:numBlock){
    TimeBlock[[j]] = (1 + shift * (j-1)) : min(shift * (j-1) + wd.length,(n-1))
  }
  
  
  gamma.grid = seq(0.001,0.05,0.001)# for set2:seq(0.001,0.08,0.001)
  rhoA.grid = seq(0.002,0.3,0.002)
  eval.mat = matrix(0,length(gamma.grid),length(rhoA.grid))
  
  for(i in 1:length(gamma.grid)){
    for(j in 1:length(rhoA.grid)){
      eval.mat[i,j]= objective.grA(var,c(gamma.grid[i],rhoA.grid[j]),TimeBlock)
    }
  }
  
  eval.min.row = apply(eval.mat,1,min)
  eval.min.col = apply(eval.mat,2,min)
  
  g.index = which(eval.min.row == min(eval.min.row )) # which.min(eval.min.row)
  rA.index = which(eval.min.col == min(eval.min.col )) # which.max(eval.min.col)
  
  gamma.hat = gamma.grid[g.index][1]
  rhoA.hat = rhoA.grid[rA.index]
  
  
  est.sim = estimate_parameters(gamma.hat = gamma.hat, rhoA.hat = rhoA.hat, 
                                var = var, shift = shift, numBlock = numBlock, 
                                wd.length = wd.length, TimeBlock=TimeBlock)
  
  ##Also estimate the mortality rates
  mort_lowess_rate = sapply(1:length(dat.sim),function(i){
    lowess((var$delD[,i]/var$H[-length(var$H[,i]),i]),f=1/16)$y})
  mort_lm_rate = sapply(1:length(dat.sim),function(i) {
    mort_lm = lm(var$delD[,i] ~ var$H[-length(var$H[,i]),i] -1)
    return(mort_lm$coefficients[1])})
  
  return(list(est.sim = est.sim,
              mortality_rate = list(mort_lowess_rate = mort_lowess_rate, mort_lm_rate = mort_lm_rate)))
  
}, mc.cores = ncore)

proc.time()-time

level <- 0.05


### Collect all the parameter estimated from each bootstrap run and write them into an .Rda file
## collect gamma
gamma.bs <- sapply(1:sim.BS.N, function(i){
  return(est.BS.list[[i]]$est.sim$gamma.hat)
})


## collect rhoA
rhoA.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$est.sim$rhoA.hat)})

## collect At
At.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$est.sim$A.hat.t)})
dim(At.bs)
At.ci.bs <- t(sapply(1:dim(At.bs)[1], function(i){
  sort(At.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

## collect Ct
Ct.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$est.sim$cum.new.infect.t)})
dim(Ct.bs)
Ct.ci.bs <- t(sapply(1:dim(Ct.bs)[1], function(i){
  sort(Ct.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

## point-wise CI for rhoH.t, phi.t
rhoH.smooth.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$est.sim$rhoH.hat.smooth)})
dim(rhoH.smooth.bs)
#range(rhoH.bs)
rhoH.ci.bs <- t(sapply(1:dim(rhoH.smooth.bs)[1], function(i){
  sort(rhoH.smooth.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

phi.smooth.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$est.sim$phi.hat.smooth)})
dim(phi.smooth.bs)
#range(phi.smooth.bs)
phi.ci.bs <- t(sapply(1:dim(phi.smooth.bs)[1], function(i){
  sort(phi.smooth.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

## point-wise CI for del.t.lowess
del.lowess.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$mortality_rate$mort_lowess_rate)})
dim(del.lowess.bs)
#range(rhoH.bs)
del.t.ci.bs <- t(sapply(1:dim(del.lowess.bs)[1], function(i){
  sort(del.lowess.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

##CI for del.t.lm
del.lm.bs <- sapply(1:sim.BS.N, function(i){return(est.BS.list[[i]]$mortality_rate$mort_lm_rate)})
##CI for alpha.t

## collect alpha.t
alpha.t.bs <- sapply(1:sim.BS.N, function(i){return(
  sapply(1:10,function(j){
    
    return((est.BS.list[[i]]$est.sim$sqrt.alpha.kappa.t[,j]/ten.kappa[1:dim(est.BS.list[[i]]$est.sim$sqrt.alpha.kappa.t)[1],j])^2)
  })
  )})
dim(alpha.t.bs)

alpha.ci.bs <- t(sapply(1:dim(alpha.t.bs)[1], function(i){
  sort(alpha.t.bs[i,])[round(c((sim.BS.N+1)*level/2, (sim.BS.N+1)*(1-level/2)))]}))

###


est.bs = list(del.lm.bs = del.lm.bs, del.lowess.bs  = del.lowess.bs , 
              phi.smooth.bs  = phi.smooth.bs , rhoH.smooth.bs = rhoH.smooth.bs,alpha.t.bs = alpha.t.bs,
              Ct.bs = Ct.bs ,At.bs= At.bs, gamma.bs = gamma.bs, rhoA.bs = rhoA.bs)
save(est.BS.list, est.bs,file = "sampling_agg_sim.RData")







