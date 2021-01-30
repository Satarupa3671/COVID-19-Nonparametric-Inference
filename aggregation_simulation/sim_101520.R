rm(list=ls())
library(tidyverse)
#### simulation study (eq (1)-(6))

## initial parameter:

Time <- 200

true.para <- list(gamma.true = 0.03,
                  rhoA.true = 0.06,
                  rhoH.t.true = c(rep(0.03,20),rep(0.04,40),rep(0.05,50),
                    rep(0.06,45),rep(0.05,45)), 
                  alpha.true = 0.2,
                  delta.t.true = c(rep(0.007,50),rep(0.005,70),rep(0.003,60),
                                   rep(0.004,20)),
                  phi.t.true = 0.001+1:Time*0.00005)

set.seed(123)
ten.initial <- data.frame(At = sample(100:300,10,replace = T),
                          St = sample((3*10^6):(5*10^6),10,replace = T),
                          Rt = sample(1:10,10,replace = T),
                          Qt = sample(50:100,10,replace = T),
                          Ht = sample(0:15,10,replace = T),
                          Dt = sample(0:7,10,replace = T),
                          Ct = sample(110:180,10,replace = T)) %>%
  mutate(Max.pop = At+St+Rt+Qt+Ht+Dt+Ct)

ten.Test <- data.frame(k1 = cumsum(c(sample(1:300,30),
                                     sample(1000:3000,50),
                                     sample(4000:8000,100),
                                     sample(7500:12000,20))),
                       k2 = cumsum(c(sample(1:200,40),
                                     sample(400:1000,50),
                                     sample(2000:6000,100),
                                     sample(3000:6000,10))),
                       k3 = cumsum(c(sample(100:300,20),
                                     sample(1200:6000,100),
                                     sample(3000:5000,60),
                                     sample(5000:7000,20))),
                       k4 = cumsum(c(sample(50:300,30),
                                     sample(2000:4000,80),
                                     sample(2000:2500,70),
                                     sample(5000:8000,20))),
                       k5 = cumsum(c(sample(40:100,30),
                                     sample(300:1300,50),
                                     sample(2000:5000,90),
                                     sample(3000:4000,30))),
                       k6 = cumsum(c(sample(1:100,20),
                                     sample(500:3000,100),
                                     sample(2000:5000,40),
                                     sample(2500:4000,40))),
                       k7 = cumsum(c(sample(30:200,40),
                                     sample(5000:1500,50),
                                     sample(900:4000,60),
                                     sample(3000:5000,50))),
                       k8 = cumsum(c(sample(10:200,30),
                                    sample(500:1500,60),
                                    sample(2000:4000,90),
                                    sample(3000:5000,20))),
                       k9 = cumsum(c(sample(50:300,40),
                                     sample(1200:5000,120),
                                     sample(2500:4000,40))),
                       k10 = cumsum(c(sample(10:400,30),
                                      sample(1000:3500,70),
                                      sample(2000:5000,40),
                                      sample(3000:4000,60))))
View(ten.Test)
ten.kappa <- data.frame(k1 = lowess(c(rep(1.2,20),rep(0.8,40),rep(0.4,30),
                                      rep(0.6,20),rep(0.9,10),rep(1.3,30),
                                      rep(0.8,20),rep(1,20),rep(1.2,10)),
                                    f=1/16)$y,
                        k2 = lowess(c(rep(1.1,20),rep(1.3,40),rep(0.8,30),
                                      rep(0.9,20),rep(1.1,10),rep(1.5,10),
                                      rep(1.2,20),rep(0.9,30),rep(1,20)),
                                    f=1/16)$y,
                        k3 = lowess(c(rep(1,30),rep(0.7,30),rep(0.4,30),
                                      rep(0.3,20),rep(0.25,10),rep(0.7,40),
                                      rep(1.2,20),rep(0.8,20)),
                                    f=1/16)$y,
                        k4 = lowess(c(rep(0.9,20),rep(1.2,40),rep(1.4,30),
                                      rep(0.8,20),rep(1,10),rep(1.1,30),
                                      rep(0.8,20),rep(0.9,20),rep(1.2,10)),
                                    f=1/16)$y,
                        k5 = lowess(c(rep(1.3,40),rep(0.6,50),rep(0.7,30),
                                      rep(0.3,20),rep(0.5,30),rep(0.8,10),
                                      rep(1,20)),
                                    f=1/16)$y,
                        k6 = lowess(c(rep(1.2,20),rep(0.6,40),rep(0.8,30),
                                      rep(1.4,20),rep(1,10),rep(0.7,30),
                                      rep(1,20),rep(0.9,30)),
                                    f=1/16)$y,
                        k7 = lowess(c(rep(1,20),rep(0.9,30),rep(0.5,30),
                                      rep(0.2,20),rep(0.4,20),rep(0.6,30),
                                      rep(0.5,20),rep(0.7,20),rep(0.9,10)),
                                    f=1/16)$y,
                        k8 = lowess(c(rep(1.5,20),rep(1.1,40),rep(0.9,30),
                                      rep(0.6,20),rep(0.9,10),rep(1,30),
                                      rep(0.8,20),rep(1,20),rep(1.2,10)),
                                    f=1/16)$y,
                        k9 = lowess(c(rep(0.8,20),rep(0.6,40),rep(0.2,30),
                                      rep(0.4,20),rep(0.5,10),rep(0.7,30),
                                      rep(0.4,20),rep(0.8,20),rep(0.9,10)),
                                    f=1/16)$y,
                        k10 = lowess(c(rep(1,30),rep(0.7,20),rep(0.3,50),
                                      rep(0.4,20),rep(0.5,10),rep(0.7,10),
                                      rep(0.8,30),rep(1,20),rep(1.2,10)),
                                    f=1/16)$y)




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

dat.sim <- lapply(1:10, function(i){
  
  return(covid_simulation(dat.int.list[[i]],true.para,Time=Time))
})
View(dat.sim[[1]]$data) ## !!! issue: At may decrease to 0 very fast and then results in NA
## maybe parameters are not reasonable?


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
numBlock = 10#8
wd.length =  50#28# 30
n = Time
shift = floor(((n-1) - wd.length)/(numBlock - 1))
#
for(j in 1:numBlock){
  TimeBlock[[j]] = (1 + shift * (j-1)) : min(shift * (j-1) + wd.length,(n-1))
}



source('fitting_function_101520.R')

#### Grid for evaluation of profile loss for (gamma, rhoA)  ###

gamma.grid = seq(0.001,0.003,0.001)
rhoA.grid = seq(0.002,0.004,0.002)
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

c(gamma.hat,rhoA.hat) #

round(eval.min.row/numBlock,2)
round(eval.min.col/numBlock,2)
gamma.hat/(gamma.hat+rhoA.hat)

plot(gamma.grid, eval.min.row,ylab='cost',
     xlab = expression(paste(gamma, ' grid')),
     cex.lab=1.5, cex.axis=1.3)
plot(rhoA.grid, eval.min.col,ylab='cost',
     xlab = expression(paste(rho[A], ' grid')),
     cex.lab=1.5, cex.axis=1.3)


source('est_plot_func_101520.R')

est.sim = estimate_parameters(gamma.hat = gamma.hat, rhoA.hat = rhoA.hat, 
                              var = var, shift = shift, numBlock = numBlock, 
                              wd.length = wd.length, TimeBlock=TimeBlock)
ggplot_est(est = est.sim,k=1) 
diagnostics(numBlock = numBlock,TimeBlock, 
            est = est.sim,var,k=1)


