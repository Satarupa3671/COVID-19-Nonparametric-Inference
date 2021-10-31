## estimation based on observed data
estimated_states = function(dat.int, true.para,Time.use){
  
  A <- dat.int$At
  S <- c(dat.int$St,rep(0,Time.use-1))
  R <- c(dat.int$Rt,rep(0,Time.use-1))
  Q <- c(dat.int$Qt,rep(0,Time.use-1))
  H <- c(dat.int$Ht,rep(0,Time.use-1))
  D <- c(dat.int$Dt,rep(0,Time.use-1))
  C <- c(dat.int$Ct,rep(0,Time.use-1))
  Test <- dat.int$Tt[1:Time.use]
  kappa <- dat.int$kappa[1:Time.use]
  
  gamma <- true.para$gamma.true
  rhoA <- true.para$rhoA.true
  rhoH.t <- true.para$rhoH.t.true
  alpha <- true.para$alpha.true
  delta.t <- true.para$delta.t.true[1:Time.use]
  phi.t <- true.para$phi.t.true
  
  theta.t <- NULL
  
  for (t in 1:(Time.use-1)){
    
    theta.t[t] <- phi.t[t]*abs(diff(Test)[t])/H[t]
    
    DelC.t = (theta.t[t] + gamma)*A[t]
    #C[t+1] = pmin((C[t]+ DelC.t),Max.pop)
    C[t+1] = abs(C[t]+ DelC.t)
    
    DelQ.t = theta.t[t]*A[t] - (gamma + rhoA)*Q[t]
    #Q[t+1] = pmin((Q[t] + DelQ.t), Max.pop)
    Q[t+1] = abs(Q[t] + DelQ.t)
    
    DelH.t = gamma*A[t] + gamma*Q[t] - rhoH.t[t]*H[t] - delta.t[t]*H[t]
    #H[t+1] = pmin((H[t] + DelH.t), C[t])
    H[t+1] = abs(H[t] + DelH.t)
    
    DelRH.t = rhoH.t[t]*H[t] 
    DelRQ.t = rhoA*Q[t]
    DelRA.t = rhoA*A[t]
    
    DelR.t = DelRH.t + DelRQ.t #+DelRA.t
    #R[t+1] = pmin((R[t] + DelR.t),C[t])
    R[t+1] = abs(R[t] + DelR.t)
    
    DelD.t = delta.t[t]*H[t]
    #D[t+1] = pmin((D[t] + DelD.t),C[t])
    D[t+1] = abs(D[t] + DelD.t)
    
    DelS.t = -kappa[t]^2 *alpha*A[t] * S[t]/(S[t]+R[t]+A[t])
    #S[t+1] = pmin((S[t] + DelS.t),Max.pop)
    S[t+1] = abs(S[t] + DelS.t)
  }
  
  # C = lowess(C, f = 1/16)$y
  # R = lowess(R, f = 1/16)$y
  # Q = lowess(Q, f = 1/16)$y
  # H = lowess(H, f = 1/16)$y
  # D = lowess(D, f = 1/16)$y
  dat.sim <- data.frame(At=A) %>%
    mutate(St = S,Rt = R,
           Qt = Q,Ht = H,
           Dt = D,Ct = C,
           Tt = Test[1:Time.use],
           kappa = kappa[1:Time.use])
  
  return(list(data=dat.sim,theta=theta.t))
  # return(list(A = A, S = S, R = R, Q = Q, H = H, D = D, C = C, Test = Test, kappa = kappa,
  #             theta.t = theta.t))
}


## estimate parameters given 10 bs process
estimation_module_agg <- function(dat,numBlock,
                                  wd.length,Time){
  
  
  ##Variables: all will become a Time * 10 matrix 
  # delC = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Ct)})
  # delR = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Rt)})
  # delH = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Ht)})
  # delD = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Dt)})
  # Q = sapply(1:length(dat),function(i){dat[[i]]$data$Qt})
  # H = sapply(1:length(dat),function(i){dat[[i]]$data$Ht})
  # delTests = sapply(1:length(dat),function(i){
  #   lowess((diff(dat[[i]]$data$Tt)), f = 1/16)$y})
  # TtoH_deno = sapply(1:length(dat),function(i){
  #   lowess((dat[[i]]$data$Ht[-dim(dat[[i]]$data)[1]]), f = 1/16)$y})
  # TtoH = delTests/TtoH_deno
  # 
  # var = list(delC = delC, delR = delR, delH = delH, 
  #            delD = delD, Q = Q, H = H, TtoH = TtoH)
  delC = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Ct)})
  delR = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Rt)})
  delH = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Ht)})
  delD = sapply(1:length(dat),function(i){diff(dat[[i]]$data$Dt)})
  Q = sapply(1:length(dat),function(i){dat[[i]]$data$Qt})
  H = sapply(1:length(dat),function(i){dat[[i]]$data$Ht})
  R = sapply(1:length(dat),function(i){dat[[i]]$data$Rt})
  D = sapply(1:length(dat),function(i){dat[[i]]$data$Dt})
  C = sapply(1:length(dat),function(i){dat[[i]]$data$Ct})
  S = sapply(1:length(dat),function(i){dat[[i]]$data$St})
  Test = sapply(1:length(dat),function(i){dat[[i]]$data$Tt})
  delTests = sapply(1:length(dat),function(i){
    lowess((diff(dat[[i]]$data$Tt)), f = 1/16)$y})
  TtoH_deno = sapply(1:length(dat),function(i){
    lowess((dat[[i]]$data$Ht[-dim(dat[[i]]$data)[1]]), f = 1/16)$y})
  TtoH = delTests/TtoH_deno
  var = list(delC = delC, delR = delR, delH = delH, 
             delD = delD, Q = Q, H = H, TtoH = TtoH,R=R,
             D=D,C=C, Test = Test,S=S)
  
  
  n = Time
  TimeBlock = NULL
  shift = floor(((n-1) - wd.length)/(numBlock - 1))
  for(j in 1:numBlock){
    TimeBlock[[j]] = (1 + shift * (j-1)) : min(shift * (j-1) + wd.length,(n-1))
  }
  TimeBlock[[numBlock]] <- TimeBlock[[numBlock]][-length(TimeBlock[[numBlock]])]
  
  gamma.grid = seq(0.001,0.05,0.001)# for set2:seq(0.001,0.08,0.001)
  rhoA.grid = seq(0.002,0.3,0.002)
  eval.mat = matrix(0,length(gamma.grid),length(rhoA.grid))
  
  for(i in 1:length(gamma.grid)){
    for(j in 1:length(rhoA.grid)){
      eval.mat[i,j] = objective.grA(var,c(gamma.grid[i],rhoA.grid[j]),TimeBlock)
    }
  }
  
  eval.gamma = apply(eval.mat,1,min)
  eval.rhoA = apply(eval.mat,2,min)
  
  g.index = which(eval.gamma == min(eval.gamma )) # which.min(eval.gamma)
  rA.index = which(eval.rhoA == min(eval.rhoA )) # which.min(eval.rhoA)
  
  gamma.hat = gamma.grid[g.index][1]
  rhoA.hat = rhoA.grid[rA.index]
  #print(c(gamma.hat, rhoA.hat))
  mort_lowess_rate = sapply(1:length(dat),function(i){
    lowess((var$delD[,i]/var$H[-length(var$H[,i]),i]),f=1/16)$y})
  mort_lm_rate = sapply(1:length(dat),function(i) {
    mort_lm = lm(var$delD[,i] ~ var$H[-length(var$H[,i]),i] -1)
    return(mort_lm$coefficients[1])})
  
  
  est1 = estimate_parameters(gamma.hat = gamma.hat, rhoA.hat = rhoA.hat[1], 
                             var = var, shift = shift, 
                             numBlock = numBlock, 
                             wd.length = wd.length, TimeBlock = TimeBlock)
  
  est.para = list(gamma.true = est1$gamma.hat,
                  rhoA.true = est1$rhoA.hat,
                  rhoH.t.true = est1$rhoH.hat.smooth,
                  alpha.true = (mean(est1$sqrt.alpha.kappa.t/dat[[1]]$data$kappa[1:(length(est1$sqrt.alpha.kappa.t)/10)]))^2,
                  delta.t.true = mort_lowess_rate,
                  phi.t.true = est1$phi.hat.smooth)
  
  dat.int.est.list <- lapply(1:10,function(i){
    re <- list(At=c(est1$A.hat.t[,i], est1$A.hat.t[dim(est1$A.hat.t)[1],i]),#sample(est1$A.hat.t[,i],1)),
               Rt = var$R[1,i],Qt = var$Q[1,i],
               Ht = var$H[1,i],Dt = var$D[1,i],
               Ct = var$C[1,i],Max.pop = var$Max.pop[i],
               St = var$S[1,i],
               Tt = var$Test[,i],
               kappa = dat.sim.new[[1]]$data$kappa)#10/17: dat.sim[[1]]$data$kappa) they are the same..
    return(re)
  })
  
  est_states = lapply(1:10, function(i){
    
    return(estimated_states(dat.int = dat.int.est.list[[i]],
                            true.para = est.para,Time.use=length(dat.int.est.list[[i]]$At)))
  })
  return(list(est_para = list(est = est1, 
                              mortality = list(mort_lm_rate = mort_lm_rate, mort_lowess_rate = mort_lowess_rate)),
              est_states = est_states, 
              tuning = list(numBlock = numBlock, wd.length = wd.length, TimeBlock = TimeBlock)
              
  ))
}


estimate_parameters = function(gamma.hat, rhoA.hat, var, shift, numBlock, wd.length, TimeBlock){
  rhoH.est = list()
  phi.est = list()
  TtoH = var$TtoH
  delC = var$delC
  
  for(l in 1:numBlock){
    TimeSet = TimeBlock[[l]]
    rhoH.est[[l]] = optim(par = 0.05,
                          fn = objective.rhoH,
                          gr = NULL,
                          var = var, grA.data = c(gamma.hat, rhoA.hat), TimeSet,
                          lower=.0001, #rhoA.hat/3,
                          upper=5*rhoA.hat,
                          control=list(ndeps=0.0001),
                          method = "L-BFGS-B")
    phi.est[[l]] = optim(par = .1,
                         fn = objective.phi,
                         gr = NULL,
                         var =var, grA.data = c(gamma.hat, rhoA.hat), TimeSet,
                         lower=0,upper=0.5/max(TtoH),
                         control=list(ndeps=0.0001),
                         method = "L-BFGS-B")
  }
  
  rhoH.hat.vec = numeric(0)
  phi.hat.vec = numeric(0)
  
  for(l in 1:numBlock){
    rhoH.hat.vec[l] = rhoH.est[[l]]$par
    phi.hat.vec[l] = phi.est[[l]]$par
  }
  round(rhoH.hat.vec,4)
  round(phi.hat.vec,4)
  
  
  ## Evaluation of phi(t) at every time point
  
  n.t = shift * (numBlock-1) + wd.length   # number of time points used in the fitting procedure
  phi.hat.mult = rep(0,n.t-1)
  
  umat = matrix(0,numBlock,n.t-1)
  for(l in 1:numBlock){
    umat[l,TimeBlock[[l]]] = rep(phi.hat.vec[l],length(TimeBlock[[l]]))#wd.length)
  }
  
  phi.hat.mult = apply(umat, 2, sum) / apply(umat !=0, 2, sum)
  phi.hat.smooth = lowess(phi.hat.mult,f=1/16)$y
  
  ## Evaluation of rhoH(t) at every time point
  rhoH.hat.mult = rep(0,n.t-1)
  
  umat2 = matrix(0,numBlock,n.t-1)
  for(j in 1:numBlock){
    umat2[j,TimeBlock[[j]]] = rep(rhoH.hat.vec[j],length(TimeBlock[[j]]))
  }
  
  rhoH.hat.mult = apply(umat2, 2, sum) / apply(umat2 !=0, 2, sum)
  rhoH.hat.smooth = lowess(rhoH.hat.mult,f=1/16)$y
  
  ################################
  ## Computation of various rate parameters and fitted value of A(t)
  
  n.t = length(phi.hat.smooth)
  range.t = 1:n.t
  gamma.t = gamma.hat[1]
  rhoA.t = rhoA.hat
  
  theta.hat.t = phi.hat.smooth * TtoH[range.t,]  # estimate of theta(t) [conifirmation fraction]
  
  delC.low = sapply(1:dim(delC)[2],function(i){
    lowess(delC[range.t,i],f=1/16)$y})  # lowess smoothing of Delta(C(t))
  CIR.t = diff(delC.low)/delC.low[-dim(delC.low)[1],]    # Crude Infection Rate =  \Delta^2 C(t)/\Delta C(t)
  
  RCCF.t = sapply(1:dim(delC)[2],function(i){
    lowess(diff(theta.hat.t[,i])/(theta.hat.t[-n.t,i]+gamma.t),f=1/6)$y})  # relative change in confirmation fraction eta(t) = theta(t)+gamma
  
  NIR.hat.t = sapply(1:dim(delC)[2],function(i){
    lowess((CIR.t[,i] - RCCF.t[,i])/(1 + RCCF.t[,i]) ,f=1/8)$y})  #estimated Net Infection Rate
  
  nu.hat.t = NIR.hat.t + theta.hat.t[-n.t,] + gamma.t + rhoA.t  # Infection rate
  
  sqrt.alpha.kappa.t = sapply(1:dim(delC)[2],function(i){
    sqrt(pmax(0,nu.hat.t[,i]))}) # estimate of sqrt{alpha} * kappa(t)
  
  A.hat.t = delC[range.t,]/(theta.hat.t[,] + gamma.t)  # fitted value of A(t)
  
  delA.hat.t = NIR.hat.t * A.hat.t[-n.t,]   # model-based growth in A(t)
  # NOTE: delA.hat.t may be different from Delta(A.hat.t),
  # which will indicate a model mismatch
  
  new.infect.t = nu.hat.t * A.hat.t[-n.t,]  # newly infected
  cum.new.infect.t = sapply(1:dim(delC)[2],function(i){
    cumsum(new.infect.t[,i])})  # cumulative number of newly infected
  est = list( gamma.hat = gamma.hat, rhoA.hat =rhoA.hat,
              rhoH.hat.vec = rhoH.hat.vec, rhoH.hat.mult = rhoH.hat.mult,
              rhoH.hat.smooth = rhoH.hat.smooth,
              phi.hat.vec = phi.hat.vec,  phi.hat.mult = phi.hat.mult,
              phi.hat.smooth = phi.hat.smooth, 
              theta.hat.t = theta.hat.t, 
              CIR.t = CIR.t, RCCF.t = RCCF.t,
              NIR.hat.t = NIR.hat.t, nu.hat.t = nu.hat.t,
              sqrt.alpha.kappa.t = sqrt.alpha.kappa.t, 
              A.hat.t = A.hat.t,
              delA.hat.t = delA.hat.t,
              new.infect.t = new.infect.t ,
              cum.new.infect.t = cum.new.infect.t,
              n.t = n.t, TimeSet = TimeSet)
  return(est)
}

## generate the residual bootstrap

resi_adj_boot = function(res, est_del_st, samp_ind){
  bw_lowess = 1/2
  resi.adj = res - mean(res)
  n = length(resi.adj)
  #s = sqrt(pmax(0,lowess((resi.adj)^2,f = bw_lowess)$y))
  ref = c(resi.adj[n:1], resi.adj, resi.adj[n:1])
  s2 = sqrt(pmax(0,lowess(ref^2, f = bw_lowess/3)$y)[(n+1):(2*n)])
  scale_factor = median(s2)/2
  resi.adj.scale = resi.adj/(s2+scale_factor)
  zeta = 0
  n = length(resi.adj.scale)
  # for(t in 2:n){
  #   zeta =+ sum((resi.adj.scale[t]*resi.adj.scale[t-1]))/sum(resi.adj.scale[t]^2)
  # }
  zeta = sum(resi.adj.scale[-1]*resi.adj.scale[-n])/sum(resi.adj.scale[-1]^2)
  zeta = ifelse( zeta <= .98, zeta, .98)
  zeta = ifelse( zeta >= -.98, zeta, -.98)
  hat.resi.adj.scale = c(resi.adj.scale[1],rep(0,(n-1)))
  for(t in 2:n){
    hat.resi.adj.scale[t] = zeta*resi.adj.scale[t-1]
  }
  resi_innov = resi.adj.scale - hat.resi.adj.scale 
  resi_innov.samp = resi_innov[samp_ind]
  resi.adj.scale.boot = resi_innov.samp + hat.resi.adj.scale
  resi.adj_boot = (s2+scale_factor)*resi.adj.scale.boot
  boot_samp = mean(res) + resi.adj_boot + est_del_st
  return(boot_samp)
}


resi_bootstrap = function(observed, estimated,theta, Time.use, Max.pop){
  
  resi.Ct = diff(sqrt(observed$Ct))[1:(Time.use-1)] - diff(sqrt(estimated$Ct))[1:(Time.use-1)] 
  resi.Rt = diff(sqrt(observed$Rt))[1:(Time.use-1)] - diff(sqrt(estimated$Rt))[1:(Time.use-1)] 
  resi.Dt = diff(sqrt(observed$Dt))[1:(Time.use-1)] - diff(sqrt(estimated$Dt))[1:(Time.use-1)] 
  resi.Ht = diff(sqrt(observed$Ht))[1:(Time.use-1)] - diff(sqrt(estimated$Ht))[1:(Time.use-1)]  ##A
  resi.Qt = diff(sqrt(observed$Qt))[1:(Time.use-1)] - diff(sqrt(estimated$Qt))[1:(Time.use-1)]  ##A
  
  samp_ind = sample(x = 1:(length(resi.Ct)), size = (length(resi.Ct)), replace = TRUE)
  r = sample(1:length(resi.Ht),1)
  
  
  delCt.adj.boot = resi_adj_boot(res = resi.Ct, est_del_st = diff(sqrt(estimated$Ct))[1:(Time.use-1)] , samp_ind)
  delDt.adj.boot = resi_adj_boot(resi.Dt, diff(sqrt(estimated$Dt))[1:(Time.use-1)] , samp_ind)
  delRt.adj.boot = resi_adj_boot(resi.Rt, diff(sqrt(estimated$Rt))[1:(Time.use-1)] , samp_ind)
  delHt.adj.boot = resi_adj_boot(resi.Ht, diff(sqrt(estimated$Ht))[1:(Time.use-1)] , samp_ind)
  delQt.adj.boot = resi_adj_boot(resi.Qt, diff(sqrt(estimated$Qt))[1:(Time.use-1)] , samp_ind)
  
  
  A = estimated$A
  C = c(observed$C[1],(cumsum(delCt.adj.boot) + sqrt(observed$C[1]))^2)
  R = c(observed$R[1],(cumsum(delRt.adj.boot) + sqrt(observed$R[1]))^2)
  D = c(observed$D[1],(cumsum(delDt.adj.boot) + sqrt(observed$D[1]))^2)
  H = c(observed$H[1],(cumsum(delHt.adj.boot) + sqrt(observed$H[1]))^2)
  Q = c(observed$Q[1],(cumsum(delQt.adj.boot) + sqrt(observed$Q[1]))^2)
  RA = estimated$R - estimated$Rt
  S = pmax((rep(Max.pop,(Time.use)) - A[1:(Time.use)] - R[1:(Time.use)] - RA[1:(Time.use)] -Q[1:(Time.use)] -H[1:(Time.use)] -D[1:(Time.use)] -C[1:(Time.use)]),1)
  
  
  
  delC = diff(C)
  delR = diff(R)
  delH = diff(H)
  delD = diff(D)
  delQ = diff(Q)
  delS = diff(S)
  
  bw_lowess = 1/16
  delTests = lowess((diff(observed$Tt)), f = bw_lowess)$y
  TtoH_deno = lowess(H, f = bw_lowess)$y
  n_eff = min(length(delTests), length(TtoH_deno))
  TtoH = delTests[1:n_eff]/TtoH_deno[1:n_eff]
  
  
  data.bs <- data.frame(At=A) %>%
    mutate(St = S , Rt = R,
           Qt = Q ,Ht = H,
           Dt = D ,Ct = C,
           Tt = observed$Tt[1:Time.use],
           kappa = observed$kappa[1:Time.use])
  
  return(list(data = data.bs, theta = theta))
}
