## 10/16/20

## modify the objective function to apply to 10 states
## at the same time
## now all states eg Ct and delC are all matrix



objective.rhoH = function(arg_para, var, grA.data=c(0.02,0.08), TimeSet){
  delC = var$delC
  delR = var$delR
  delH = var$delH
  delD = var$delD
  Q = var$Q
  H = var$H
  TtoH = var$TtoH
  
  rhoH  = arg_para 
  gamma = grA.data[1]
  rhoA = grA.data[2]
  
  respR = sqrt(pmax(delR[TimeSet,],0)) 
  predR = sqrt(pmax(rhoH * H[TimeSet,] + rhoA * Q[TimeSet,],0)) 
  
  mseR = mean(c((respR - predR)^2))# use c() to vectorize from matrix
  return(mseR)
}


objective.phi = function(arg_para, var, grA.data = c(0.02,0.08), TimeSet){
  
  delC = var$delC
  delR = var$delR
  delH = var$delH
  delD = var$delD
  Q = var$Q
  H = var$H
  TtoH = var$TtoH
  
  phi   = arg_para 
  gamma = grA.data[1]
  rhoA = grA.data[2]
  
  respH = sqrt(pmax(delH[TimeSet,] + delD[TimeSet,] + delR[TimeSet,],0)) 
  predH = sqrt(pmax((gamma + rhoA) * Q[TimeSet,] + gamma * delC[TimeSet,] 
                    /(gamma + phi * TtoH[TimeSet,]),0)) 
  
  mseH = mean(c((respH - predH)^2))# use c() to vectorize from matrix
  return(mseH)
}

objective.grA = function(var, eval_para, TimeBlock){
  
  gamma = eval_para[1]
  rhoA = eval_para[2]
  
  crit = 0
  for(k in 1:length(TimeBlock)){
    TimeSet = TimeBlock[[k]]
    crit = crit + optim(par=0.05,
                        fn=objective.rhoH,
                        gr=NULL,
                        var, c(gamma,rhoA), TimeSet,
                        lower=rhoA/3,upper=2*rhoA,
                        method = "L-BFGS-B",
                        control=list(ndeps=0.0001))$value
    
    crit = crit + optim(par=0.05,
                        fn=objective.phi,
                        gr=NULL, var, c(gamma,rhoA), TimeSet,
                        lower=0.001,upper=0.5/max(var$TtoH),
                        method = "L-BFGS-B",
                        control=list(ndeps=0.0001))$value
  }
   #print(paste("gamma = ",gamma,", rhoA =",rhoA, ", crit = ",crit))
  return(crit)
}
