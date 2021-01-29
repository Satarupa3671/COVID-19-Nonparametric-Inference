##Code for optimizing the profile loss function to estimate gamma and rhoA.

## Objective function for estimating rhoH given (gamma, rhoA) within each time block (kernel support)
bw_lowess = 1/32
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
  
  respR = sqrt(pmax(delR[TimeSet],0))
  respR[is.na(respR)] = 0
  predR = sqrt(pmax(rhoH * H[TimeSet] + rhoA * Q[TimeSet],0))
  predR[is.na(predR)] = 0
  mseR = mean((respR - predR)^2)
  return(mseR)
}


## Objective function for estimating phi given (gamma, rhoA)  within each time block (kernel support)

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
  
  respH = sqrt(pmax(delH[TimeSet] + delD[TimeSet] + delR[TimeSet],0))
  predH = sqrt(pmax((gamma + rhoA) * Q[TimeSet] + gamma * delC[TimeSet]
                    /(gamma + phi * TtoH[TimeSet]),0))
  respH[is.na(respH)] = 0
  predH[is.na(predH)] = 0
  mseH = mean((respH - predH)^2)
  return(mseH)
}

## Profile loss function for (gamma,rhoA) summing contributions across numBlock time blocks (kernel supports)

objective.grA = function(var, eval_para,TimeBlock){
  
  gamma = eval_para[1]
  rhoA = eval_para[2]
  
  crit = 0
  for(k in 1:length(TimeBlock)){
    TimeSet = TimeBlock[[k]]
    crit = crit + optim(par=0.05,
                        fn=objective.rhoH,
                        gr=NULL,
                        var, c(gamma,rhoA), TimeSet,
                        lower=rhoA/3,upper=rhoA,
                        method = "L-BFGS-B",
                        control=list(ndeps=0.0001))$value
    
    crit = crit + optim(par=0.05,
                        fn=objective.phi,
                        gr=NULL, var, c(gamma,rhoA), TimeSet,
                        lower=0.001,upper=0.3/max(var$TtoH),
                        method = "L-BFGS-B",
                        control=list(ndeps=0.0001))$value
  }
  #print(paste("gamma = ",gamma,", rhoA =",rhoA, ", crit = ",crit))
  return(crit)
}
print(bw_lowess)