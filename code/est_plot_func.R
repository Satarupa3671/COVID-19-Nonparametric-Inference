###Code for estimation of the parameters from the observed data using fitting_functions_0827.R. Also contains the ggplot codes for plotting the estimated values (displays the goodness of the fits).
##Estimate the other parameters for a given value of gamma and rhoA.

library(lemon)
bw_lowess = bw_lowess
estimate_parameters = function(gamma.hat, rhoA.hat, var, shift, numBlock, wd.length, TimeBlock){
  rhoH.est = list()
  phi.est = list()
  TtoH = abs(var$TtoH)
  delC = var$delC
  Q = var$Q
  delH = var$delH
  H = var$H
  delD = var$delD
  
  for(l in 1:numBlock){
    TimeSet = TimeBlock[[l]]
    rhoH.est[[l]] = optim(par = 0.05,
                          fn = objective.rhoH,
                          gr = NULL,
                          var = var, grA.data = c(gamma.hat, rhoA.hat), TimeSet,
                          lower=.0001, #rhoA.hat/3,
                          upper=2*rhoA.hat,
                          control=list(ndeps=0.0001),
                          method = "L-BFGS-B")
    phi.est[[l]] = optim(par = .05,
                         fn = objective.phi,
                         gr = NULL,
                         var =var, grA.data = c(gamma.hat, rhoA.hat), TimeSet,
                         lower=0,upper=0.3/max(TtoH),
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
  phi.hat.mult = rep(0,n.t)
  
  
  umat = matrix(0,numBlock,n.t)
  for(l in 1:numBlock){
    umat[l,TimeBlock[[l]]] = rep(phi.hat.vec[l],wd.length)
  }
  
  phi.hat.mult = apply(umat, 2, sum) / pmax(1,(apply(umat !=0, 2, sum)))
  phi.hat.smooth = lowess(phi.hat.mult,f=bw_lowess)$y
  
  ## Evaluation of rhoH(t) at every time point
  rhoH.hat.mult = rep(0,n.t)
  
  umat2 = matrix(0,numBlock,n.t)
  for(j in 1:numBlock){
    umat2[j,TimeBlock[[j]]] = rep(rhoH.hat.vec[j],wd.length)
  }
  
  rhoH.hat.mult = apply(umat2, 2, sum) / pmax(1, (apply(umat2 !=0, 2, sum)))
  rhoH.hat.smooth = lowess(rhoH.hat.mult,f=bw_lowess)$y
  
  
  ################################
  ## Computation of various rate parameters and fitted value of A(t)
  
  n.t = length(phi.hat.smooth)
  range.t = 1:(n.t-1)
  gamma.t = gamma.hat[1]
  rhoA.t = rhoA.hat
  
  theta.hat.t = pmax(phi.hat.smooth[1:(length(TtoH))] * TtoH,0.01) # estimate of theta(t) [confirmation fraction]
  
  delC.low = lowess(delC,f=bw_lowess)$y  # lowess smoothing of Delta(C(t))
  CIR.t = diff(delC.low)/delC.low[-length(delC.low)]    # Crude Infection Rate =  \Delta^2 C(t)/\Delta C(t)
  
  RCCF.t = lowess(diff(theta.hat.t)/(theta.hat.t[1:length(diff(theta.hat.t))]+gamma.t),f=1/8)$y  # relative change in confirmation fraction eta(t) = theta(t)+gamma
  
  NIR.hat.t = lowess((CIR.t - RCCF.t[1:length(CIR.t)])/(1 + RCCF.t[1:length(CIR.t)]) ,f=1/8)$y  #estimated Net Infection Rate
  
  nu.hat.t = abs(NIR.hat.t + theta.hat.t[1:length( NIR.hat.t)] + gamma.t + rhoA.t) # Infection rate
  
  sqrt.alpha.kappa.t = sqrt(pmax(0,nu.hat.t)) # estimate of sqrt{alpha} * kappa(t)
  
  A.hat.t = delC/pmax((theta.hat.t[1:length(delC)] + gamma.t),.01)  # fitted value of A(t)
  
  delA.hat.t = NIR.hat.t * A.hat.t[1:length(NIR.hat.t)]   # model-based growth in A(t)
  # NOTE: delA.hat.t may be different from Delta(A.hat.t),
  # which will indicate a model mismatch
  
  new.infect.t = nu.hat.t * A.hat.t[1:length(nu.hat.t)]  # newly infected
  cum.new.infect.t = cumsum(new.infect.t)  # cumulative number of newly infected
  
  mort_lowess_rate = lowess((delD/H[-length(H)]), f=bw_lowess)$y
  rhoH.final = pmax(0,
                    ((gamma.hat * (A.hat.t + Q[1:length(A.hat.t)]) - 
                        delH - mort_lowess_rate *H[1:length(mort_lowess_rate)])/H[1:length(mort_lowess_rate)]))
  
  est = list(gamma.hat = gamma.hat, rhoA.hat =rhoA.hat,
             rhoH.hat.vec = rhoH.hat.vec, rhoH.hat.mult = rhoH.hat.mult,
             rhoH.hat.smooth = rhoH.hat.smooth,
             rhoH.final = rhoH.final,
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



## Creates the ggplots for displaying the fits
plot_est = function(var,  raw_data , est_st, dates){
  
  est1 = est_st$est_para$est ##list containing various estimated epi parameters
  mort = est_st$est_para$mortality ##list containing estimated mortality paramteres
  
  
  ##Save the estimates of the parameters
  sink(sprintf("../Results for USA States/%s/%s_opt_gamma_rhoA.txt", st, st))
  print(paste("gamma.hat = ", est1$gamma.hat, "rhoA.hat = ", est1$rhoA.hat))
  print(paste("gamma.hat/(gamma.hat+rhoA.hat)=", est1$gamma.hat/(est1$gamma.hat + est1$rhoA.hat)))
  print(paste("rhoH.hat = ", round(est1$rhoH.hat.vec,4)))
  print(paste("phi.hat.smooth = ", round(est1$phi.hat.smooth),4))
  print(paste("rhoH.final = ", round(est1$rhoH.final),4))
  print(paste("mort_lm_rate = ", round(mort$mort_lm_rate,4)))
  print(paste("mort_lowess_rate = ", round(mort$mort_lowess_rate,4)))
  sink()
  
  
  TimeSet = est1$TimeSet
  n.t = est1$n.t
  Time.t = n.t + 10
  offH = 6#lags.para$offH   ## time offset for H (not important)
  
  plot.data <- as.data.frame(
    cbind(Time = (4+offH):(Time.t+offH-8),
          phi.hat.mult = est1$phi.hat.mult[1:length(est1$CIR.t)],
          phi.hat.smooth = est1$phi.hat.smooth[1:length(est1$CIR.t)],
          rhoH.final = pmax(0,
                            ((est1$gamma.hat * (est1$A.hat.t + var$Q[1:length(est1$A.hat.t)]) - 
                                var$delH - mort$mort_lowess_rate * var$H[1:length(mort$mort_lowess_rate)])/
                               var$H[1:length(mort$mort_lowess_rate)]))[1:length(est1$CIR.t)],
          
          theta.hat.t = est1$theta.hat.t[1:length(est1$CIR.t)],
          RCCF.t = est1$RCCF.t[1:length(est1$CIR.t)],
          CIR.t = est1$CIR.t,
          CIR.t.smooth = lowess(est1$CIR.t,f = 1/8)$y,
          NIR.hat.t = est1$NIR.hat.t[1:length(est1$CIR.t)],
          A.hat.t = pmax(est1$A.hat.t[1:length(est1$CIR.t)],0),
          diff.A.hat.t = diff(est1$A.hat.t),
          new.infect.t = pmax(est1$new.infect.t[1:length(est1$CIR.t)],0),
          delC = pmax(var$delC[1:length(est1$CIR.t)],0),
          cum.new.infect.t = pmax(est1$cum.new.infect.t,0),
          sqrt.alpha.kappa.t = est1$sqrt.alpha.kappa.t[1:length(est1$CIR.t)]))
  plot.data$ddates = as.Date(dates[1:length(est1$CIR.t)], format = "%m/%d/%Y") 
  
  
  
  common.length = length(est_st$est_states$C)-1
  df_CRHD = data.frame(Days = 1:common.length,
                       
                       raw_Obs_C = raw_data$C[1:common.length],
                       Obs_C = var$C[1:common.length],
                       Est_C = est_st$est_states$C[1:common.length],
                       
                       raw_Obs_DelC = raw_data$delC[1:common.length],
                       Obs_DelC = var$delC[1:common.length],
                       Est_DelC = lowess(diff(est_st$est_states$C),f = bw_lowess)$y,
                       
                       raw_Obs_D = raw_data$D[1:common.length],
                       Obs_D = var$D[1:common.length],
                       Est_D = est_st$est_states$D[1:common.length],
                       
                       raw_Obs_DelD = raw_data$delD[1:common.length],
                       Obs_DelD = var$delD[1:common.length],
                       Est_DelD = diff(est_st$est_states$D),
                       
                       raw_Obs_R = raw_data$R[1:common.length],
                       Obs_R = var$R[1:common.length],
                       Est_R = est_st$est_states$R.Reported[1:common.length],
                       
                       raw_Obs_DelR = raw_data$delR[1:common.length],
                       Obs_DelR = var$delR[1:common.length],
                       Est_DelR = diff(est_st$est_states$R.Reported),
                       
                       raw_Obs_H = raw_data$H[1:common.length],
                       Obs_H = var$H[1:common.length],
                       Est_H = est_st$est_states$H[1:common.length],
                       
                       raw_Obs_Q = raw_data$Q[1:common.length],
                       Obs_Q = var$Q[1:common.length],
                       Est_Q = est_st$est_states$Q[1:common.length],
                       
                       TtoH = var$TtoH[1:common.length],
                       raw_TtoH = raw_data$TtoH[1:common.length],
                       Test = var$Test[1:common.length],
                       raw_Test = raw_data$Test[1:common.length])
  
  df_CRHD$dates1 = as.Date(dates, format = "%m/%d/%Y")[1:common.length]
  
  
  
  ###Validation H
  
  break.vec <- c(df_CRHD$dates1[1],
                 seq(from = df_CRHD$dates1[1] , to = df_CRHD$dates1[length(df_CRHD$dates1)], 
                     by="month"), df_CRHD$dates1[length(df_CRHD$dates1)])
  
  validation_H = ggplot(data = df_CRHD)+
    geom_pointline(mapping = aes(y = raw_Obs_H, x = dates1, color = "Observed"), size =.8 ) +
    #geom_point(mapping = aes(y = raw_Obs_H, x = dates1, color = "Observed"), size =.4 ) +
    #geom_line(mapping = aes(y = raw_Obs_H, x = dates1, color = "Observed"), size =.4 ) +
    geom_line(mapping = aes(y = Est_H ,x = dates1, color = "Fitted"), size = 1) +
    scale_color_manual(values = c('Fitted' = 'red',
                                  'Observed' = 'black')) +
    labs(color = 'Y series', x = "Time", y = 'Reported Current Hosptitalization') +
    theme_bw() +
    scale_x_date(breaks =  break.vec,  #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text = element_text(size = 20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title = element_text(size = 26)) +
    theme(legend.position = "none", #c(.5,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text = element_text(size = 24))
  
  #validation_H
  
  ggsave(file = sprintf("../Results for USA States/%s/%s_validaton_H.pdf", st, st))
  
  
  ###RCCF
  
  
  # break.vec <- c( plot.data$ddates[1],
  #                seq(from = plot.data$ddates[1] , to = plot.data$ddates[length(plot.data$ddates)], 
  #                    by="month"), plot.data$ddates[length(plot.data$ddates)])
  
  RCCF = ggplot(plot.data,aes(x=ddates,y=RCCF.t))+
    geom_pointline(size = 1) +
    # geom_point(size=0.7)+
    # geom_line(y=lowess(est1$RCCF.t,f = 1/8)$y)+
     labs(y = 'RCCF(t)' ,x = "Time")+
    theme_bw()+
    scale_x_date(breaks = break.vec,# date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position= "top",
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  
  #RCCF
  ggsave(file =sprintf("../Results for USA States/%s/%s_RCCF.pdf", st, st))
  
  
  ###CIR NIR
  df_CIR_NIR = as.data.frame(cbind(Time = plot.data$Time, 
                                   CIR.t = plot.data$CIR.t, 
                                   NIR.t = plot.data$NIR.hat.t,
                                   CIR.t.sm = lowess(est1$CIR.t,f = bw_lowess)$y))
  
  df_CIR_NIR$dates1 = as.Date(dates[1:length(plot.data$CIR.t)], format = "%m/%d/%Y") 
  
  CIR_NIR = ggplot(df_CIR_NIR, aes(x=dates1, y=CIR.t))+
    #geom_point(size = .7)+ 
    geom_line(aes(y=CIR.t.sm,x= dates1,color="CIR(t)"),size=1 ) +
    geom_line(aes(y=NIR.t,x= dates1,color="NIR(t)"),size=1) +
    scale_color_manual(values = c(
      'CIR(t)' = 'black',
      'NIR(t)' = 'red')) +
    labs(color = 'Y series', x = "Time", y = 'CIR(t) and NIR(t)')+
    theme_bw()+
    scale_x_date(breaks = break.vec,# date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position= "none",# c(.4,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  #CIR_NIR
  ggsave(file =sprintf("../Results for USA States/%s/%s_NIR_CIR.pdf", st, st))
  
  
  
  ##CFR
  df_CFR = data.frame( obs_CFR = cumsum(100*var$delD[2:(n.t-1)])/ 
                                          cumsum(var$delC[2:(n.t -1)]),
                       est_CFR = cumsum(100*var$delD[2:(n.t-1)])/ 
                                          cumsum(est1$sqrt.alpha.kappa.t[2:(n.t-1)]*
                                                   est1$A.hat.t[2:(n.t-1)]))
  
  df_CFR$dates1 = as.Date(dates[2:(n.t-1)], format = "%m/%d/%Y") 
  
  
  CFR_plot = ggplot(df_CFR)+
    #geom_point(size = .7)+ 
    geom_line(aes(y = obs_CFR , x = dates1, color = "Observed"),size=1) +
    geom_line(aes(y = est_CFR , x = dates1, color = "Estimated"),size=1) +
    scale_color_manual(values = c(
      'Observed' = 'black',
      'Estimated' = 'red')) +
    labs(color = 'Y series', x = "Time", y = 'Case Fatality Rate')+
    theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position= "none",# c(.4,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  #CFR_plot
  ggsave(file =sprintf("../Results for USA States/%s/%s_CFR.pdf", st, st))
  
  
  ##Plot for R0 analogue
  R0_numer = (est_st$est_para$est$sqrt.alpha.kappa.t)^2 * unique(var$pop19)
  R0_denom = (est_st$est_states$S + est_st$est_states$A + est_st$est_states$R) * 
    (est_st$est_para$est$theta.hat.t + est_st$est_para$est$gamma.hat + est_st$est_para$est$rhoA.hat)
  R0_df = data.frame (Days = 1:length(R0_numer), R0 = R0_numer/R0_denom[1:length(R0_numer)])
  R0_df$dates1 = as.Date(dates, format = "%m/%d/%Y")[1:length(R0_numer)]
  
  R0_plot = ggplot(data = R0_df)+
    geom_line(mapping = aes(y = R0, x = dates1), size = 1 ) +
    geom_hline(yintercept = 1, col ='red') +
    labs(x = "Time", y = expression(R[0])) +
    theme_bw() +
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text = element_text(size = 20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title = element_text(size = 26))+
    theme(legend.position = "none", #c(.5,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text = element_text(size = 24))
  
  #R0_plot
  ggsave(file = sprintf("../Results for USA States/%s/%s_R0_analogue.pdf", st, st))
  
  
  ###A_hat
  
  Ahat_plot = ggplot(plot.data,aes(x=ddates,y=A.hat.t))+
    geom_pointline(size=.8)+
    #geom_line(y=lowess(est1$A.hat.t,f = 1/8)$y[1:(est1$n.t -1)])+
    geom_vline(xintercept = range(TimeSet)+10, linetype="dashed",
               color = "blue", size=.6)+
    labs(x = "Time", y = expression(hat(A[t])))+
    theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position= c(.5,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  
  #Ahat_plot
  ggsave(file =sprintf("../Results for USA States/%s/%s_Ahat.pdf", st, st))
  
  
  ###New Infections
  
  df_new.infect = as.data.frame(cbind( Time = plot.data$Time,
                                       new.infect.t = pmax(plot.data$new.infect.t,0),
                                       new.infect.sm = pmax(lowess(est1$new.infect.t,f = bw_lowess)$y,0),
                                       delC = pmax(var$delC[2:n.t],0),
                                       delC.sm = pmax(lowess(var$delC[2:n.t],f = bw_lowess)$y,0)))
  
  df_new.infect$dates1 = as.Date(dates[1:length(plot.data$Time)], format = "%m/%d/%Y")
  New_Infect = ggplot(data = df_new.infect)+
    geom_pointline(mapping = aes(x = dates1,y= new.infect.t, color = "Fitted"), size =.8 ) +
    #geom_point(aes(x = dates1,y= new.infect.t, color = "Fitted"), size = .7)+
    geom_pointline(aes(x=dates1,y= delC, color = "Observed" ), size = .8)+
    #geom_line(mapping=aes(y = new.infect.sm ,x= dates1, color = "Fitted"),size=.4 ) +
    #geom_line(mapping=aes(y = delC.sm,x= dates1,color = "Observed"),size=.4) +
    scale_color_manual(values = c('Fitted' = 'red',
                                  'Observed' = 'black')) +
    labs(color = 'Y series', x = "Time", y = 'New Infection')+
    theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position = "none", #c(.5,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  
  #New_Infect
  ggsave(file =sprintf("../Results for USA States/%s/%s_NewInfect.pdf", st, st))
  
  ##Doubling Rate
  
  df_doub = data.frame(
    doub_rt.obs = lowess(diff(log(df_CRHD$raw_Obs_C))[-n.t], f = bw_lowess)$y, 
    doub_rt.est = lowess(diff(log(plot.data$cum.new.infect.t)), f = bw_lowess)$y)
  df_doub$dates1 = as.Date(dates[1:(n.t-2)], format = "%m/%d/%Y")
  
  doub_rt_plot = ggplot(df_doub, aes(dates1)) + 
    geom_line(aes(y = doub_rt.obs , x = dates1, color = "Observed"),size= 1 ) +
    geom_line(aes(y = doub_rt.est , x = dates1, color = "Estimated"),size=1) +
    scale_color_manual(values = c(
      'Observed' = 'black',
      'Estimated' = 'red')) +
    labs(color = 'Y series', x = "Time", y = 'Doubling Rate')+
    theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position= "none",# c(.4,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  #doub_rt_plot
  ggsave(file =sprintf("../Results for USA States/%s/%s_Doubling_Rate.pdf", st, st))
  
  TtoH_plot = ggplot(df_CRHD, aes(dates1)) + 
    geom_pointline(aes(y = TtoH , x = dates1) , size=.8 ) +
    #geom_point(aes(y = TtoH , x = dates1) , size=.4 ) +
    labs( x = "Time", y = 'Testing per Hospitalization')+
    theme_bw()+
    scale_x_date(breaks = break.vec, #date_breaks("1 months"),
                 labels = date_format("%b %d")) +
    theme(axis.text=element_text(size=20), 
          axis.text.x = element_text(size = 20, angle = 90, hjust = 1),
          axis.title=element_text(size=26))+
    theme(legend.position= "none",# c(.4,.85),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=24))
  #TtoH_plot
  ggsave(file =sprintf("../Results for USA States/%s/%s_TtoH_plot.pdf", st, st))
  
}





#####################
#### Diagnostics ###

## Creates the diagnostics plots displaying the fits
diagnostics = function(whichBlock, TimeBlock , est,var){
  
  n.t = est$n.t
  j = whichBlock
  
  TimeSet = TimeBlock[[j]]
  
  phi.t = est$phi.hat.vec[j]#est$rhoH.hat.vec[j]
  gamma.t = est$gamma.hat
  rhoA.t = est$rhoA.hat
  rhoH.t = est$rhoH.hat.vec[j]
  
  respH.t = sqrt(pmax(var$delH + var$delD + var$delR,0))[1:n.t]
  predH.t = sqrt(pmax((gamma.t + rhoA.t) * lowess(var$Q[1:n.t],f = 1/8)$y + 
                        gamma.t * lowess(var$delC[1:n.t]/(gamma.t + phi.t * var$TtoH[1:n.t]), f=1/8)$y,0))
  
  respR.t = sqrt(pmax(var$delR,0))[1:n.t]
  predR.t = sqrt(pmax(rhoH.t * lowess(var$H[1:n.t],f=1/8)$y + 
                        rhoA.t * lowess(var$Q[1:n.t],f=1/8)$y,0))
  
  
  #quartz()
  diag.data <- as.data.frame(cbind(
    Time = 1:n.t,
    respH.t = respH.t,
    predH.t = predH.t,
    resH.t = respH.t - predH.t,
    respR.t = respR.t,
    predR.t = predR.t,
    resR.t = respR.t-predR.t
  ))
  
  p10 <- ggplot(diag.data,
                aes(x= Time,y=respH.t))+
    geom_point(size=0.7)+
    geom_line(y=lowess(respH.t,f = 1/8)$y)+
    geom_line(y=lowess(predH.t,f = 1/8)$y,color='red')+
    geom_vline(xintercept = range(TimeSet), linetype="dashed",
               color = "blue", size=.5)+
    labs(title =paste("Block ",j,": H-rate fits"))+
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))+
    theme(legend.position= c(.2,.9),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=10))
  
  p11 <- ggplot(diag.data,
                aes(x= Time,y=resH.t))+
    geom_point(size=0.7)+
    geom_line(aes(y=lowess(resH.t,f = 1/8)$y),color='red')+
    geom_vline(xintercept = range(TimeSet), linetype="dashed",
               color = "blue", size=.5)+      
    geom_hline(yintercept = 0, linetype="dashed",
               color = "blue", size=.5)+
    labs(title =paste("Block ",j,": H-rate residuals"),
         y='respH.t-predH.t')+
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))+
    theme(legend.position= c(.2,.9),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=10))
  
  p12 <- ggplot(diag.data,
                aes(x= Time,y=respR.t))+
    geom_point(size=0.7)+
    geom_line(aes(y=lowess(respR.t,f = 1/8)$y))+
    geom_line(aes(y=lowess(predR.t,f = 1/8)$y),color='red')+
    geom_vline(xintercept = range(TimeSet), linetype="dashed",
               color = "blue", size=.5)+
    labs(title =paste("Block ",j,": R-rate fits"))+theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))+
    theme(legend.position= c(.2,.9),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=10))
  
  p13 <- ggplot(diag.data,
                aes(x= Time,y=resR.t))+
    geom_point(size=0.7)+
    geom_line(aes(y=lowess(resR.t,f = 1/8)$y),color='red')+
    geom_vline(xintercept = range(TimeSet), linetype="dashed",
               color = "blue", size=.5)+      
    geom_hline(yintercept = 0, linetype="dashed",
               color = "blue", size=.5)+
    labs(title =paste("Block ",j,": R-rate residulas"),
         y='respR.t-predR.t')+
    theme_bw()+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12))+
    theme(legend.position= c(.2,.9),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.text=element_text(size=10))
  
  grid.arrange(p10,p11,p12,p13,nrow=2)+theme(plot.margin = unit(c(3,3,3,3), "lines"))
}


print(bw_lowess)