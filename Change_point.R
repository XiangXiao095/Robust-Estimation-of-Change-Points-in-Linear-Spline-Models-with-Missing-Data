library(parallel)
library(foreach)
library(doParallel)



Ncores<-max(1,detectCores()-2)
cl<-makeCluster(Ncores)
registerDoParallel(cl)


#tau stands for true value, while tau_int specifies initial point for two-step algorithm.
Changepoint<-function(seed,tau,beta,alpha0,tau_int,n)
{
  
  #Generate data#
  set.seed(20240319+seed)
  X_full<-runif(n,0,100)
  e_full<-rnorm(n,0,2)
  V_full<-rnorm(n,mean=-0.5*X_full,sd=2)
  Y_full<-30+0.125*X_full+beta*as.numeric(X_full>tau)*(X_full-tau)+0.75*V_full+e_full
  C_full<-rbinom(n,size=1,prob=1/(1+exp(-(alpha0+0.1*X_full+0.2*V_full+0.01*V_full^2))))
  data.full<-data.frame(Y_full,X_full,C_full,V_full)
  
  
  
  #Observed data#
  X<-data.full$X_full[which(C_full==1)]
  Y<-data.full$Y_full[which(C_full==1)]
  V<-data.full$V_full[which(C_full==1)]
  I_full<-rep(1,length(X_full))
  I<-rep(1,length(X))
  Y_full[which(C_full==0)]<-0
  
  
  
  #Correctly specified propensity score model#
  model.pi_C<-glm(C_full~X_full+V_full+I(V_full^2),family='binomial',data=data.full)
  pi_full_C<-predict(model.pi_C,type='response')
  wi_full_C<-1/pi_full_C
  wi_C<-wi_full_C[which(C_full==1)]
  wi_full_C[which(C_full==0)]<-0
  pi_C<-pi_full_C[which(C_full==1)]
  
  
  
  #Mis-specified missing probability model#
  model.pi_M<-glm(C_full~X_full+V_full,family='binomial',data=data.full)
  pi_full_M<-predict(model.pi_M,type='response')
  wi_full_M<-1/pi_full_M
  wi_M<-wi_full_M[which(C_full==1)]
  wi_full_M[which(C_full==0)]<-0
  pi_M<-pi_full_M[which(C_full==1)]
  
  
  
  #Correctly specified outcome regression model#
  tau_aug_C<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    X_aug_C_indicator<-X>tau_aug_C
    X_aug_C<-(X-tau_aug_C)*X_aug_C_indicator
    lin.fit_aug_C<-lm(Y~X+X_aug_C+V)
    theta_aug_C<-coef(lin.fit_aug_C)
    Y_res_aug_C<-Y-predict(lin.fit_aug_C)
    U<--sum(Y_res_aug_C*X_aug_C_indicator)
    J<-theta_aug_C[3]*sum(X_aug_C_indicator)
    eta<-U/J
    tau_aug_C<-tau_aug_C+eta
    j<-j+1
  }
  X_full_aug_C_indicator<-X_full>tau_aug_C
  X_full_aug_C<-(X_full-tau_aug_C)*X_full_aug_C_indicator
  Y_full_aug_C<-predict(lin.fit_aug_C,newdata=data.frame(X=X_full,X_aug_C=X_full_aug_C,V=V_full))
  
  
  
  #Mis-specified outcome regression model#
  tau_aug_M<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    X_aug_M_indicator<-X>tau_aug_M
    X_aug_M<-(X-tau_aug_M)*X_aug_M_indicator
    lin.fit_aug_M<-lm(Y~X_aug_M+V)
    theta_aug_M<-coef(lin.fit_aug_M)
    Y_res_aug_M<-Y-predict(lin.fit_aug_M)
    U<--sum(Y_res_aug_M*X_aug_M_indicator)
    J<-theta_aug_M[2]*sum(X_aug_M_indicator)
    eta<-U/J
    tau_aug_M<-tau_aug_M+eta
    j<-j+1
  }
  X_full_aug_M_indicator<-X_full>tau_aug_M
  X_full_aug_M<-(X_full-tau_aug_M)*X_full_aug_M_indicator
  Y_full_aug_M<-predict(lin.fit_aug_M,newdata=data.frame(X_aug_M=X_full_aug_M,V=V_full))
  
  
  
  ##Complete case##
  tau_cc<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    X_cc_indicator<-X>tau_cc
    X_cc<-(X-tau_cc)*X_cc_indicator
    lin.fit_cc<-lm(Y~X+X_cc)
    theta_cc<-coef(lin.fit_cc)
    Y_res_cc<-Y-predict(lin.fit_cc)
    U<--sum(Y_res_cc*X_cc_indicator)
    J<-theta_cc[3]*sum(X_cc_indicator)
    eta<-U/J
    tau_cc<-tau_cc+eta
    j<-j+1
  }
  H_cc<-rbind(I,X,X_cc,-theta_cc[3]*X_cc_indicator)
  U_derivative_cc_inv<-solve(H_cc%*%t(H_cc))
  U_cc<-t(Y_res_cc*t(H_cc))
  cov_cc<-U_derivative_cc_inv%*%U_cc
  cov_cc<-cov_cc%*%t(cov_cc)
  cover_cc<-((tau_cc-1.96*sqrt(cov_cc[4,4]))<=tau&tau<=(tau_cc+1.96*sqrt(cov_cc[4,4])))
  
  
  
  ##IPW(M)##
  tau_ipw_M<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    X_ipw_M_indicator<-X>tau_ipw_M
    X_ipw_M<-(X-tau_ipw_M)*X_ipw_M_indicator
    lin.fit_ipw_M<-lm(Y~X+X_ipw_M,weights=wi_M)
    theta_ipw_M<-coef(lin.fit_ipw_M)
    Y_res_ipw_M<-Y-predict(lin.fit_ipw_M)
    U<--sum(wi_M*Y_res_ipw_M*X_ipw_M_indicator)
    J<-theta_ipw_M[3]*sum(wi_M*X_ipw_M_indicator)
    eta<-U/J
    tau_ipw_M<-tau_ipw_M+eta
    j<-j+1
  }
  X_full_ipw_M_indicator<-X_full>tau_ipw_M
  X_full_ipw_M<-(X_full-tau_ipw_M)*X_full_ipw_M_indicator
  H_ipw_M<-rbind(I_full,X_full,X_full_ipw_M,-theta_ipw_M[3]*X_full_ipw_M_indicator)
  U_derivative_ipw_M<-H_ipw_M%*%(wi_full_M*t(H_ipw_M))
  U_derivative_ipw_M_inv<-solve(U_derivative_ipw_M)
  Y_full_res_ipw_M<-Y_full-predict(lin.fit_ipw_M,newdata=data.frame(X=X_full,X_ipw_M=X_full_ipw_M))
  U_ipw_M<-t(wi_full_M*Y_full_res_ipw_M*t(H_ipw_M))
  L_ipw_M<-rbind(I_full,X_full,V_full)
  Q_ipw_M<-t((C_full-pi_full_M)*t(L_ipw_M))
  U_ipw_M_alpha<-H_ipw_M%*%((1/pi_full_M-1)*C_full*Y_full_res_ipw_M*t(L_ipw_M))
  Q_ipw_M_alpha<-L_ipw_M%*%(pi_full_M*(1-pi_full_M)*t(L_ipw_M))
  S_ipw_M<-U_ipw_M-U_ipw_M_alpha%*%solve(Q_ipw_M_alpha)%*%Q_ipw_M
  cov_ipw_M<-U_derivative_ipw_M_inv%*%S_ipw_M
  cov_ipw_M<-cov_ipw_M%*%t(cov_ipw_M)
  cover_ipw_M<-((tau_ipw_M-1.96*sqrt(cov_ipw_M[4,4]))<=tau&tau<=(tau_ipw_M+1.96*sqrt(cov_ipw_M[4,4])))

  
  
  ##IPW(C)##
  tau_ipw_C<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    X_ipw_C_indicator<-X>tau_ipw_C
    X_ipw_C<-(X-tau_ipw_C)*X_ipw_C_indicator
    lin.fit_ipw_C<-lm(Y~X+X_ipw_C,weights=wi_C)
    theta_ipw_C<-coef(lin.fit_ipw_C)
    Y_res_ipw_C<-Y-predict(lin.fit_ipw_C)
    U<--sum(wi_C*Y_res_ipw_C*X_ipw_C_indicator)
    J<-theta_ipw_C[3]*sum(wi_C*X_ipw_C_indicator)
    eta<-U/J
    tau_ipw_C<-tau_ipw_C+eta
    j<-j+1
  }
  X_full_ipw_C_indicator<-X_full>tau_ipw_C
  X_full_ipw_C<-(X_full-tau_ipw_C)*X_full_ipw_C_indicator
  H_ipw_C<-rbind(I_full,X_full,X_full_ipw_C,-theta_ipw_C[3]*X_full_ipw_C_indicator)
  U_derivative_ipw_C<-H_ipw_C%*%(wi_full_C*t(H_ipw_C))
  U_derivative_ipw_C_inv<-solve(U_derivative_ipw_C)
  Y_full_res_ipw_C<-Y_full-predict(lin.fit_ipw_C,newdata=data.frame(X=X_full,X_ipw_C=X_full_ipw_C))
  U_ipw_C<-t(wi_full_C*Y_full_res_ipw_C*t(H_ipw_C))
  L_ipw_C<-rbind(I_full,X_full,V_full,V_full^2)
  Q_ipw_C<-t((C_full-pi_full_C)*t(L_ipw_C))
  U_ipw_C_alpha<-H_ipw_C%*%((1/pi_full_C-1)*C_full*Y_full_res_ipw_C*t(L_ipw_C))
  Q_ipw_C_alpha<-L_ipw_C%*%(pi_full_C*(1-pi_full_C)*t(L_ipw_C))
  S_ipw_C<-U_ipw_C-U_ipw_C_alpha%*%solve(Q_ipw_C_alpha)%*%Q_ipw_C
  cov_ipw_C<-U_derivative_ipw_C_inv%*%S_ipw_C
  cov_ipw_C<-cov_ipw_C%*%t(cov_ipw_C)
  cover_ipw_C<-((tau_ipw_C-1.96*sqrt(cov_ipw_C[4,4]))<=tau&tau<=(tau_ipw_C+1.96*sqrt(cov_ipw_C[4,4])))
  
  
  
  ##AIPW(MM)##
  tau_aipw_MM<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    Y_full_res_aug_MM<-Y_full+(wi_full_M-1)*(Y_full-Y_full_aug_M)
    X_full_aipw_MM_indicator<-X_full>tau_aipw_MM
    X_full_aipw_MM<-(X_full-tau_aipw_MM)*X_full_aipw_MM_indicator
    lin_fit_aipw_MM<-lm(Y_full_res_aug_MM~X_full+X_full_aipw_MM)
    theta_aipw_MM<-coef(lin_fit_aipw_MM)
    Y_full_res_aipw_MM<-Y_full_res_aug_MM-predict(lin_fit_aipw_MM)
    U<--sum(Y_full_res_aipw_MM*X_full_aipw_MM_indicator)
    J<-theta_aipw_MM[3]*sum(X_full_aipw_MM_indicator)
    eta<-U/J
    tau_aipw_MM<-tau_aipw_MM+eta
    j<-j+1
  }
  H_aipw_MM<-rbind(I_full,X_full,X_full_aipw_MM,-theta_aipw_MM[3]*X_full_aipw_MM_indicator)
  U_derivative_aipw_MM<-H_aipw_MM%*%t(H_aipw_MM)
  U_derivative_aipw_MM_inv<-solve(U_derivative_aipw_MM)
  U_aipw_MM<-t(Y_full_res_aipw_MM*t(H_aipw_MM))
  K_aipw_MM<-rbind(I_full,X_full_aug_M,-theta_aug_M[2]*X_full_aug_M_indicator,V_full)
  L_aipw_MM<-rbind(I_full,X_full,V_full)
  R_aipw_MM<-t(C_full*(Y_full-Y_full_aug_M)*t(K_aipw_MM))
  Q_aipw_MM<-t((C_full-pi_full_M)*t(L_aipw_MM))
  U_aipw_MM_zeta<-H_aipw_MM%*%((C_full/pi_full_M-1)*t(K_aipw_MM))
  R_aipw_MM_zeta<-K_aipw_MM%*%(C_full*t(K_aipw_MM))
  U_aipw_MM_alpha<-H_aipw_MM%*%((1/pi_full_M-1)*C_full*(Y_full-Y_full_aug_M)*t(L_aipw_MM))
  Q_aipw_MM_alpha<-L_aipw_MM%*%(pi_full_M*(1-pi_full_M)*t(L_aipw_MM))
  S_aipw_MM<-U_aipw_MM-U_aipw_MM_zeta%*%solve(R_aipw_MM_zeta)%*%R_aipw_MM-U_aipw_MM_alpha%*%solve(Q_aipw_MM_alpha)%*%Q_aipw_MM
  cov_aipw_MM<-U_derivative_aipw_MM_inv%*%S_aipw_MM
  cov_aipw_MM<-cov_aipw_MM%*%t(cov_aipw_MM)
  cover_aipw_MM<-((tau_aipw_MM-1.96*sqrt(cov_aipw_MM[4,4]))<=tau&tau<=(tau_aipw_MM+1.96*sqrt(cov_aipw_MM[4,4])))
  
  
  
  ##AIPW(MC)##
  tau_aipw_MC<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    Y_full_res_aug_MC<-Y_full+(wi_full_M-1)*(Y_full-Y_full_aug_C)
    X_full_aipw_MC_indicator<-X_full>tau_aipw_MC
    X_full_aipw_MC<-(X_full-tau_aipw_MC)*X_full_aipw_MC_indicator
    lin_fit_aipw_MC<-lm(Y_full_res_aug_MC~X_full+X_full_aipw_MC)
    theta_aipw_MC<-coef(lin_fit_aipw_MC)
    Y_full_res_aipw_MC<-Y_full_res_aug_MC-predict(lin_fit_aipw_MC)
    U<--sum(Y_full_res_aipw_MC*X_full_aipw_MC_indicator)
    J<-theta_aipw_MC[3]*sum(X_full_aipw_MC_indicator)
    eta<-U/J
    tau_aipw_MC<-tau_aipw_MC+eta
    j<-j+1
  }
  H_aipw_MC<-rbind(I_full,X_full,X_full_aipw_MC,-theta_aipw_MC[3]*X_full_aipw_MC_indicator)
  U_derivative_aipw_MC<-H_aipw_MC%*%t(H_aipw_MC)
  U_derivative_aipw_MC_inv<-solve(U_derivative_aipw_MC)
  U_aipw_MC<-t(Y_full_res_aipw_MC*t(H_aipw_MC))
  K_aipw_MC<-rbind(I_full,X_full,X_full_aug_C,-theta_aug_C[3]*X_full_aug_C_indicator,V_full)
  L_aipw_MC<-rbind(I_full,X_full,V_full)
  R_aipw_MC<-t(C_full*(Y_full-Y_full_aug_C)*t(K_aipw_MC))
  Q_aipw_MC<-t((C_full-pi_full_M)*t(L_aipw_MC))
  U_aipw_MC_zeta<-H_aipw_MC%*%((C_full/pi_full_M-1)*t(K_aipw_MC))
  R_aipw_MC_zeta<-K_aipw_MC%*%(C_full*t(K_aipw_MC))
  U_aipw_MC_alpha<-H_aipw_MC%*%((1/pi_full_M-1)*C_full*(Y_full-Y_full_aug_C)*t(L_aipw_MC))
  Q_aipw_MC_alpha<-L_aipw_MC%*%(pi_full_M*(1-pi_full_M)*t(L_aipw_MC))
  S_aipw_MC<-U_aipw_MC-U_aipw_MC_zeta%*%solve(R_aipw_MC_zeta)%*%R_aipw_MC-U_aipw_MC_alpha%*%solve(Q_aipw_MC_alpha)%*%Q_aipw_MC
  cov_aipw_MC<-U_derivative_aipw_MC_inv%*%S_aipw_MC
  cov_aipw_MC<-cov_aipw_MC%*%t(cov_aipw_MC)
  cover_aipw_MC<-((tau_aipw_MC-1.96*sqrt(cov_aipw_MC[4,4]))<=tau&tau<=(tau_aipw_MC+1.96*sqrt(cov_aipw_MC[4,4])))
  
  
  
  ##AIPW(CM)##
  tau_aipw_CM<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    Y_full_res_aug_CM<-Y_full+(wi_full_C-1)*(Y_full-Y_full_aug_M)
    X_full_aipw_CM_indicator<-X_full>tau_aipw_CM
    X_full_aipw_CM<-(X_full-tau_aipw_CM)*X_full_aipw_CM_indicator
    lin_fit_aipw_CM<-lm(Y_full_res_aug_CM~X_full+X_full_aipw_CM)
    theta_aipw_CM<-coef(lin_fit_aipw_CM)
    Y_full_res_aipw_CM<-Y_full_res_aug_CM-predict(lin_fit_aipw_CM)
    U<--sum(Y_full_res_aipw_CM*X_full_aipw_CM_indicator)
    J<-theta_aipw_CM[3]*sum(X_full_aipw_CM_indicator)
    eta<-U/J
    tau_aipw_CM<-tau_aipw_CM+eta
    j<-j+1
  }
  H_aipw_CM<-rbind(I_full,X_full,X_full_aipw_CM,-theta_aipw_CM[3]*X_full_aipw_CM_indicator)
  U_derivative_aipw_CM<-H_aipw_CM%*%t(H_aipw_CM)
  U_derivative_aipw_CM_inv<-solve(U_derivative_aipw_CM)
  U_aipw_CM<-t(Y_full_res_aipw_CM*t(H_aipw_CM))
  K_aipw_CM<-rbind(I_full,X_full_aug_M,-theta_aug_M[2]*X_full_aug_M_indicator,V_full)
  L_aipw_CM<-rbind(I_full,X_full,V_full,V_full^2)
  R_aipw_CM<-t(C_full*(Y_full-Y_full_aug_M)*t(K_aipw_CM))
  Q_aipw_CM<-t((C_full-pi_full_C)*t(L_aipw_CM))
  U_aipw_CM_zeta<-H_aipw_CM%*%((C_full/pi_full_C-1)*t(K_aipw_CM))
  R_aipw_CM_zeta<-K_aipw_CM%*%(C_full*t(K_aipw_CM))
  U_aipw_CM_alpha<-H_aipw_CM%*%((1/pi_full_C-1)*C_full*(Y_full-Y_full_aug_M)*t(L_aipw_CM))
  Q_aipw_CM_alpha<-L_aipw_CM%*%(pi_full_C*(1-pi_full_C)*t(L_aipw_CM))
  S_aipw_CM<-U_aipw_CM-U_aipw_CM_zeta%*%solve(R_aipw_CM_zeta)%*%R_aipw_CM-U_aipw_CM_alpha%*%solve(Q_aipw_CM_alpha)%*%Q_aipw_CM
  cov_aipw_CM<-U_derivative_aipw_CM_inv%*%S_aipw_CM
  cov_aipw_CM<-cov_aipw_CM%*%t(cov_aipw_CM)
  cover_aipw_CM<-((tau_aipw_CM-1.96*sqrt(cov_aipw_CM[4,4]))<=tau&tau<=(tau_aipw_CM+1.96*sqrt(cov_aipw_CM[4,4])))
  
  
  
  ##AIPW(CC)##
  tau_aipw_CC<-tau_int
  eta<-1
  j<-1
  while (abs(eta)>1e-04&j<=300&!is.na(eta))
  {
    Y_full_res_aug_CC<-Y_full+(wi_full_C-1)*(Y_full-Y_full_aug_C)
    X_full_aipw_CC_indicator<-X_full>tau_aipw_CC
    X_full_aipw_CC<-(X_full-tau_aipw_CC)*X_full_aipw_CC_indicator
    lin_fit_aipw_CC<-lm(Y_full_res_aug_CC~X_full+X_full_aipw_CC)
    theta_aipw_CC<-coef(lin_fit_aipw_CC)
    Y_full_res_aipw_CC<-Y_full_res_aug_CC-predict(lin_fit_aipw_CC)
    U<--sum(Y_full_res_aipw_CC*X_full_aipw_CC_indicator)
    J<-theta_aipw_CC[3]*sum(X_full_aipw_CC_indicator)
    eta<-U/J
    tau_aipw_CC<-tau_aipw_CC+eta
    j<-j+1
  }
  H_aipw_CC<-rbind(I_full,X_full,X_full_aipw_CC,-theta_aipw_CC[3]*X_full_aipw_CC_indicator)
  U_derivative_aipw_CC<-H_aipw_CC%*%t(H_aipw_CC)
  U_derivative_aipw_CC_inv<-solve(U_derivative_aipw_CC)
  U_aipw_CC<-t(Y_full_res_aipw_CC*t(H_aipw_CC))
  K_aipw_CC<-rbind(I_full,X_full,X_full_aug_C,-theta_aug_C[3]*X_full_aug_C_indicator,V_full)
  L_aipw_CC<-rbind(I_full,X_full,V_full,V_full^2)
  R_aipw_CC<-t(C_full*(Y_full-Y_full_aug_C)*t(K_aipw_CC))
  Q_aipw_CC<-t((C_full-pi_full_C)*t(L_aipw_CC))
  U_aipw_CC_zeta<-H_aipw_CC%*%((C_full/pi_full_C-1)*t(K_aipw_CC))
  R_aipw_CC_zeta<-K_aipw_CC%*%(C_full*t(K_aipw_CC))
  U_aipw_CC_alpha<-H_aipw_CC%*%((1/pi_full_C-1)*C_full*(Y_full-Y_full_aug_C)*t(L_aipw_CC))
  Q_aipw_CC_alpha<-L_aipw_CC%*%(pi_full_C*(1-pi_full_C)*t(L_aipw_CC))
  S_aipw_CC<-U_aipw_CC-U_aipw_CC_zeta%*%solve(R_aipw_CC_zeta)%*%R_aipw_CC-U_aipw_CC_alpha%*%solve(Q_aipw_CC_alpha)%*%Q_aipw_CC
  cov_aipw_CC<-U_derivative_aipw_CC_inv%*%S_aipw_CC
  cov_aipw_CC<-cov_aipw_CC%*%t(cov_aipw_CC)
  cover_aipw_CC<-((tau_aipw_CC-1.96*sqrt(cov_aipw_CC[4,4]))<=tau&tau<=(tau_aipw_CC+1.96*sqrt(cov_aipw_CC[4,4])))
  
  
  
  ##results##
  tau_est<-c(tau_cc,tau_ipw_M,tau_ipw_C,tau_aipw_MM,tau_aipw_MC,tau_aipw_CM,tau_aipw_CC)-tau
  cov_est<-c(cov_cc[4,4],cov_ipw_M[4,4],cov_ipw_C[4,4],cov_aipw_MM[4,4],cov_aipw_MC[4,4],cov_aipw_CM[4,4],cov_aipw_CC[4,4])
  cover<-c(cover_cc,cover_ipw_M,cover_ipw_C,cover_aipw_MM,cover_aipw_MC,cover_aipw_CM,cover_aipw_CC)
  result_matrix<-rbind(tau_est,cov_est,cover)
  
  return(result_matrix)}



#Replications
seed.values<-c(1:1000)
results_1<-foreach(seed=seed.values) %dopar% {
  tryCatch({result<-Changepoint(seed,30,0.5,-1.25,30,500)
  if (!is.matrix(result)||anyNA(result)){
    cat("Invalid result (NA or not a matrix) in task:",seed,"\n")
    return(NULL)}
  return(result)},error=function(e){
    cat("Error in task:",seed,"\n")
    cat("Error message:",conditionMessage(e),"\n")
    return(NULL)
  })
}
#Leads to 20% missing rate. Set alpha0 to -2.5 for 30% missing rate.
results_2<-foreach(seed=seed.values) %dopar% {
  tryCatch({result<-Changepoint(seed,30,0.5,-1.25,30,1000)
  if (!is.matrix(result)||anyNA(result)){
    cat("Invalid result (NA or not a matrix) in task:",seed,"\n")
    return(NULL)}
  return(result)},error=function(e){
    cat("Error in task:",seed,"\n")
    cat("Error message:",conditionMessage(e),"\n")
    return(NULL)
  })
}
results_3<-foreach(seed=seed.values) %dopar% {
  tryCatch({result<-Changepoint(seed,30,0.5,-1.25,30,5000)
  if (!is.matrix(result)||anyNA(result)){
    cat("Invalid result (NA or not a matrix) in task:",seed,"\n")
    return(NULL)}
  return(result)},error=function(e){
    cat("Error in task:",seed,"\n")
    cat("Error message:",conditionMessage(e),"\n")
    return(NULL)
  })
}
results_4<-foreach(seed=seed.values) %dopar% {
  tryCatch({result<-Changepoint(seed,30,-0.5,-1.25,30,500)
  if (!is.matrix(result)||anyNA(result)){
    cat("Invalid result (NA or not a matrix) in task:",seed,"\n")
    return(NULL)}
  return(result)},error=function(e){
    cat("Error in task:",seed,"\n")
    cat("Error message:",conditionMessage(e),"\n")
    return(NULL)
  })
}
results_5<-foreach(seed=seed.values) %dopar% {
  tryCatch({result<-Changepoint(seed,30,-0.5,-1.25,30,1000)
  if (!is.matrix(result)||anyNA(result)){
    cat("Invalid result (NA or not a matrix) in task:",seed,"\n")
    return(NULL)}
  return(result)},error=function(e){
    cat("Error in task:",seed,"\n")
    cat("Error message:",conditionMessage(e),"\n")
    return(NULL)
  })
}
results_6<-foreach(seed=seed.values) %dopar% {
  tryCatch({result<-Changepoint(seed,30,-0.5,-1.25,30,5000)
  if (!is.matrix(result)||anyNA(result)){
    cat("Invalid result (NA or not a matrix) in task:",seed,"\n")
    return(NULL)}
  return(result)},error=function(e){
    cat("Error in task:",seed,"\n")
    cat("Error message:",conditionMessage(e),"\n")
    return(NULL)
  })
}
stopCluster(cl)



#Result processing
results_list<-list(results_1,results_2,results_3,results_4,results_5,results_6)

process_result<-function(results) {
  results<-results[!sapply(results,is.null)]
  results_array<-simplify2array(results)
  results_mean<-apply(results_array,c(1,2),mean)
  montecarlo_var<-apply(as.matrix(results_array[1,,]),1,var)
  output<-rbind(results_mean,montecarlo_var)
  output<-t(output)
  output<-cbind(output[,c(1,2,4)],output[,3])
  rownames(output)<-c("CC","IPW(M)","IPW(C)","AIPW(MM)","AIPW(MC)","AIPW(CM)","AIPW_(CC)")
  colnames(output)<-c("bias","estimated variance","monte_carlo variance","cover rate")
  return(100*output)
}

outputs_list<-lapply(results_list,process_result)
for (i in seq_along(outputs_list)) {assign(paste0("output_",i),outputs_list[[i]])}
output_all<-rbind(cbind(output_1,output_2,output_3),cbind(output_4,output_5,output_6))


save.image("Change_point.Rdata")