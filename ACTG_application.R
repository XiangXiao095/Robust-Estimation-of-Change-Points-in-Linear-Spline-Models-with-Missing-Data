library(speff2trial)
library(ggplot2)
library(dplyr)


#Data processing
data(ACTG175)
data.used<-ACTG175
trt1<-data.used$arms==1
trt2<-data.used$arms==2
trt3<-data.used$arms==3
Y_full<-data.used$cd496
X_full<-data.used$cd40
C_full<-1-is.na(Y_full)
data.used<-cbind(data.used,trt1,trt2,trt3)
Cov_full<-as.matrix(data.used[c("trt1","trt2","trt3","age","wtkg","race","gender","str2","offtrt")])
Aux_full<-as.matrix(data.used[c("cd420","cd80","cd820")])
fit.pro<-glm(C_full~.,data=as.data.frame(cbind(C_full,X_full,Cov_full,Aux_full)),family=binomial('logit'))
pi_full<-predict(fit.pro,type='response')
pi<-pi_full[which(C_full==1)]
wi_full<-1/pi_full
wi_full[which(C_full==0)]<-0
wi<-wi_full[which(C_full==1)]

Y<-Y_full[which(C_full==1)]
X<-X_full[which(C_full==1)]
Cov<-Cov_full[which(C_full==1),]
Aux<-Aux_full[which(C_full==1),]
I<-rep(1,length(X))
Y_full[which(C_full==0)]<-0  
tau_initial<-500



#Complete case
tau_cc<-tau_initial
eps<-10
k<-1

while (abs(eps)>0.1&k<=1000&!is.na(eps)){
  X_tau_cc<-(X-tau_cc)*as.numeric(X>tau_cc)
  data.cc<-cbind(Y,X,X_tau_cc,Cov)
  data.cc<-as.data.frame(data.cc)
  lin.fit_cc<-lm(Y~.,data=data.cc)
  coe_cc<-as.numeric(coef(lin.fit_cc))
  sig<-sigma(lin.fit_cc)
  constant_cc<-coe_cc[1]
  slope_cc<-coe_cc[2]
  theta_cc<-coe_cc[3]
  eta_cc<-coe_cc[-c(1:3)]
  y_pred_cc<-predict(lin.fit_cc)
  Un=-sum((Y-y_pred_cc)*as.numeric(X>tau_cc))
  Jn<-theta_cc*sum(as.numeric(X>tau_cc))
  eps<-1/Jn*Un
  tau_up<-tau_cc+eps
  tau_cc<-tau_up
  k<-k+1
}

Y_res_cc<-as.vector(Y-constant_cc-slope_cc*X-theta_cc*X_tau_cc-Cov%*%eta_cc)
U_derivative_cc<-rbind(I,X,X_tau_cc,-theta_cc*as.numeric(X>tau_cc),t(Cov))
U_derivative_cc_inv<-solve(U_derivative_cc%*%t(U_derivative_cc))

U_cc<-t(Y_res_cc*t(U_derivative_cc))
cov_cc<-U_derivative_cc_inv%*%U_cc
cov_cc<-cov_cc%*%t(cov_cc)



#IPW approach#
tau_ipw<-tau_initial
eps<-10
k<-1

while (abs(eps)>0.1&k<=1000&!is.na(eps)){
  X_tau_ipw<-(X-tau_ipw)*as.numeric(X>tau_ipw)
  data.ipw<-cbind(Y,X,X_tau_ipw,Cov)
  data.ipw<-as.data.frame(data.ipw)
  lin.fit_ipw<-lm(Y~.,weights=wi,data=data.ipw)
  coe_ipw<-as.numeric(coef(lin.fit_ipw))
  constant_ipw<-coe_ipw[1]
  slope_ipw<-coe_ipw[2]
  theta_ipw<-coe_ipw[3]
  eta_ipw<-coe_ipw[-c(1:3)]
  y_pred_ipw<-predict(lin.fit_ipw)
  Un=-sum(wi*(Y-y_pred_ipw)*as.numeric(X>tau_ipw))
  Jn<-theta_ipw*sum(wi*as.numeric(X>tau_ipw))
  eps<-1/Jn*Un
  tau_up<-tau_ipw+eps
  tau_ipw<-tau_up
  k<-k+1
}

X_full_indicator<-as.numeric(X_full>tau_ipw)
X_full_tau_ipw<-(X_full-tau_ipw)*X_full_indicator
H_ipw<-rbind(rep(1,length(X_full)),X_full,X_full_tau_ipw,-theta_ipw*X_full_indicator,t(Cov_full))
U_derivative_ipw<-H_ipw%*%(wi_full*t(H_ipw))
U_derivative_ipw_inv<-solve(U_derivative_ipw)
y_res_ipw_full<-as.vector(Y_full-constant_ipw-slope_ipw*X_full-theta_ipw*X_full_tau_ipw-Cov_full%*%eta_ipw)
U_ipw<-t(wi_full*y_res_ipw_full*t(H_ipw))
L_ipw<-rbind(rep(1,length(X_full)),X_full,t(Cov_full),t(Aux_full))
Q_ipw<-t((C_full-pi_full)*t(L_ipw))
U_ipw_alpha<-H_ipw%*%((1/pi_full-1)*C_full*y_res_ipw_full*t(L_ipw))
Q_ipw_alpha<-L_ipw%*%(pi_full*(1-pi_full)*t(L_ipw))
S_ipw<-U_ipw-U_ipw_alpha%*%solve(Q_ipw_alpha)%*%Q_ipw
cov_ipw<-U_derivative_ipw_inv%*%(S_ipw%*%t(S_ipw))%*%U_derivative_ipw_inv



##AIPW approach##
tau_aug<-tau_initial
eps<-10
k<-1

while (abs(eps)>0.1&k<=1000&!is.na(eps))
{
  X_tau_aug<-as.numeric(X>tau_aug)*(X-tau_aug)
  data.aug<-cbind(Y,X,X_tau_aug,Cov,Aux)
  data.aug<-as.data.frame(data.aug)
  lin.fit_aug<-lm(Y~.,data=data.aug)
  coe_aug<-as.numeric(coef(lin.fit_aug))
  y_pred_aug<-predict(lin.fit_aug)
  constant_aug<-coe_aug[1]
  slope_aug<-coe_aug[2]
  theta_aug<-coe_aug[3]
  eta_aug<-coe_aug[-c(1:3)]
  Un=-sum((Y-y_pred_aug)*as.numeric(X>tau_aug))
  Jn=theta_aug*sum(as.numeric(X>tau_aug))
  eps=1/Jn*Un
  tau_up<-tau_aug+eps
  tau_aug<-tau_up
  k<-k+1
}

y_full_aug<-as.vector(constant_aug+slope_aug*X_full+theta_aug*(X_full-tau_aug)*as.numeric(X_full>tau_aug)+cbind(Cov_full,Aux_full)%*%eta_aug)
tau_aipw<-tau_initial
eps<-10
k<-1

while (abs(eps)>0.1&k<=1000&!is.na(eps))
{
  Y_res_aug<-Y_full+(wi_full-1)*(Y_full-y_full_aug)
  X_tau_aipw_full<-(X_full-tau_aipw)*as.numeric(X_full>tau_aipw)
  data.aipw<-cbind(Y_res_aug,X_full,X_tau_aipw_full,Cov_full)
  data.aipw<-as.data.frame(data.aipw)
  lin_fit_aipw<-lm(Y_res_aug~.,data=data.aipw)
  coe_aipw<-as.numeric(coef(lin_fit_aipw))
  constant_aipw<-coe_aipw[1]
  slope_aipw<-coe_aipw[2]
  theta_aipw<-coe_aipw[3]
  eta_aipw<-coe_aipw[-c(1:3)]
  y_pred_aipw_full<-predict(lin_fit_aipw,data=data.aipw)
  Un<--sum((Y_res_aug-y_pred_aipw_full)*as.numeric(X_full>tau_aipw))
  Jn<-theta_aipw*sum(as.numeric(X_full>tau_aipw))
  eps<-1/Jn*Un
  tau_up<-tau_aipw+eps
  tau_aipw<-tau_up
  k<-k+1
}

X_full_indicator<-as.numeric(X_full>tau_aipw)
X_full_tau_aipw<-(X_full-tau_aipw)*X_full_indicator
H_aipw<-rbind(rep(1,length(X_full)),X_full,X_full_tau_aipw,-theta_aipw*X_full_indicator,t(Cov_full))
U_derivative_aipw<-H_aipw%*%t(H_aipw)
U_derivative_aipw_inv<-solve(U_derivative_aipw)
y_pred_aipw_full<-as.vector(constant_aipw+slope_aipw*X_full+theta_aipw*X_full_tau_aipw+Cov_full%*%eta_aipw)
y_res_aipw_full<-as.vector(Y_full-y_pred_aipw_full+(wi_full-1)*(Y_full-y_full_aug))
U_aipw<-t(y_res_aipw_full*t(H_aipw))
K_aipw<-rbind(rep(1,length(X_full)),X_full,(X_full-tau_aug)*as.numeric(X_full>tau_aug),-theta_aug*as.numeric(X_full>tau_aug),t(Cov_full),t(Aux_full))
L_aipw<-rbind(rep(1,length(X_full)),X_full,t(Cov_full),t(Aux_full))
R_aipw<-t(C_full*(Y_full-y_full_aug)*t(K_aipw))
Q_aipw<-t((C_full-pi_full)*t(L_aipw))
U_aipw_zeta<-H_aipw%*%((C_full/pi_full-1)*t(K_aipw))
R_aipw_zeta<-K_aipw%*%(C_full*t(K_aipw))
U_aipw_alpha<-H_aipw%*%((1/pi_full-1)*C_full*(Y_full-y_full_aug)*t(L_aipw))
Q_aipw_alpha<-L_aipw%*%(pi_full*(1-pi_full)*t(L_aipw))
S_aipw<-U_aipw-U_aipw_zeta%*%solve(R_aipw_zeta)%*%R_aipw-U_aipw_alpha%*%solve(Q_aipw_alpha)%*%Q_aipw
cov_aipw<-U_derivative_aipw_inv%*%(S_aipw%*%t(S_aipw))%*%U_derivative_aipw_inv



#Results
output_cc<-c(constant_cc,slope_cc,theta_cc,tau_cc,cov_cc[4,4])
output_ipw<-c(constant_ipw,slope_ipw,theta_ipw,tau_ipw,cov_ipw[4,4])
output_aipw<-c(constant_aipw,slope_aipw,theta_aipw,tau_aipw,cov_aipw[4,4])
output<-rbind(output_cc,output_ipw,output_aipw)
rownames(output)<-c("CC","IPW","DR-AIPW")
colnames(output)<-c("constant","slope","theta","tau","var est")
Cov_mean<-apply(Cov_full,2,mean)
Add<-crossprod(eta_aipw,Cov_mean)
data.plot<-data.frame(X,Y)



#Plot
method_colors<-c("DR-AIPW"="#0072B2","CC"="#E69F00")
method_fills<-method_colors
data.est<-rbind(
  data.frame(method="DR-AIPW",x=c(0,tau_aipw,900),
             y=c(constant_aipw+Add,
                 constant_aipw+slope_aipw*tau_aipw+Add,
                 constant_aipw+slope_aipw*900+theta_aipw*(900-tau_aipw)+Add)),
  data.frame(method="CC",x=c(0,tau_cc,900),
             y=c(constant_cc+Add,
                 constant_cc+slope_cc*tau_cc+Add,
                 constant_cc+slope_cc*900+theta_cc*(900-tau_cc)+Add))
) %>%
  group_by(method) %>%
  mutate(xend=x,yend=y,x=lag(x),y=lag(y)) %>%
  filter(!is.na(x)) %>% ungroup()

data.rib<-rbind(
  data.frame(method="DR-AIPW",X=c(tau_aipw-1.96*sqrt(cov_aipw[4,4]),
                                  tau_aipw+1.96*sqrt(cov_aipw[4,4]))),
  data.frame(method="CC",X=c(tau_cc-1.96*sqrt(cov_cc[4,4]),
                             tau_cc+1.96*sqrt(cov_cc[4,4])))
)

vline_data<-data.frame(method=c("DR-AIPW","CC"),xintercept=c(tau_aipw,tau_cc))

g<-ggplot()+
  geom_point(data=data.plot,aes(x=X,y=Y),size=1,shape=20,alpha=0.5)+
  coord_fixed(ratio=0.5)+
  theme_bw(base_size=12)+
  theme(panel.grid=element_blank(),legend.position="top")+
  scale_x_continuous(breaks=seq(0,900,by=100))+
  scale_y_continuous(breaks=seq(0,1200,by=200))+
  
  geom_segment(data=data.est,aes(x=x,y=y,xend=xend,yend=yend,color=method),linewidth=1)+
  
  geom_vline(data=vline_data,aes(xintercept=xintercept,color=method),linewidth=0.7)+
  
  geom_vline(xintercept=tau_aipw-1.96*sqrt(cov_aipw[4,4]),linetype=6,color=method_colors["DR-AIPW"],linewidth=0.4)+
  geom_vline(xintercept=tau_aipw+1.96*sqrt(cov_aipw[4,4]),linetype=6,color=method_colors["DR-AIPW"],linewidth=0.4)+
  geom_vline(xintercept=tau_cc-1.96*sqrt(cov_cc[4,4]),linetype=6,color=method_colors["CC"],linewidth=0.4)+
  geom_vline(xintercept=tau_cc+1.96*sqrt(cov_cc[4,4]),linetype=6,color=method_colors["CC"],linewidth=0.4)+

  geom_ribbon(data=subset(data.rib,method=="DR-AIPW"),
  aes(x=X,ymin=-Inf,ymax=Inf,fill=method),alpha=0.12)+
  geom_ribbon(data=subset(data.rib,method=="CC"),
  aes(x=X,ymin=-Inf,ymax=Inf,fill=method),alpha=0.12)+

  scale_color_manual(values=method_colors,name=NULL)+
  scale_fill_manual(values=method_fills,name=NULL)+
  xlab(expression(CD4[0]))+
  ylab(expression(CD4[96]))+
  theme(
  panel.grid=element_blank(),
  legend.position=c(0.95,0.95),
  legend.justification=c(1,1),
  legend.background=element_rect(fill=alpha("white",0.7),color=NA),
  legend.title=element_text(size=10),
  legend.text=element_text(size=9))

g

save.image("ACTG_application.Rdata")
