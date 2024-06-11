#本代码计算实用不同模型拟合、按不同模型计算的似然值



##################
#####准备部分#####
##################

#导入函数文件
source('D:\\paper\\2CRMDH\\Code_mine\\Empirical\\Likelihood_functions.R')
library(beepr)
library(mfGARCH)

###估计初始值及上下界###
#GARCH
inig<- c(0,0.1,0.1,0.2)
lowerg<-c(-Inf,0,0,0)
upperg<-c(Inf,Inf,Inf,Inf)

#RGARCH
ini0<- c(0,0.01,0.1,0.2,0.000,0.000,0.00,0.1,1)
lower0<-c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0)
upper0<-c(Inf, Inf ,Inf ,Inf ,Inf ,Inf ,Inf ,Inf,Inf)

#RMDH
ini1_bound<- c(0,1e-2,0.1,0.2,1e-3,0,0,0,1,1,0.5,0.1,0.8)
lower1_bound<-c(-Inf,0,0,0,0,-Inf,-Inf,-Inf,0,0,0,0,0)
upper1_bound<-c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)

#2CRMDH
ini3<- c(0,1e-2,1e-2,0.05,0.05,0.05,0.05,0.1,0.1,1e-3,1e-3,0,0,0,1,1,1,0.5,0.1,1)
lower3<-c(-Inf,0,0,0,0,0,0,0,0,0,0,-Inf,-Inf,-Inf,0,0,0,0,0,0)
upper3<-c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)

#GARCH-MIDAS
inigm<- c(0,0.1,0.8,4,0.1,0,1.5)
lowergm<-c(-Inf,-Inf,-Inf,-Inf,-Inf,1,1)
uppergm<-c(Inf,Inf,Inf,Inf,Inf,300,300)

#合并初始参数
ini_list<-list(inig,ini0,ini1_bound,ini3,inigm)
lower_list<-list(lowerg,lower0,lower1_bound,lower3,lowergm)
upper_list<-list(upperg,upper0,upper1_bound,upper3,uppergm)

#文件名
name_list<-c("000001","000002","000651","000858","600000","600016","600019","600028",
             "600104","600276","600519","600690","AAPL","BAC","GE",
             "IBM","INTC","JPM","KO","MCD","MSFT","WMT")

#合并函数名
llf_list<-list(llfg,llf0,llf1,llf3,llfgm)
get_llf_list<-list(geth_llfg,geth_llf0,geth_llf1,geth_llf3,geth_llfgm)
pllf_list<-list(pllfg,pllf0)



################
######估计######
################

#各模型的参数长度
s_list<-c(4,9,13,20,7)
K=20

#储存结果
FitResult<-matrix(0,length(name_list),7) #7个似然
FitParameter<-matrix(0,length(name_list)*length(llf_list),max(s_list)) #储存模型拟合参数
# date_mat=matrix(0,length(name_list),2)

#正式估计
time1=Sys.time()
for (q in 1:length(name_list)){
  
  #读取数据
  symbol <- name_list[q]
  last_namr<-paste(symbol,"_otc_RV_5min.csv",sep="")
  family_name <-paste("D:\\paper\\2CRMDH\\Peng\\5min_regdata_otc_final\\",last_namr,sep="")
  data<-read.csv(family_name)[,3:5]
  data_RV = data$V3
  t<-nrow(data)
  
  # date_mat[q,1]=read.csv(family_name)[1,2]
  # date_mat[q,2]=read.csv(family_name)[t,2]
  
  # new_data<-read.csv(family_name)[,2:3]
  # new_data$V1=as.Date(as.character(new_data$V1),format="%Y%m%d")
  # colnames(new_data)[1]='date'
  # new_data$V2=new_data$V2*100
  
  #####使用GARCH模型拟合，计算GARCH模型似然#####
  i=1
  
  #GARCH模型拟合出的参数
  fit_RV<-nlminb(ini_list[[i]], llf_list[[i]], data=data, lower = lower_list[[i]], upper = upper_list[[i]])
  theta<-fit_RV[["par"]]
  FitParameter[5*(q-1)+i,1:s_list[i]]<-theta
  
  #传递参数信息
  h = get_llf_list[[i]](theta,data) #计算在此参数下的波动率
  mu = theta[1]
  
  #计算GARCH似然
  QL = pllfg(mu,data,h)
  FitResult[q,1]=QL
  
  
  #####使用RCARCH模型拟合，计算GARCH模型和RGARCH模型似然#####
  i=2
  
  #RCARCH模型拟合出的参数
  fit_RV<-nlminb(ini_list[[i]], llf_list[[i]], data=data, lower = lower_list[[i]], upper = upper_list[[i]])
  theta<-fit_RV[["par"]]
  FitParameter[5*(q-1)+i,1:s_list[i]]<-theta
  
  #传递参数信息
  h = get_llf_list[[i]](theta,data)
  mu = theta[1]
  param = numeric(6)
  param[1]=theta[1] #mu
  param[2]=theta[5] #xi
  param[3]=theta[9] #phi
  param[4:6]=theta[6:8] #tau1,tau2,sigmau
  
  #计算GARCH似然
  QL = pllfg(mu,data,h)
  FitResult[q,2]=QL
  
  #计算RGARCH似然
  QL = pllf0(param,data,h)
  FitResult[q,5]=QL
  
  
  #####使用RMDH模型拟合，计算GARCH模型和RGARCH模型似然#####
  i=3
  
  #RMDH模型拟合出的参数,control=list(trace=TRUE)
  fit_RV<-nlminb(ini_list[[i]], llf_list[[i]], data=data, lower = lower_list[[i]], upper = upper_list[[i]])
  theta<-fit_RV[["par"]]
  FitParameter[5*(q-1)+i,1:s_list[i]]<-theta
  
  #传递参数信息
  h = get_llf_list[[i]](theta,data)
  mu = theta[1]
  param=numeric(6)
  param[1]=theta[1] #mu
  param[2]=theta[6] #xi
  param[3]=theta[13] #phi
  param[4:5]=theta[7:8] #tau1和tau2
  param[6]=theta[12] #sigmau
  
  #计算GARCH似然
  QL = pllfg(mu,data,h)
  FitResult[q,3]=QL
  
  #计算RGARCH似然
  QL = pllf0(param,data,h)
  FitResult[q,6]=QL
  
  
  #####使用2CRMDH模型拟合，计算GARCH模型和RGARCH模型似然#####
  i=4
  
  #2CRMDH模型拟合出的参数
  fit_RV<-nlminb(ini_list[[i]], llf_list[[i]], data=data, lower = lower_list[[i]], upper = upper_list[[i]])
  theta<-fit_RV[["par"]]
  FitParameter[5*(q-1)+i,1:s_list[i]]<-theta
  
  #传递参数信息
  h = get_llf_list[[i]](theta,data)
  mu = theta[1]
  param=numeric(6)
  param[1]=theta[1] #mu
  param[2]=theta[12] #xi
  param[3]=theta[20] #phi
  param[4:5]=theta[13:14] #tau1和tau2
  param[6]=theta[19] #sigmau
  
  #计算GARCH似然
  QL = pllfg(mu,data,h)
  FitResult[q,4]=QL
  
  #计算RGARCH似然
  QL = pllf0(param,data,h)
  FitResult[q,7]=QL
  
  
  #####使用GARCH-MIDAS模型拟合，计算GARCH模型和RGARCH模型似然#####
  i=5
  
  #GARCH-MIDAS模型拟合出的参数
  fit_RV<-nlminb(ini_list[[i]], llf_list[[i]], data=data, lower = lower_list[[i]], upper = upper_list[[i]])
  theta<-fit_RV[["par"]]
  FitParameter[5*(q-1)+i,1:s_list[i]]<-theta
  
  #传递参数信息
  h = get_llf_list[[i]](theta,data,data_RV,index0,K)
  mu = theta[1]
  
  #计算GARCH似然
  QL = pllfg(mu,data,h)
  FitResult[q,4]=QL

  
  
  #####输出循环信息#####
  print(symbol)
  
}
time2=Sys.time()

time2-time1
beep(sound=8)

#储存估计结果
write.csv(FitParameter,"D:\\paper\\2CRMDH\\Result\\Empirical\\Para_result_7_II.csv")
write.csv(FitResult,"D:\\paper\\2CRMDH\\Result\\Empirical\\likelihood_result_7_II.csv")

##############
#####补丁#####
##############
###根据上述储存结果计算论文中的Table1###

#读取数据
ll_mat<-as.matrix(read.csv("D:\\paper\\2CRMDH\\Result\\Empirical\\likelihood_result_7_II.csv",header=TRUE)[,-1])

#文件名
name_list<-c("000001","000002","000651","000858","600000","600016","600019","600028",
             "600104","600276","600519","600690","AAPL","BAC","GE",
             "IBM","INTC","JPM","KO","MCD","MSFT","WMT")

sv_mat=matrix(0,nrow=(length(name_list)+2),ncol=7) #7个似然

#中国组数据
for (ch_id in 1:12){
  
  #完全数据录入
  sv_mat[ch_id,]=ll_mat[ch_id,]
  
  #秩计算
  sv_mat[13,1:4]=sv_mat[13,1:4]+rank(-ll_mat[ch_id,1:4])/12
  sv_mat[13,5:7]=sv_mat[13,5:7]+rank(-ll_mat[ch_id,5:7])/12
  
}

#美国组数据
for (am_id in 13:22){
  
  #完全数据录入
  sv_mat[(am_id+1),]=ll_mat[am_id,]
  
  #秩计算
  sv_mat[24,1:4]=sv_mat[24,1:4]+rank(-ll_mat[am_id,1:4])/10
  sv_mat[24,5:7]=sv_mat[24,5:7]+rank(-ll_mat[am_id,5:7])/10
  
}

#保存数据
write.csv(sv_mat,"D:\\paper\\2CRMDH\\Result\\Empirical\\likelihood_result_show_II.csv")

#latex代码
library(stargazer)
stargazer(sv_mat, summary=FALSE, rownames=T,digits = 1) 
