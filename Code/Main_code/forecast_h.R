#本代码预测out-of-sample的波动率，为接下来计算各种error做准备



##################
#####准备部分#####
##################

#导入函数文件
Code_path<-'D:\\Academic\\Projects\\2CRMDH\\Code\\Base_code'
Data_path<-'D:\\Academic\\Projects\\2CRMDH\\Data\\Past_Data'
library(gmm)
library(dplyr)
setwd(Code_path)
source('Likelihood_functions.R')

###估计初始值及上下界###
#GARCH
inig=c(0,0.1,0.1,0.2)
lowerg=c(-Inf,0,0,0)
upperg=c(Inf,Inf,Inf,Inf)

#RGARCH
ini0=c(0,0.01,0.1,0.2,0.000,0.000,0.00,0.1,1)
lower0=c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0)
upper0=c(Inf, Inf ,Inf ,Inf ,Inf ,Inf ,Inf ,Inf,Inf)

#RMDH
ini1_bound=c(0,1e-2,0.1,0.2,1e-3,0,0,0,1,1,0.5,0.1,0.8)
lower1_bound=c(-Inf,0,  0,  0,  0,  -Inf,-Inf,-Inf,0,  0,  0,  0,  0)
upper1_bound=c(Inf, Inf,Inf,Inf,Inf,Inf, Inf, Inf, Inf,Inf,Inf,Inf,Inf)

#2CRMDH
ini3= c(0,1e-2,1e-2,0.05,0.05,0.05,0.05,0.1,0.1,1e-3,1e-3,0,0,0,1,1,1,0.5,0.1,1)
lower3=c(-Inf,0,0,0,0,0,0,0,0,0,0,-Inf,-Inf,-Inf,0,0,0,0,0,0)
upper3=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)

#GARCH-MIDAS
K=20
inigm<- c(0,0.1,0.8,0.01,0.1,0,1.5)
lowergm<-c(0,0,0,0,0,0,0)
uppergm<-c(Inf,Inf,Inf,Inf,Inf,300,300)

#合并初始参数
ini_list=list(inig,ini0,ini1_bound,ini3)
lower_list=list(lowerg,lower0,lower1_bound,lower3)
upper_list=list(upperg,upper0,upper1_bound,upper3)

#合并函数名
llf_list=list(llfg,llf0,llf1,llf3) #llf3是2CRMDH的似然函数
iter_list=list(iterg,iter0,iter1,iter3) # GARCH11,RGARCH,MDHRGARCH,2CRMDH

#文件名
name_list=c("000001","000002","000651","000858","600000","600016","600019","600028",
             "600104","600276","600519","600690","AAPL","BAC","GE",
             "IBM","INTC","JPM","KO","MCD","MSFT","WMT")

#定义必要参数
s_list=c(4,9,13,20) #几种GARCH模型的参数长度
fore_num=100 #需要预测波动率的时间长度 
Lw_vec=numeric(length(name_list))



################
######估计######
################

time1=Sys.time()

for (q in 1:length(name_list)){
  
  #读取所选中股票数据
  symbol=name_list[q]
  last_name=paste(symbol,"_otc_RV_5min",".csv",sep="")
  family_name=paste(Data_path,last_name,sep="\\")
  data=read.csv(family_name)[,2:5]
  
  #筛选数据
  days_num=nrow(data)
  start_index=which.max(data[,1]-20120900>0) #依论文所说，筛选9月份之后的数据
  Lw_vec[q]=days_num-fore_num-start_index+1
  Lw=days_num-fore_num-start_index+1
  
  #预测结果储存所需
  forecast_GARCH=matrix(0,fore_num,length(iter_list)) #行为预测的每一天，列为4种GARCH模型
  forecast_HAR=numeric(fore_num) # 预测HAR模型
  forecast_MMDH=numeric(fore_num) # 预测MMDH模型
  
  #为MMDH模型估计获取初始值
  data_tmp=select(data[start_index:(start_index+Lw-1),2:4],-'V3')
  fit_tmp=gmm(g=g1,x=data_tmp,t0=c(rbar=0,omega=0,beta=0.7,sigmau=0.4,c=0.2,m0=1,m1=1),
              type='iterative',vcov='HAC',optfct="nlminb", lower=c(-Inf,-Inf,0.1,1e-1,1e-4,0,0),
             upper=c(Inf, Inf, 1, Inf,Inf,Inf,Inf))
  theta_tmp=fit_tmp[["coefficients"]]
  mu_tmp=theta_tmp[2]/(1-theta_tmp[3])
  sigmasq=theta_tmp[4]^2/(1-theta_tmp[3]^2)
  K0=exp(mu_tmp + sigmasq/2)
  for (j in 1:(Lw-1)){
    K0=K0^theta_tmp[3]*exp(theta_tmp[2]+theta_tmp[4]^2/2) 
  }
  
  #针对每一被预测天来看
  for (day_id in 1:fore_num){
    
    #######################################################
    #####针对GARCH，RGARCH，RMDH，2CRMDH模型计算估计值#####
    
    #获取数据
    data_est=data[(start_index+day_id-1):(start_index+day_id-1+Lw-1),2:4] #用以估计参数的数据
    
    #针对模型来看
    for (md_id in 1:length(iter_list)){
      
      # time_a=Sys.time()
      #获取估计参数
      fit_RV=nlminb(ini_list[[md_id]], 
                    llf_list[[md_id]], 
                    data=data_est, 
                    lower=lower_list[[md_id]], 
                    upper=upper_list[[md_id]])
      theta=fit_RV[["par"]]
      
      # time_b=Sys.time()
      # time_b-time_a
      # #拟合发现：GARCH拟合耗时4.537792 secs
      # #----------RGARCH拟合耗时1.048534 mins
      # #----------RMDH拟合耗时2.371101 mins
      # #----------2CRMDH拟合耗时4.741553 mins
      
      #获取模型预测值
      forecast_GARCH[day_id,md_id]=iter_list[[md_id]](theta,data_est)
      
    }
    
    
    ##################################
    #####针对HAR-RV模型计算估计值#####
    
    #准备填充数据
    data_HAR=matrix(0,Lw,4)
    
    for (k in 1:Lw){
      
      #这是在选取什么数据
      data_HAR[k,1]=data[(start_index+k+day_id-1-1),3]
      data_HAR[k,2]=data[(start_index+k+day_id-1-2),3]
      data_HAR[k,3]=mean(data[(start_index+k+day_id-1-6):(start_index+k+day_id-1-3),3])
      data_HAR[k,4]=mean(data[(start_index+k+day_id-1-23):(start_index+k+day_id-1-7),3])
      
    }
    
    #获取预测结果
    coef=lm(data_HAR[,1]~data_HAR[,2:4])[["coefficients"]]
    forecast_HAR[day_id]=coef[1]+coef[2]*data_HAR[Lw,2]+coef[3]*data_HAR[Lw,3]+coef[4]*data_HAR[Lw,4]
    
    
    ######################################
    #####针对MMDH模型(SARV)计算估计值#####
    fit_result=gmm(g=g1,x=select(data_est,-'V3'),t0=c(rbar=0,omega=0,beta=0.7,sigmau=0.4,c=0.2,m0=1,m1=1),type='iterative',vcov='HAC',optfct="nlminb", lower=c(-Inf,-Inf,0.1,1e-1,1e-4,0,0),
               upper=c(Inf, Inf, 1, Inf,Inf,Inf,Inf))    
    theta=fit_result[["coefficients"]]
    K0=K0^theta[3]*exp(theta[2]+theta[4]^2/2)
    forecast_MMDH[day_id]= K0
      
    
    ######################
    #####输出循环信息#####
    print(day_id)
    
  }
  
  #合并预测结果
  forecast_RV1=cbind(forecast_GARCH,forecast_HAR[1:fore_num],forecast_MMDH[1:fore_num])
  
  #储存估计结果
  IP_10=paste("D:\\paper\\2CRMDH\\Result\\Empirical\\forecast\\RV\\new_250\\",symbol,"_forecast_RV1.csv",sep="")
  write.csv(forecast_RV1,IP_10)
  
  #输出循环信息
  print(symbol)
  
}
time2=Sys.time()

time2-time1



################################
######补充---关于MIDAS模型######
################################

time1=Sys.time()
for (q in 1:length(name_list)){
  
  #读取所选中股票数据
  symbol=name_list[q]
  
  #读取RV数据
  last_name=paste(symbol,"_otc_RV_5min",".csv",sep="")
  family_name=paste("D:\\paper\\2CRMDH\\Peng\\5min_regdata_otc_final",last_name,sep="\\")
  data=read.csv(family_name)[,2:5]
  data_RV = data$V3
  
  #筛选数据
  days_num=nrow(data)
  start_index=which.max(data[,1]-20120900>0) #依论文所说，筛选9月份之后的数据
  Lw=days_num-fore_num-start_index+1
  
  #预测结果储存所需
  forecast_MIDAS=numeric(fore_num)
  
  #针对每一被预测天来看
  for (day_id in 1:fore_num){
    
    #获取数据
    data_est=data[(start_index+day_id-1):(start_index+day_id-1+Lw-1),2] #用以估计参数的数据
    
    #针对模型来看
    fit_RV=nlminb(inigm, llfgm_with_penalty, data_5min=data_est, lower=lowergm, upper=uppergm,K=K,data_RV=data_RV,index0=start_index)
    theta=fit_RV[["par"]]
    
    #获取模型预测值
    forecast_MIDAS[day_id]=itergm(theta,data_est,data_RV,start_index,K)
    
    
    ######################
    #####输出循环信息#####
    print(day_id)
    
  }
  
  #储存估计结果
  IP_10=paste("D:\\paper\\2CRMDH\\Result\\Empirical\\forecast\\RV\\new_250\\MIDAS\\",symbol,"_forecast_RV1.csv",sep="")
  write.csv(forecast_RV1,IP_10)
  
  #输出循环信息
  print(symbol)
  
}
time2=Sys.time()

time2-time1

