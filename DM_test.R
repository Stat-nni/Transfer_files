#计算预测数据与真实数据的各种损失函数，及计算第三节的Mincer-Zarnowitz regression检验



################
#####预定义#####
################
library(forecast)
fore_num=100

name_list=c("000001","000002","000651","000858","600000","600016","600019","600028",
            "600104","600276","600519","600690","AAPL","BAC","GE",
            "IBM","INTC","JPM","KO","MCD","MSFT","WMT")

md_num=6 #自左往右为GARCH,RGARCH,RMDH,2CRMDH,HAR-RV和MMDH的预测结果
DM_mat_Cor=matrix(0,length(name_list),md_num-1)
DM_mat_unCor=matrix(0,length(name_list),md_num-1)

#计算
for (i in 1:length(name_list)){
  
  
  ##################
  #####读取数据#####
  #读取真实数据
  symbol = name_list[i]
  last_name=paste(symbol,"_otc_RV_5min",".csv",sep="")
  family_name=paste("D:\\paper\\2CRMDH\\Peng\\5min_regdata_otc_final\\",last_name,sep="")
  data=read.csv(family_name)[,2:5]
  data_RV = data$V3
  
  #筛选波动率数据
  start_index=which.max(data[,1]-20120900>0)
  days_num=nrow(data)
  Lw=days_num-fore_num-start_index+1
  real_mat=data[(start_index+Lw):(start_index+Lw+fore_num-1),3] #第3列是波动率？认为是的
  
  #读取预测数据
  last_name=paste(symbol,"_forecast_RV1.csv",sep="")
  family_name=paste("D:\\paper\\2CRMDH\\Result\\Empirical\\forecast\\RV\\combine\\",last_name,sep="")
  hat_mat=read.csv(family_name)[,2:(md_num+1)]
  
  ####################
  #####进行DM检验#####
  #计算残差
  res_GARCH=abs(hat_mat[,1]-real_mat)
  res_RGARCH=abs(hat_mat[,2]-real_mat)
  res_RMDH=abs(hat_mat[,3]-real_mat)
  res_2CRMDH=abs(hat_mat[,4]-real_mat)
  res_HARRV=abs(hat_mat[,5]-real_mat)
  res_MMDH=abs(hat_mat[,6]-real_mat)
  
  #进行检验
  h_len=length(res_GARCH)^(1/3)+1
  DM_mat_Cor[i,1]=round(dm.test(res_2CRMDH,res_GARCH,alternative = "less",varestimator="acf",h=h_len,power=1)$p.value,2)
  DM_mat_Cor[i,2]=round(dm.test(res_2CRMDH,res_RGARCH,alternative = "less",varestimator="acf",h=h_len,power=1)$p.value,2)
  DM_mat_Cor[i,3]=round(dm.test(res_2CRMDH,res_RMDH,alternative = "less",varestimator="acf",h=h_len,power=1)$p.value,2)
  DM_mat_Cor[i,4]=round(dm.test(res_2CRMDH,res_HARRV,alternative = "less",varestimator="acf",h=h_len,power=1)$p.value,2)
  DM_mat_Cor[i,5]=round(dm.test(res_2CRMDH,res_MMDH,alternative = "less",varestimator="acf",h=h_len,power=1)$p.value,2)
  
  DM_mat_unCor[i,1]=round(dm.test(res_2CRMDH,res_GARCH,alternative = "less",varestimator="bartlett",h=h_len,power=1)$p.value,2)
  DM_mat_unCor[i,2]=round(dm.test(res_2CRMDH,res_RGARCH,alternative = "less",varestimator="bartlett",h=h_len,power=1)$p.value,2)
  DM_mat_unCor[i,3]=round(dm.test(res_2CRMDH,res_RMDH,alternative = "less",varestimator="bartlett",h=h_len,power=1)$p.value,2)
  DM_mat_unCor[i,4]=round(dm.test(res_2CRMDH,res_HARRV,alternative = "less",varestimator="bartlett",h=h_len,power=1)$p.value,2)
  DM_mat_unCor[i,5]=round(dm.test(res_2CRMDH,res_MMDH,alternative = "less",varestimator="bartlett",h=h_len,power=1)$p.value,2)
  
}

#生成LateX代码
library(stargazer)
stargazer(DM_mat_Cor, summary=FALSE, rownames=T,digits = 2) 
stargazer(DM_mat_unCor, summary=FALSE, rownames=T,digits = 2) 