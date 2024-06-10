#计算预测数据与真实数据的各种损失函数，及计算第三节的Mincer-Zarnowitz regression检验



################
#####预定义#####
################
fore_num=100

name_list=c("000001","000002","000651","000858","600000","600016","600019","600028",
             "600104","600276","600519","600690","AAPL","BAC","GE",
             "IBM","INTC","JPM","KO","MCD","MSFT","WMT")



##################################################################
#####具体计算各loss function及Mincer-Zarnowitz regression检验#####
##################################################################
#定义数组
md_num=7 #md_num为6应当是自左往右为GARCH,RGARCH,RMDH,2CRMDH,HAR-RV和MMDH以及MIDAS模型的预测结果

MAE_mat=matrix(0,(length(name_list)+2),md_num) #+2是为储存各自的平均值
RMSE_mat=matrix(0,(length(name_list)+2),md_num)
HMAE_mat=matrix(0,(length(name_list)+2),md_num)
HRMSE_mat=matrix(0,(length(name_list)+2),md_num)
AMAPE_mat=matrix(0,(length(name_list)+2),md_num)

alpha_mat_RMDH=matrix(0,length(name_list),md_num)
alpha_mat_2CRMDH=matrix(0,length(name_list),md_num)

#计算中国组数据
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
  # hat_mat=read.csv("D:\\paper\\2CRMDH\\000001_forecast_gm.csv")
  
  
  ###########################
  #####计算loss function#####
  if (i <= 12){
    #计算中国组的损失函数
    for (k in 1:md_num){
      MAE_mat[i,k]=mean(abs(hat_mat[,k]-real_mat))
      RMSE_mat[i,k]=sqrt(mean((hat_mat[,k]-real_mat)^2))
      HMAE_mat[i,k]=mean(abs(hat_mat[,k]/real_mat-1))
      HRMSE_mat[i,k]=sqrt(mean((hat_mat[,k]/real_mat-1)^2))
      AMAPE_mat[i,k]=mean(abs((hat_mat[,k]-real_mat)/(hat_mat[,k]+real_mat)))
    }
  }
  else{
    #计算美国组的损失函数
    for (k in 1:md_num){
      MAE_mat[(i+1),k]=mean(abs(hat_mat[,k]-real_mat))
      RMSE_mat[(i+1),k]=sqrt(mean((hat_mat[,k]-real_mat)^2))
      HMAE_mat[(i+1),k]=mean(abs(hat_mat[,k]/real_mat-1))
      HRMSE_mat[(i+1),k]=sqrt(mean((hat_mat[,k]/real_mat-1)^2))
      AMAPE_mat[(i+1),k]=mean(abs((hat_mat[,k]-real_mat)/(hat_mat[,k]+real_mat)))
    }
  }
  
  
  #############################################
  #####进行Mincer-Zarnowitz regression检验#####
  ###针对RMDH模型的拟合结果进行检验###
  fit_RMDH=lm(real_mat~hat_mat[,3])
  
  #计算alpha0&alpha1
  fit_coef=coef(fit_RMDH)
  alpha_mat_RMDH[i,1:2]=round(fit_coef,2)
  
  #计算p-value
  fit_cov=vcov(fit_RMDH)
  fit_s=(fit_coef-c(0, 1)) %*% solve(fit_cov) %*% c(fit_coef-c(0, 1))
  alpha_mat_RMDH[i,3]=round(pchisq(fit_s, 2, lower.tail = FALSE),3)
  
  #计算R2
  alpha_mat_RMDH[i,4]=round(summary(fit_RMDH)$r.squared,2)
  
  #储存alpha系数的p-value
  alpha_mat_RMDH[i,5:6]=round(summary(fit_RMDH)$coef[,4],2)
  
  ###针对2CRMDH模型的拟合结果进行检验###
  fit_2CRMDH=lm(real_mat~hat_mat[,4])
  
  #计算alpha0&alpha1
  fit_coef=coef(fit_2CRMDH)
  alpha_mat_2CRMDH[i,1:2]=round(fit_coef,2)
  
  #计算p-value
  fit_cov=vcov(fit_2CRMDH)
  fit_s=(fit_coef-c(0, 1)) %*% solve(fit_cov) %*% c(fit_coef-c(0, 1))
  alpha_mat_2CRMDH[i,3]=round(pchisq(fit_s, 2, lower.tail = FALSE),3)
  
  #计算R2
  alpha_mat_2CRMDH[i,4]=round(summary(fit_2CRMDH)$r.squared,2)
  
  #储存alpha系数的p-value
  alpha_mat_2CRMDH[i,5:6]=round(summary(fit_2CRMDH)$coef[,4],2)
}


####################
#####计算平均值#####
#计算中国组平均值数据
for (ch_id in 1:12){
  
  MAE_mat[13,]=MAE_mat[13,]+rank(MAE_mat[ch_id,])/12
  RMSE_mat[13,]=RMSE_mat[13,]+rank(RMSE_mat[ch_id,])/12
  HMAE_mat[13,]=HMAE_mat[13,]+rank(HMAE_mat[ch_id,])/12
  HRMSE_mat[13,]=HRMSE_mat[13,]+rank(HRMSE_mat[ch_id,])/12
  AMAPE_mat[13,]=AMAPE_mat[13,]+rank(AMAPE_mat[ch_id,])/12
  
}

#计算美国组平均值数据
for (am_id in 14:23){
  
  MAE_mat[24,]=MAE_mat[24,]+rank(MAE_mat[am_id,])/10
  HMAE_mat[24,]=HMAE_mat[24,]+rank(HMAE_mat[am_id,])/10
  HRMSE_mat[24,]=HRMSE_mat[24,]+rank(HRMSE_mat[am_id,])/10
  AMAPE_mat[24,]=AMAPE_mat[24,]+rank(AMAPE_mat[am_id,])/10
  
  if (am_id==20){
    RMSE_mat[24,]=RMSE_mat[24,]+c(4,1,3,7,5,6,1)/10
  }
  else{
    RMSE_mat[24,]=RMSE_mat[24,]+rank(RMSE_mat[am_id,])/10
  }
  
}

######################
#####储存计算结果#####
######################
#储存列名
colnames(MAE_mat)=c('GARCH','RGARCH','RMDH','2CRMDH','HAR-RV','MMDH','MIDAS')
colnames(RMSE_mat)=c('GARCH','RGARCH','RMDH','2CRMDH','HAR-RV','MMDH','MIDAS')
colnames(HMAE_mat)=c('GARCH','RGARCH','RMDH','2CRMDH','HAR-RV','MMDH','MIDAS')
colnames(HRMSE_mat)=c('GARCH','RGARCH','RMDH','2CRMDH','HAR-RV','MMDH','MIDAS')
colnames(AMAPE_mat)=c('GARCH','RGARCH','RMDH','2CRMDH','HAR-RV','MMDH','MIDAS')

#储存行名
row_name=c('PING AN BANK','CHINA VANKE','GREE','WULIANGYE','SPD BANK','MINSHENG BANK','BAOSTEEL','SINOPEC','SAIC MOTOR','HENGRUI MEDICINE','KWEICHOW MOUTAI','HAIER','Averaged Ranks (Chinese)','AAPL','BAC','GE','IBM','INTC','JPM','KO','MCD','MSFT','WMT','Averaged Ranks (US)')
rownames(MAE_mat)=row_name
rownames(RMSE_mat)=row_name
rownames(HMAE_mat)=row_name
rownames(HRMSE_mat)=row_name
rownames(AMAPE_mat)=row_name

#转换顺序
use_od=c('GARCH','RGARCH','HAR-RV','MMDH','RMDH','MIDAS','2CRMDH')

MAE_mat=MAE_mat[,use_od]
RMSE_mat=RMSE_mat[,use_od]
HMAE_mat=HMAE_mat[,use_od]
HRMSE_mat=HRMSE_mat[,use_od]
AMAPE_mat=AMAPE_mat[,use_od]

#储存
write.csv(MAE_mat,"D:\\paper\\2CRMDH\\Result\\Empirical\\MAE.csv")
write.csv(RMSE_mat,"D:\\paper\\2CRMDH\\Result\\Empirical\\RMSE.csv")
write.csv(HMAE_mat,"D:\\paper\\2CRMDH\\Result\\Empirical\\HMAE.csv")
write.csv(HRMSE_mat,"D:\\paper\\2CRMDH\\Result\\Empirical\\HRMSE.csv")
write.csv(AMAPE_mat,"D:\\paper\\2CRMDH\\Result\\Empirical\\AMAPE.csv")

write.csv(alpha_mat_RMDH,"D:\\paper\\2CRMDH\\Result\\Empirical\\RMDH_alpha.csv")
write.csv(alpha_mat_2CRMDH,"D:\\paper\\2CRMDH\\Result\\Empirical\\2CRMDH_alpha.csv")

#生成LateX代码
library(stargazer)
stargazer(RMSE_mat, summary=FALSE, rownames=T,digits = 4) 
