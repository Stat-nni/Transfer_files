
# 数据处理
library(gmm)
library(dplyr)

Code_path<-'D:\\Academic\\Projects\\2CRMDH\\Code\\Base_code'
Data_path<-'D:\\Academic\\Projects\\2CRMDH\\Data\\Raw_Data'
New_Data_path<-'D:\\Academic\\Projects\\2CRMDH\\Data\\New_Data'

setwd(Data_path)
data_namevec<-list.files()
data_namevec<-data_namevec[substr(data_namevec, 
                                  nchar(data_namevec) - 3, nchar(data_namevec)) 
                           == ".csv"]


for(i in 1:length(data_namevec)){
  setwd(Data_path)
  data<-read.csv(data_namevec[i])
  if(ncol(data)>=5){
    data<-data[,2:5]
  }else{
    data<-data
  }
  colnames(data)<-c('date','time','price','volume')
  data_corrected<-matrix(NA,nrow=length(unique(data$time)),ncol=length(unique(data$date)))
  volume_data_corrected<-matrix(NA,nrow=length(unique(data$time)),ncol=length(unique(data$date)))
  
  for(time_iter in 1:length(unique(data$time))){# 价格数据
    data_aux<-filter(data,time==unique(data$time)[time_iter])
    data_filled<-data_aux$price
    data_corrected[time_iter,]<-data_filled
  }
  for(time_iter in 1:length(unique(data$time))){# 交易量数据
    volume_data_aux<-filter(data,time==unique(data$time)[time_iter])
    volume_data_filled<-volume_data_aux$volume
    volume_data_corrected[time_iter,]<-volume_data_filled
  }
      ## 去除缺失值
      colnames(data_corrected)<-unique(data$date)
      colnames(volume_data_corrected)<-unique(data$date)
      data_corrected<-data_corrected[,apply(data_corrected,2,function(x) !is.na(sum(x)))]
      volume_data_corrected<-
        volume_data_corrected[,apply(volume_data_corrected,2,function(x) !is.na(sum(x)))]
      date_intersection<-intersect(colnames(data_corrected),colnames(volume_data_corrected))
      volume_data_corrected<-volume_data_corrected[,date_intersection]
      data_corrected<-data_corrected[,date_intersection]
      data_corrected_new<-data_corrected
  
  # 去除隔夜效应
      
  # diffs <- sapply(1:(ncol(data_corrected) - 1), function(date) {
  #   diff <- data_corrected[nrow(data_corrected), date] - data_corrected[1, (date + 1)]
  #   return(diff)
  # })|>as.numeric()|>cumsum()
  # 
  # data_corrected_new<-apply(data_corrected,1,function(x){
  #   x[-1]<-x[-1]+diffs
  #   return(x)
  # })|>t() 
  
  # 构建新的数据序列
  data_new<-matrix(NA,ncol=4,nrow=ncol(data_corrected_new))
  colnames(data_new)<-c('date','return','RV','volume')
  
  # 导入日期序列
  data_new[,'date']<-colnames(data_corrected_new)
  
  # 导入收益率
    ## 第一天是open to close
    data_new[1,'return']<-log(data_corrected_new)[nrow(data_corrected_new),1]-log(data_corrected_new)[1,1]
    ## 之后都是close to close
    data_new[-1,'return']<-diff(log(data_corrected_new[nrow(data_corrected_new),]))
    
  # 导入RV
    ## RV只能是open to close
  return_aux<-apply(log(data_corrected_new),2,diff)
  data_new[,'RV']<-apply(return_aux^2,2,sum)
  
  # 导入交易量
      ## 构造交易量
      volume_data<-apply(volume_data_corrected,2,sum)
      plot(apply(volume_data_corrected,1,mean),type='l',main=data_namevec[i])
      
      ## 进行单边去除趋势处理吧，用之前240天的均值作为去趋势
      volume_data_new<-rep(NA,length(volume_data))
      for(i2 in 1:length(volume_data)){
        if(i2<=240){
          volume_data_new[i2]<-volume_data[i2]/mean(volume_data[1:240])
        }else{
          volume_data_new[i2]<-volume_data[i2]/mean(volume_data[(i2-240):(i2-1)])
        }
      }
      
      ## 导入交易量
      data_new[,'volume']<-volume_data_new
      ## 保留新的数据
      setwd(New_Data_path)
      write.csv(data_new,file=data_namevec[i])
}
