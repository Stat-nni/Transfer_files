#进行GARCH—MIDAS模型估计



#---------------------Define upper bound and lower bound----------




#-------Recalculate the data-------------
#-------Put the weighted RV to the fourth column--------------


name_list<-c("000001","000002","000651","000858","600000","600016","600019","600028",
             "600104","600276","600519","600690","AAPL","BAC","GE",
             "IBM","INTC","JPM","KO","MCD","MSFT","WMT")



m<-100 # Sample size of out of sample 
delta<-1 # Every delta day a new parameter is estimated


for (q in 1:length(name_list)){
  symbol <- name_list[q]
  IP_2<-paste(symbol,"_otc_RV_5min",".csv",sep="")
  IP_4 <-paste("D:\\paper\\2CRMDH\\Peng\\5min_regdata_otc_final",IP_2,sep="\\")
  data<-read.csv(IP_4)[,2:5]
  data_RV = data$V3
  t<-length(data[,1])
  
  
  index0<-which.max(data[,1]-20120900>0)
  # index1<-which.max(data[,1]-20030700>0)
  win<-t-m-index0+1
  
  
  
  
  forecast_RV1<-matrix(0,m,1)
  param_RV1<-matrix(0,(m/delta)*1,7)
  
  i=1
  
  for (l in 1:(m/delta)){
    data1<-data[(index0+(l-1)*delta):(index0+(l-1)*delta+win-1),2:4]
    t<-length(data1[,1])
    tmp1 = nlminb(inigm,llfgm,data= data1,lower=lowergm,upper= uppergm,K=K,data_RV=data_RV,index0=index0)   
    theta<-tmp1[["par"]]
    param_RV1[l+(m/delta)*(i-1),1:7]<-theta
    #---ÿһ????????ʱ??????K0-------------
    for (k in 1:delta){
      data2<-data[(index0+(l-1)*delta+k-1):(index0+(l-1)*delta+win+k-2),2:4]
      forecast_RV1[delta*(l-1)+k,i]<-iter2(theta,data2,data_RV,index0,K)
    }
    print(c(i,l))
  }
  
  IP_10<-paste("D:\\thesis_project\\project2\\thesis py\\real\\",symbol,"_forecast_gm.csv",sep="")
  write.csv(forecast_RV1,IP_10,row.names=F)
  IP_12<-paste("D:\\thesis_project\\project2\\thesis py\\real\\",symbol,"_param_gm.csv",sep="")
  write.csv(param_RV1,IP_12,row.names=F)
  print(symbol)
}