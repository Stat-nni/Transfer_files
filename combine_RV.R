#合并已有RV与MIDAS的结果

name_list=c("000001","000002","000651","000858","600000","600016","600019","600028",
            "600104","600276","600519","600690","AAPL","BAC","GE",
            "IBM","INTC","JPM","KO","MCD","MSFT","WMT")

for (q in 1:length(name_list)){
  
  #读取已有RV
  symbol=name_list[q]
  last_name=paste(symbol,"_forecast_RV1",".csv",sep="")
  family_name=paste("D:\\paper\\2CRMDH\\Result\\Empirical\\forecast\\RV",last_name,sep="\\")
  data=read.csv(family_name)[,2:7]
  
  #读取MIDAS-RV
  midas_name=paste("D:\\paper\\2CRMDH\\Result\\Empirical\\forecast\\RV\\MIDAS",last_name,sep="\\")
  data_md=read.csv(midas_name)[,2]
  
  #合并
  mer_data=cbind(data,data_md)
  colnames(mer_data)[7]='V7'
  
  #储存数据
  save_name=paste("D:\\paper\\2CRMDH\\Result\\Empirical\\forecast\\RV\\combine",last_name,sep="\\")
  write.csv(mer_data,save_name)
  
}
