#导入函数文件
Code_path<-'D:\\Academic\\Projects\\Corr_CLT'
Data_A_path<-'D:\\Academic\\Projects\\Corr_CLT\\price_I'
Data_B_path<-'D:\\Academic\\Projects\\Corr_CLT\\price_II'

setwd(Data_A_path)
file_A_names<-list.files()# A资产的一千条轨道，每条轨道100天*48个数据
setwd(Data_B_path)
file_B_names<-list.files()# B资产的一千条轨道，每条轨道100天*48个数据

# 参数设置
l=5 # 小的估计窗宽
n=48
Delta=1/n # Delta
kappa= seq(0,1,length.out = 48)
N=100
   
# 读取数据
path_iternum=1
data_A<-read.table(paste(Data_A_path,file_A_names[path_iternum],sep='\\'))
data_A<-matrix(diff(log(data_A[,1])),ncol=100) # 得到收益率数据

data_B<-read.table(paste(Data_B_path,file_B_names[path_iternum],sep='\\'))
data_B<-matrix(diff(log(data_B[,1])),ncol=100) # 得到收益率数据


# 曲线估计
B_i_kappa<-matrix(NA,nrow=length(kappa),ncol=N)
sigma_A_i_kappa<-matrix(NA,nrow=length(kappa),ncol=N)
sigma_B_i_kappa<-matrix(NA,nrow=length(kappa),ncol=N)

for(i in 1:N){
  if(i == 1){
    for(kappa_iternum in 1:length(kappa)){
      j_kappa<-floor(kappa[kappa_iternum]*n)
      j<-j_kappa-l+1
      
      ## local spot covariation
      if(j>=1){
        B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]*data_B[(j:j_kappa),i])
      }else{
        B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
          sum(data_A[(1:l),i]*data_B[(1:l),i])
        )
      }
      
      ## volatility estimators 
      if(j>=1){
        sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]^2)
        sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_B[(j:j_kappa),i]^2)
      }else{
        sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
          sum(data_A[(1:l),i]^2)
        )
        sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
          sum(data_B[(1:l),i]^2)
        )
      }
    }
  }
  else{
    # 对于i>=2的天数，可用前一天的数据进行估计
    for(kappa_iternum in 1:length(kappa)){
      
      j_kappa<-floor(kappa[kappa_iternum]*n)
      j<-j_kappa-l+1
      
      ## local spot covariation
      if(j>=1){
        B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]*data_B[(j:j_kappa),i])
      }else{
        B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
          sum(data_A[(1:j_kappa),i]*data_B[(1:j_kappa),i])+ # 此为当天
            sum(data_A[((n+j):n),i-1]*data_B[((n+j):n),i-1]) # 此为前一天
        )
      }
      
      ## volatility estimators 
      if(j>=1){
        sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_A[(j:j_kappa),i]^2)
        sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*sum(data_B[(j:j_kappa),i]^2)
      }else{
        sigma_A_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
          sum(data_A[(1:j_kappa),i]^2)+ # 此为当天
            sum(data_A[((n+j):n),i-1]^2) # 此为前一天
        )
        sigma_B_i_kappa[kappa_iternum,i]<-1/(l*Delta)*(
          sum(data_B[(1:j_kappa),i]^2)+ # 此为当天
            sum(data_B[((n+j):n),i-1]^2) # 此为前一天
        )
      }
    }
  }
}

# correlation diurnal pattern function estimate
g_rho_kappa<-rep(NA,length(kappa))
g_rho_kappa<-apply(B_i_kappa,1,mean)/(
  (apply(sigma_A_i_kappa,1,mean)^1/2)*(apply(sigma_B_i_kappa,1,mean)^1/2)
                                      )

yita<-sqrt(Delta*sum(g_rho_kappa^2))

f_rho_kappa<-g_rho_kappa/yita


