#本文件储存了需要使用的似然函数



################################
#####GARCH(1,1)的Likelihood#####
################################
llfg=function(param,data){#未知波动率
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  
  #计算初值
  h[1]=1 #omega/(1-beta-gamma)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  
  #计算似然函数值
  QL=0
  llvalue1=0
  for (i in 2:t){
    
    #迭代计算
    h[i]=omega+beta*h[i-1]+gamma*(z[i-1]*sqrt(h[i-1]))^2
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    
    #计算似然
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2
    QL=QL+llvalue1
  }
  
  #返回结果
  return(-QL)
}

geth_llfg=function(param,data){
  #根据估计的参数值和可观测值返回h全时期的值
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  
  #计算h
  h[1]=1 #omega/(1-beta-gamma)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  for (i in 2:t){
    h[i]=omega+beta*h[i-1]+gamma*(z[i-1]*sqrt(h[i-1]))^2
    z[i]=(data[i,1]-mu)/sqrt(h[i])
  }
  
  #返回结果
  return(h)
}

pllfg=function(mu,data,ht){
  #计算GARCH的似然---需知晓波动率
  
  t=nrow(data)
  h=ht
  z=numeric(t)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  QL=0
  llvalue1=0
  for (i in 2:t){
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(2*pi)
    QL=QL+llvalue1
  }
  return(QL)
}

iterg=function(param,data){
  #根据参数和数据计算下期一期的波动率
  
  n=nrow(data)
  h=numeric(n+1)
  h[1]=1 #param[2]/(1-param[3]-param[4]) #取1是不是太武断了？？？
  for (i in 1:n){
    h[i+1]=param[2]+param[3]*h[i]+param[4]*(data[i,1]-param[1])^2
  }
  x_pred = h[n+1]
  return(x_pred)
}



####################################
#####Realized GARCH的Likelihood#####
####################################
llf0=function(param,data){#未知波动率
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  xi=param[5]
  tau1=param[6]
  tau2=param[7]
  sigmau=param[8]
  phi=param[9]
  
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  tau=numeric(t)
  g=numeric(t)
  
  #计算初值
  h[1]=1 #exp((omega+gamma*phi)/(1-beta-gamma*phi))
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
  g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
  
  #计算似然函数值
  QL=0
  llvalue1=0
  llvalue2=0
  for (i in 2:t){
    
    #迭代计算
    h[i]=exp(omega+beta*log(h[i-1])+gamma*log(data[i-1,2]))
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
    g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
    
    #计算似然
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
    llvalue2=-(g[i])^2/(2*sigmau^2)
    QL=QL+llvalue1+llvalue2
  }
  
  #返回结果
  return(-QL)
}

geth_llf0=function(param,data){
  #根据估计的参数值和可观测值返回h全时期的值
  
  #获取参数
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  phi=param[9]
  t=nrow(data)
  h=numeric(t)
  
  #计算h
  h[1]=1 #exp((omega+gamma*phi)/(1-beta-gamma*phi))
  for (i in 2:t){
    h[i]=exp(omega+beta*log(h[i-1])+gamma*log(data[i-1,2]))
  }
  
  #返回结果
  return(h)
}

pllf0=function(param,data,ht){
  #计算RGARCH的偏似然---需知晓波动率
  
  mu=param[1]
  xi=param[2]
  phi=param[3]
  tau1=param[4]
  tau2=param[5]
  sigmau=param[6]
  t=nrow(data)
  h=ht
  z=numeric(t)
  tau=numeric(t)
  g=numeric(t)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
  g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
  QL=0
  llvalue1=0
  llvalue2=0
  for (i in 2:t){
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
    g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
    llvalue2=-(g[i])^2/(2*sigmau^2)-log(2*pi)
    QL=QL+llvalue1+llvalue2
  }
  return(QL)
}

iter0=function(param,data){
  #根据参数和数据计算下期一期的波动率
  
  n=nrow(data)
  h=numeric(n+1)
  h[1]=1 #exp((param[2]+param[4]*param[9])/(1-param[3]-param[4]*param[9]))
  for (i in 1:n){
    h[i+1]=exp(param[2]+param[3]*log(h[i])+param[4]*log(data[i,2]))
  }
  x_pred = h[n+1]
  return(x_pred)
}



###################################
#####MDH_RealGARCH的Likelihood#####
###################################
llf1=function(param,data){#未知波动率
  
  #获取参数
  mu=param[1]
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  zeta=param[5]
  xi=param[6]
  tau1=param[7]
  tau2=param[8]
  m0=param[9]
  m1=param[10]
  c=param[11]
  sigmau=param[12]
  phi=param[13]
  
  t=nrow(data)
  h=numeric(t)
  z=numeric(t)
  tau=numeric(t)
  g=numeric(t)
  
  #计算初值
  h[1]=1 #(omega+xi*m0+gamma*exp(xi))/(1-beta-xi*m1)
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
  g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
  
  #计算似然函数值
  QL=0
  llvalue1=0
  llvalue2=0
  llvalue3=0
  for (i in 2:t){
    
    #迭代计算
    h[i]=omega+beta*h[i-1]+gamma*data[i-1,2]+zeta*data[i-1,3]
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
    g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
    
    #计算似然
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
    llvalue2=-(g[i])^2/(2*sigmau^2)
    llvalue3=-m0-m1*h[i]+data[i,3]/c*log(m1*h[i]+m0)-lgamma(data[i,3]/c+1)-log(c)
    QL=QL+llvalue1+llvalue2+llvalue3
  }
  
  #返回结果
  return(-QL)
}

geth_llf1=function(param,data){
  #根据估计的参数值和可观测值返回h全时期的值
  
  #获取参数
  omega=param[2]
  beta=param[3]
  gamma=param[4]
  zeta=param[5]
  xi=param[6]
  m0=param[9]
  m1=param[10]
  t=nrow(data)
  h=numeric(t)
  
  #计算h
  h[1]=1 #(omega+xi*m0+gamma*exp(xi))/(1-beta-xi*m1)
  for (i in 2:t){
    h[i]=omega+beta*h[i-1]+gamma*data[i-1,2]+zeta*data[i-1,3]
  }
  
  #返回结果
  return(h)
}

iter1=function(param,data){
  #根据参数和数据计算下期一期的波动率
  
  n=nrow(data)
  h=numeric(n+1)
  h[1]=1 #(param[2]+param[6]*param[9]+param[4]*exp(param[6]))/(1-param[3]-param[6]*param[10])
  for (i in 1:n){
    h[i+1]=param[2]+param[3]*h[i]+param[4]*data[i,2]+param[5]*data[i,3]
  }
  x_pred=h[n+1]
  return(x_pred)
}



##################################
#####2CRMDH-GARCH的likelihood#####
##################################
llf3=function(param,data){
  
  #导入参数
  mu=param[1]
  omega1=param[2]
  omega2=param[3]
  beta11=param[4]
  beta12=param[5]
  beta21=param[6]
  beta22=param[7]
  gamma1=param[8]
  gamma2=param[9]
  zeta1=param[10]
  zeta2=param[11]
  xi=param[12]
  tau1=param[13]
  tau2=param[14]
  m0=param[15]
  m1=param[16]
  m2=param[17]
  c=param[18]
  sigmau=param[19]
  phi=param[20]
  
  #计算初值
  t=nrow(data)
  h1=numeric(t)
  h2=numeric(t)
  h=numeric(t)
  z=numeric(t)
  tau=numeric(t)
  g=numeric(t)
  h1[1]=1 #(omega1+omega2+(gamma1+gamma2)*exp(xi)+(zeta1+zeta2)*c*m0)/(2-(beta11+beta12+beta21+beta22)-(zeta1+zeta2)*c*(m1+m2))
  h2[1]=1
  h[1]=h1[1]+h2[1]
  z[1]=(data[1,1]-mu)/sqrt(h[1])
  tau[1]=tau1*z[1]+tau2*(z[1]^2-1)
  g[1]=log(data[1,2])-xi-tau[1]-phi*log(h[1])
  
  #计算似然函数值
  QL=0
  llvalue1=0
  llvalue2=0
  llvalue3=0
  for (i in 2:t){
    
    #计算辅助函数
    h1[i]=omega1+beta11*h1[i-1]+beta12*h2[i-1]+gamma1*data[i-1,2]+zeta1*data[i-1,3]
    # 所以导入data的第二列是RV，第三列是交易量（去趋势后的）
    
    h2[i]=omega2+beta21*h1[i-1]+beta22*h2[i-1]+gamma2*data[i-1,2]+zeta2*data[i-1,3]
    
    
    h[i]=h1[i]+h2[i]
    z[i]=(data[i,1]-mu)/sqrt(h[i])
    tau[i]=tau1*z[i]+tau2*(z[i]^2-1)
    g[i]=log(data[i,2])-xi-tau[i]-phi*log(h[i])
    
    #计算似然
    llvalue1=-1/2*log(h[i])-1/2*z[i]^2-1/2*log(sigmau^2)
    llvalue2=-(g[i])^2/(2*sigmau^2)
    llvalue3=-m0-m1*h1[i]-m2*h2[i]+(data[i,3]/c)*log(m1*h1[i]+m2*h2[i]+m0)-lgamma((data[i,3]/c)+1)-
      log(c)
    QL=QL+llvalue1+llvalue2+llvalue3
  }
  QL=QL/(t-1)
  return(-QL)
}

geth_llf3=function(param,data){
  #根据估计的参数值和可观测值返回h全时期的值
  
  #获取参数
  omega1=param[2]
  omega2=param[3]
  beta11=param[4]
  beta12=param[5]
  beta21=param[6]
  beta22=param[7]
  gamma1=param[8]
  gamma2=param[9]
  zeta1=param[10]
  zeta2=param[11]
  xi=param[12]
  m0=param[15]
  m1=param[16]
  m2=param[17]
  c=param[18]
  t=nrow(data)
  h1=numeric(t)
  h2=numeric(t)
  h=numeric(t)
  
  #计算h
  h1[1]=1 #(omega1+omega2+(gamma1+gamma2)*exp(xi)+(zeta1+zeta2)*c*m0)/(2-(beta11+beta12+beta21+beta22)-(zeta1+zeta2)*c*(m1+m2))
  h2[1]=h1[1]
  h[1]=2*h1[1]
  for (i in 2:t){
    h1[i]=omega1+beta11*h1[i-1]+beta12*h2[i-1]+gamma1*data[i-1,2]+zeta1*data[i-1,3]
    h2[i]=omega2+beta21*h1[i-1]+beta22*h2[i-1]+gamma2*data[i-1,2]+zeta2*data[i-1,3]
    h[i]=h1[i]+h2[i]
    
  }
  
  #返回结果
  return(h)
}

iter3=function(param,data){
  #根据参数和数据计算下期一期的波动率
  
  n=nrow(data)
  h1=numeric(n+1)
  h2=numeric(n+1)
  h=numeric(n+1)
  h1[1]=1 #(param[2]+param[3]+(param[8]+param[9])*exp(param[12])+(param[10]+param[11])*param[18]*param[15])/(2-(param[4]+param[5]+param[6]+param[7])-(param[10]+param[11])*param[18]*(param[16]+param[17]))
  h2[1]=h1[1]
  h[1]=2*h1[1]
  for (i in 1:n){
    h1[i+1]=param[2]+param[4]*h1[i]+param[5]*h2[i]+param[8]*data[i,2]+param[10]*data[i,3]
    h2[i+1]=param[3]+param[6]*h1[i]+param[7]*h2[i]+param[9]*data[i,2]+param[11]*data[i,3]
    h[i+1]=h1[i+1]+h2[i+1]
  }
  x_pred=h[n+1]
  return(x_pred)
}



#################################
#####GARCH-MIDAS的likelihood#####
#################################
kernel_phi_vec<-function(K,omega1,omega2){
  
  #需要的核函数
  vec = numeric(K)
  sum = 0
  for (i in 1:K){
    sum = sum + ((i/K)**(omega1-1))*((1-i/K)**(omega2-1))
  }
  for (i in 1:K){
    tmp1 = ((i/K)**(omega1-1))*((1-i/K)**(omega2-1))
    tmp2 = tmp1/sum
    vec[i] = tmp2
  }
  return(rev(vec))
}


llfgm<-function(param,data_5min,data_RV,index0,K){
  #GARCH-MIDAS
  #似然函数本尊
  mu<-param[1]
  alpha = param[2]
  beta = param[3]
  mbar = param[4]
  theta = param[5]
  omega1 = param[6]
  omega2 = param[7]
  vec_tmp = kernel_phi_vec(K,param[6],param[7])
  t_day = length(data_RV)
  g = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  g[1,K+1] = 1#初始化这个即可
  m = numeric(t_day)
  m[K+1] = mbar + theta * vec_tmp%*%data_RV[(index0):(index0-1+K)]#从K+1开始算，
  #最后得到(250-K)个数据
  h = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  h[1,K+1] = g[1,K+1]*m[1]
  QL<-0
  llvalue1<-0
  for(t in (K+1):t_day){#这样给了多少天，会有前K天的数据只用来计算别的
    # t=K+1
    data=data_5min[,t]#取出第t天的数据
    for (i in 1:length(data)){#这个都是固定了天数t得到的结果，对观测点i
      # i=1
      m[t]<- (mbar + theta * vec_tmp%*%data_RV[(index0+t-K-1):(index0-2+t)])
      if(i==1){
        if(t==K+1){
          g[i,t]=g[1,K+1]
        }else{
          g[i,t]<-(1-alpha-beta) + alpha * ((data_5min[length(data),(t-1)]-mu)**2)/m[t] + 
            beta * g[length(data),(t-1)]
        }
      }else{
        g[i,t]<- (1-alpha-beta) + alpha * ((data[i-1]-mu)**2)/m[t] + beta * g[i-1,t]
      }
      h[i,t]<- m[t]*g[i,t]#
      z = (data[i]-mu)/sqrt(h[i,t])
      llvalue1<-llvalue1+(-1/2*log(h[i,t])-1/2*z^2)
    }
    QL<-QL+llvalue1
  }
  return(-QL)
}

penalty_function <- function(param,data_5min,data_RV,index0,K) {
  # 如果 x 不满足约束条件，返回一个大的惩罚值
  if ((param[2]+param[3]) >= 1) {
    return(1e6)  # 可以根据具体问题设置适当的惩罚值
  } else {
    return(0)
  }
}

llfgm_with_penalty<-function(param,data_5min,data_RV,index0,K){
  llfgm(param,data_5min,data_RV,index0,K)+ penalty_function(param,data_5min,data_RV,index0,K)
}

itergm<-function(param,data_5min,data_RV,index0,K){
  param<-theta_result
  data_5min<-data_est
  mu<-param[1]
  alpha = param[2]
  beta = param[3]
  mbar = param[4]
  theta = param[5]
  omega1 = param[6]
  omega2 = param[7]
  vec_tmp = kernel_phi_vec(K,param[6],param[7])
  t_day = length(data_RV)
  g = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  g[1,K+1] = 1#初始化这个即可
  m = numeric(t_day+1)
  m[K+1] = mbar + theta * vec_tmp%*%data_RV[(index0):(index0-1+K)]#从K+1开始算，
  #最后得到(250-K)个数据
  h = matrix(0,ncol=t_day,nrow = nrow(data_5min))
  h[1,K+1] = g[1,K+1]*m[1]
  for(t in (K+1):t_day){#这样给了多少天，会有前K天的数据只用来计算别的
    # t=K+1
    data=data_5min[,t]#取出第t天的数据
    for (i in 1:length(data)){#这个都是固定了天数t得到的结果，对观测点i
      # i=2
      m[t]<- (mbar + theta * vec_tmp%*%data_RV[(index0+t-K-1):(index0-2+t)])
      if(i==1){
        if(t==K+1){
          g[i,t]=g[1,K+1]
        }else{
          g[i,t]<-(1-alpha-beta) + alpha * ((data_5min[length(data),(t-1)]-mu)**2)/m[t] + 
            beta * g[length(data),(t-1)]
        }
      }else{
        g[i,t]<- (1-alpha-beta) + alpha * ((data[i-1]-mu)**2)/m[t] + beta * g[i-1,t]
      }
      h[i,t]<- m[t]*g[i,t]#
    }
  }
  #跑完循环后，g矩阵已经全部得到了
  m[t_day+1]<-mbar + theta * vec_tmp%*%data_RV[(index0+t_day-K):(index0-1+t_day)]
  result_day_vol<-m[t_day+1]*(48+(1-((alpha+beta)^48)/(1-(alpha+beta)))*g[48,t_day])
  return(result_day_vol)
}



##############################
#####MMDH模型估计时的函数#####
##############################
g1 = function(param,x)
{
  #传递参数
  rbar = param[1]
  omega = param[2]
  beta = param[3]
  sigmau = param[4]
  c = param[5]
  m0 = param[6]
  m1 =  param[7]
  mu = omega/(1-beta)
  sigmasq = sigmau^2/(1-beta^2)
  
  #---????K?ľ?
  K_half = exp(mu/2+sigmasq/8)
  K_one = exp(mu + sigmasq/2)
  K_oahalf = exp(3*mu/2 + 9*sigmasq/8)
  K_two = exp(2*mu + 2*sigmasq)
  K_three = exp(3*mu + 9*sigmasq/2)
  #---????V?ľ?--
  V_one = c*m0 + c*m1*K_one
  V_two = c*V_one + c^2 * m1^2 * (K_two - K_one^2)
  V_three = c^2 * V_one + 3 * c^3 * m1^2 * (K_two - K_one^2) +
    c^3*m1^3*(K_three - 3 * K_two * K_one + 2 * K_one^3)
  # V_three = c^3 * (m0^3+m1^3*K_three+3*m0*m1^2*K_two+3*m1*m0^2*K_one + 3*(m0^2+2*m0*m1*K_one+m1^2*K_two)+m0+m1*K_one)-
  #   3*V_one^2*c-V_one^3- 3*V_one*c^2*m1^2*(K_two - K_one^2)
  #---????R??V?????ľ?
  RV_one_nc = rbar * V_one
  RV_one = c * sqrt(2/pi) * m1 * (K_oahalf - K_half*K_one)
  # RV_one = c * sqrt(2/pi) * m1 * (K_oahalf - K_half)
  RV_twoone = V_one * K_one + m1 * c * (K_two - K_one^2)
  # RV_twoone = V_one * K_one + m1 * (K_two - K_one^2)
  # RV_two = c * K_one * V_one + c^2 * m1 * (K_two - K_one^2) +
  #   c^2 * m1^2 * (K_three - 3 * K_two * K_one + 2 * K_one^3 - K_one * (K_two - K_one^2))
  RV_two = c^2*((m0^2+m0-2*V_one*m0/c+V_one^2/c^2)*K_one + (2*m0*m1+m1-2*V_one*m1/c)*K_two + m1^2*K_three)
  #---????R???ͺ???
  R_lagone = (2/pi) * K_half^2 *exp(beta*sigmasq/4)
  R_sqlagone = K_one^2 * exp(beta * sigmasq)
  R_lagtwo = (2/pi) * K_half^2 *exp(beta^2 *sigmasq/4)
  R_sqlagtwo = K_one^2 * exp(beta^2 * sigmasq)
  R_lagthree = (2/pi) * K_half^2 *exp(beta^3 *sigmasq/4)
  R_sqlagthree = K_one^2 * exp(beta^3 * sigmasq)
  R_lagfour = (2/pi) * K_half^2 *exp(beta^4 *sigmasq/4)
  R_sqlagfour = K_one^2 * exp(beta^4 * sigmasq)
  R_lagfive = (2/pi) * K_half^2 *exp(beta^5 *sigmasq/4)
  R_sqlagfive = K_one^2 * exp(beta^5 * sigmasq)
  #---????V???ͺ????ľ?
  V_lagone = c^2*m1^2 * (R_sqlagone - K_one^2)
  V_lagtwo = c^2*m1^2 * (R_sqlagtwo - K_one^2)
  V_lagthree = c^2*m1^2 * (R_sqlagthree - K_one^2)
  V_lagfour = c^2*m1^2 * (R_sqlagfour - K_one^2)
  V_lagfive = c^2*m1^2 * (R_sqlagfive - K_one^2)
  #---????R??ƽ????V???ͺ???(V?ͺ?)
  RV_lagone =  c *(m0*K_one + m1 * R_sqlagone)
  RV_lagtwo =  c *(m0*K_one + m1 * R_sqlagtwo)
  RV_lagthree =  c *(m0*K_one + m1 * R_sqlagthree)
  RV_lagfour =  c *(m0*K_one + m1 * R_sqlagfour)
  RV_lagfive =  c *(m0*K_one + m1 * R_sqlagfive)
  #-------??ʼ????---------
  a1 = rbar - x[,1]
  a2 = sqrt(2/pi)*K_half - abs(x[,1]-rbar)
  a3 = K_one - (x[,1]-rbar)^2
  a4 = 2*sqrt(2/pi)*K_oahalf - abs(x[,1]-rbar)^3
  a5 = 3*K_two - (x[,1]-rbar)^4
  b1 = V_one - x[,2]
  b2 = V_two - (x[,2] - V_one)^2
  b3 = V_three - (x[,2] - V_one)^3
  c1 = RV_one_nc - x[,1]*x[,2]
  c2 = RV_one - abs(x[,1]-rbar)*(x[,2]-V_one)
  c3 = RV_twoone - (x[,1]-rbar)^2 * x[,2]
  c4 = RV_two - (x[,1]-rbar)^2 * (x[,2]-V_one)^2
  t = length(x[,1])
  #------------???????ͺ???
  x_lagone = x[2:t,1]
  x_lagone = c(x_lagone,NA)
  d1 = R_lagone - abs(x[,1]-rbar) * abs(x_lagone-rbar)
  d2 = R_sqlagone - (x[,1]-rbar)^2 * (x_lagone-rbar)^2
  x_lagtwo = x[3:t,1]
  x_lagtwo = c(x_lagtwo,rep(NA,2))
  d3 = R_lagtwo - abs(x[,1]-rbar) * abs(x_lagtwo-rbar)
  d4 = R_sqlagtwo - (x[,1]-rbar)^2 * (x_lagtwo-rbar)^2
  x_lagthree = x[4:t,1]
  x_lagthree = c(x_lagthree,rep(NA,3))
  d5 = R_lagthree - abs(x[,1]-rbar) * abs(x_lagthree-rbar)
  d6 = R_sqlagthree - (x[,1]-rbar)^2 * (x_lagthree-rbar)^2
  x_lagfour = x[5:t,1]
  x_lagfour = c(x_lagfour,rep(NA,4))
  d7 = R_lagfour - abs(x[,1]-rbar) * abs(x_lagfour-rbar)
  d8 = R_sqlagfour - (x[,1]-rbar)^2 * (x_lagfour-rbar)^2
  x_lagfive = x[6:t,1]
  x_lagfive = c( x_lagfive,rep(NA,5))
  d9 = R_lagfive - abs(x[,1]-rbar) * abs(x_lagfive-rbar)
  d10 = R_sqlagfive - (x[,1]-rbar)^2 * (x_lagfive-rbar)^2
  #----?ɽ?��?ͺ???
  y_lagone = x[2:t,2]
  y_lagone = c(y_lagone,NA)
  e1 = V_lagone - (x[,2]-V_one)*(y_lagone-V_one)
  y_lagtwo = x[3:t,2]
  y_lagtwo = c(y_lagtwo,rep(NA,2))
  e2 = V_lagtwo - (x[,2]-V_one)*(y_lagtwo-V_one)
  y_lagthree = x[4:t,2]
  y_lagthree = c(y_lagthree,rep(NA,3))
  e3 = V_lagthree - (x[,2]-V_one)*(y_lagthree-V_one)
  y_lagfour = x[5:t,2]
  y_lagfour = c(y_lagfour,rep(NA,4))
  e4 = V_lagfour - (x[,2]-V_one)*(y_lagfour-V_one)
  y_lagfive = x[6:t,2]
  y_lagfive = c(y_lagfive,rep(NA,5))
  e5 = V_lagfive - (x[,2]-V_one)*(y_lagfive-V_one)
  f1 = RV_lagone - y_lagone * (x[,1]-rbar)^2
  f2 = RV_lagtwo - y_lagtwo * (x[,1]-rbar)^2
  f3 = RV_lagthree - y_lagthree * (x[,1]-rbar)^2
  f4 = RV_lagfour - y_lagfour * (x[,1]-rbar)^2
  f5 = RV_lagfive - y_lagfive * (x[,1]-rbar)^2
  result = cbind(a1,a2,a3,a4,a5,b1,b2,b3,c1,c2,c3,c4,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,e1,e2,e3,e4,e5,f1,f2,f3,f4,f5)
  result = na.omit(result)
  return(result)
}