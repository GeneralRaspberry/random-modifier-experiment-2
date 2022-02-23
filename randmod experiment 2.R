
library("parallel")
library("spatstat")
library("dplyr")
library("ggplot2")
library("ggthemes")
library("RColorBrewer")
library("beepr")
library("reshape2")
library("purrr")
library("writexl")

## tau-leap Gillespie algorithm function
tauLeapG <- function(beta, # transmission rate
                     theta, # dispersal scale
                     b=1, # kernel shape parameter, 1 for exponential
                     sigma=0, # asymptomatic period, used for outputing the time series
                     q0=0, # starting incidence if ppp is without marks
                     q.end=1, # stoping condition 1: incidence lvl
                     t.end=Inf, # stoping condition 2: time after first simulated time step
                     area.host=1, # surface area occupied by one host
                     delta.t=10, # time step
                     ppp, # point pattern as a ppp object, optinally with marks 1/0 for infected/healthy
                     dist.mat=NULL){ # matrix distance if its computation is to be avoided here (for e.g. repeated calls)
  
  ## if the point pattern has no marks, generate some randomly that fits q0
  if (is.null(marks(ppp))){
    inf.start <- max(1, round(ppp$n * q0))
    marks(ppp) <- sample(c(rep(FALSE, ppp$n-inf.start), rep(TRUE, inf.start)))
  }
  
  ## compute distance matrix if not provided
  if (is.null(dist.mat)){ 
    ## add the kernel computation that can be added directly on the dist matrix to reduce comp time
    dist.mat <- exp(-pairdist(ppp)^b / theta^b)
    diag(dist.mat) <- NA
  }
  
  ## function that compute infection event probability, based on the dispersal kernel
  k.norm <- beta * area.host * (b/(2*pi*theta^2*gamma(2/b))) # constant part of the exponential power kernel
  infection <- function(infected, dist){
    inf <-  matrix(k.norm * dist[infected,!infected],
                   nrow=sum(infected), byrow=FALSE)
    inf[is.na(inf)] <- 0
    inf
  }
  
  ## starting time
  time <- 0
  ## inititate the heavy dataframe that will aggregate all changes
  df.big <- data.frame(time=0, who=which(ppp$marks), t(ppp$marks))
  
  ## computation loop
  while (any(!ppp$marks) & time <= t.end & mean(ppp$marks) < q.end){
    ## infection event probaility
    events <- infection(ppp$marks, dist=dist.mat)
    ## random proisson draws
    new.infected <- which(!ppp$marks)[rpois(n=sum(!ppp$marks), lambda=apply(events, 2, sum) * delta.t) > 0]
    ## change marks of newly infected
    ppp$marks[new.infected] <- TRUE
    ## increment time
    time <- time + delta.t
    ## if some infection, increment the big dataframe
    if (length(new.infected) > 0){
      df.big <- rbind(df.big, data.frame(time=time, who=new.infected, t(ppp$marks)))
    }
    ## print a dot per new infection
    # cat(paste0(rep('.', length(new.infected)), collapse = '')) ## comment for quiet
  }
  
  ## make compact, time only, version of the big dataframe
  times.i <- unique(df.big[,1])
  times.d <- times.i + sigma
  times <- sort(unique(c(times.i, times.d)))
  infected <- sapply(times, FUN=function(t) sum(t >= df.big[,1]))
  detectable <- sapply(times, FUN=function(t) sum(t >= df.big[,1] + sigma))
  df.small <- data.frame(time=times, infected=infected, detectable=detectable)
  
  ## out put the simplified time series, and the big one
  list(df.small[df.small$time <= max(df.big$time),], df.big) 
} 




## meta parameters
delta.t <- 100 # time step (ALEX-THIS IS BIGGER THAN THE EXPERIMENT BELOW BECAUSE IT IS TAKING SO MUCH LONGER!)

## epidemic parameters

betavalues <- seq(from=50, to=50,by=50)
thetavalues<-seq (from=40,to=40,by=4)##The data I sent you, which is called data in R is the 1000 realisations of these parameters
theta <- 80
randmodvalues<-seq(from=0, to=1, by=.1)
b <- 1
area.host<-1
infbegin<-1
iter<-1000

##################################add a timer##############################################################

ts<-proc.time()

###########################################################################################################
##Concatenating a list of metric values
##-----------------------------------------

dim<-1000
hosts<-900


tempbind<-c()

sim_par <- function(i=NULL){
  for (j in betavalues){
    for(k in randmodvalues){
      #for (l in 1:(length(betavalues)*iter)){
      # print(l)
      #l<-1  
      
      set.seed(seed=NULL)
      
      radiusCluster<-50
      lambdaParent<-.05
      lambdaDaughter<-25
      randmod<-k
      
      numbparents<-rpois(1,lambdaParent*dim)
      
      xxParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
      yyParent<-runif(numbparents,0+radiusCluster,dim-radiusCluster)
      
      numbdaughter<-rpois(numbparents,(lambdaDaughter))
      sumdaughter<-sum(numbdaughter)
      
      
      
      thetaLandscape<-2*pi*runif(sumdaughter)
      
      rho<-radiusCluster*sqrt(runif(sumdaughter))
      
      
      
      xx0=rho*cos(thetaLandscape)
      yy0=rho*sin(thetaLandscape)
      
      
      xx<-rep(xxParent,numbdaughter)
      yy<-rep(yyParent,numbdaughter)
      
      xx<-xx+xx0
      
      yy<-yy+yy0
      cds<-data.frame(xx,yy)
      is_outlier<-function(x){
        x > dim| x < 0
      }
      cds<-cds[!(is_outlier(cds$xx)|is_outlier(cds$yy)),]
      while (nrow(cds)<hosts){
        dif<-hosts-nrow(cds)
        extraparentxx<-sample(xxParent,dif,replace = TRUE)
        extraparentyy<-sample(yyParent,dif,replace = TRUE)
        extrathetaLandscape<-2*pi*runif(dif)
        extrarho<-radiusCluster*sqrt(runif(dif))
        newextracoodsxx<-extrarho*cos(extrathetaLandscape)
        newextracoodsyy<-extrarho*sin(extrathetaLandscape)
        extraxx<-extraparentxx+newextracoodsxx
        extrayy<-extraparentyy+newextracoodsyy
        cdsextra<-data.frame(xx=extraxx,yy=extrayy)
        cds<-rbind(cds,cdsextra)
      }
      
      sampleselect<-sample(1:nrow(cds),hosts,replace=F)
      cds<-cds%>%slice(sampleselect)
      
      randfunction<-function(x){
        x<-runif(length(x),0,dim)
      }
      randselect<-sample(1:nrow(cds),floor(hosts*randmod),replace=F)
      cds[randselect,]<-apply(cds[randselect,],1,randfunction)
      
      landscape<-ppp(x=cds$xx,y=cds$yy,window=owin(xrange=c(0,dim),yrange=c(0,dim)))
      
      print(length(landscape$marks))
      
      data <- data.frame(x=landscape$x, y=landscape$y, id=1:hosts)
      
      ## design a function that will be called
      
      
      set.seed(seed=NULL)
      marks(landscape2)<- sample(c(rep(TRUE,infbegin), rep(FALSE, hosts-infbegin)))
      
      
      output <- tauLeapG(beta = j, theta = theta, b = b,
                         sigma = 0, delta.t = delta.t,
                         ppp = landscape)
      
      temp <- output[[2]][,1:2][order(output[[2]][,2]),]
      temp<-cbind(temp,beta=j,theta=theta, randmod=k)#,sim=l)
      tempbind<-rbind(tempbind,temp)
      #l<-l+1
    }
  }
  
  datatest1<-data.frame(time=tempbind$time, who=tempbind$who, x=landscape$x[tempbind$who], y=landscape$y[tempbind$who],beta=tempbind$beta,theta=tempbind$theta,randmod=tempbind$randmod)#,sim=tempbind$sim)
  
}


#}



## create a cluster with the set number of cores, say nmax-1
cl <- makeCluster(mc <- getOption("cl.cores", 3))
## call the library loading function in them
clusterCall(cl, function() library("spatstat"))
clusterCall(cl,function() library("ggplot2"))
clusterCall(cl,function() library("tidyverse"))
## export all to the nodes, that's dirty, so run this with a clean environement otherwise your memory will be flooded
clusterExport(cl=cl, varlist=ls())
## call the function in a parallel lapply
par_results <- parLapply(1:iter, fun=sim_par, cl=cl) ## test with 10 first, but then replace 10 by 1000
clusterEvalQ(cl,sim_par)
#simtest<-clusterSplit(cl,seq=1:(iter*length(betavalues)))
#par_results<-cbind(par_results,simtest)
## stop the cluster
stopCluster(cl)
## call cbind on your list of lines to find the matrix you expect
data <- do.call("rbind", par_results)

simtest2<-rep((1:(iter*length(betavalues)*length(randmodvalues))),hosts)
simtest2<-simtest2[order(simtest2)]

data<-cbind(data,sim=simtest2)