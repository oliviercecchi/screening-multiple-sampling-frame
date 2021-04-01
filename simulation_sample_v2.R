#####################################################################################
# functions
#####################################################################################
# calcul stratified sample variance
se_kirsh_2<-function(sample, population, p){
  # Kish et al., formula 3.3.4, p.82 (pdf p. 96)
  N_h<-population
  N  <- sum(population)
  n_h<- sample
  p_h<- p
  
  W_h=N_h/N  
  f=n_h/N
  
  variance<-sum(
    W_h^2     
    *(1-f)*    
      p_h    *
      (1-p_h) /
      (n_h - 1)
  )
  
  se<- sqrt(variance)
  se
}
# calcul Z
zc <- function(cl){qnorm(1-(1-cl)/2)}

# calcul sample size uf RD
Ssize<-function (x,A=0.95,p=0.5,E=0.05) {ceiling((qchisq(A,1)*x*p*(1-p)) / (E^2*(x-1)+qchisq(A,df=1)*p*(1-p)))}
# x = pop
# A = confidence level
# p = proportion in sample
# E = error margin

calc_sample<-function(sampling_frame,ratio_sample,sample_size,p=0.5,cl=0.95){
  sample<-c(
    ceiling(sample_size*(ratio_sample)),
    ceiling(sample_size*((1-ratio_sample)))
  )
  # population size
  population<-c(
    sampling_frame$yb,
    sampling_frame$yab
  )
  # calculate the standard error and the error margin
  se_kirsh_2(sample, population, p) * zc(cl)
}



#####################################################################################
## set parameters
#####################################################################################
# minimum precision needed
cl<-seq(0.95,0.95,0.01)

# error margin
min_em<-0.1
# proportion in sample
p=.5

# cost 
cost_rdd<-30
cost_list<-15

# sample size rang needed
sample_size<-seq(Ssize(500, min(cl),p, min_em),800,1)

# ratio range between list and RDD 
ratio_list<-seq(0.1,0.9,0.01)


#####################################################################################
## prepare the files
#####################################################################################

library(readxl);library(dplyr)

sampling_frame<-read_excel("Screening multiple frame example.xlsx","Overview_per strata")
# prepare the data
sampling_frame$psu<-paste0(sampling_frame$admin2,"_",sampling_frame$strata)

# calculate the number of hh that are in RDD and not in the list
sampling_frame$yab<-sampling_frame$ya-sampling_frame$yb


#####################################################################################
## implement - quite wasteful since it is making a lot of caculation to assess the best option 

param_list<-expand.grid(cl,sample_size,ratio_list) %>% as.data.frame()
names(param_list)<-c("cl","sample_size","ratio_list")

# calculate error margin if rds
splf<-split.data.frame(sampling_frame,sampling_frame$psu)

result<-lapply(splf,function(ds,param_list,cost_rdd,cost_list,min_em,p){
  em<-lapply(1:nrow(param_list),function(i,param_list,ds,p){
    calc_sample(ds,param_list$ratio_list[i],param_list$sample_size[i],p=p,cl=param_list$cl[i])
  },param_list=param_list,ds=ds,p=p)  
  
  out<-param_list %>% mutate(
    population=ds$ya,
    sample_size_if_rd=Ssize(ds$ya,cl,p,min_em),
    min_em=min_em,
    em=unlist(em),
    name=ds$psu,
    survey_in_list=ceiling(sample_size*(ratio_list)),
    survey_in_rdd=ceiling(sample_size*(1-ratio_list)),
    cost=survey_in_rdd*cost_rdd+ survey_in_list*cost_list
  )
  
  out<-out[which(out$em<=min_em),]
  return(out)
  
},param_list=param_list,cost_list=cost_list,cost_rdd=cost_rdd,min_em=min_em,p=p) %>% bind_rows() %>% as.data.frame



# decision 
result<-result %>% mutate(
  reference_cost=sample_size_if_rd*cost_rdd,
  cheaper=cost<=reference_cost)

result<- result %>% split.data.frame(.,.$name) %>% 
  lapply(.,
        function(x){
          is_cheaper<-which(x$cheaper)
          if(length(is_cheaper)==0){
            x$diff<-x$cost-x$reference_cost
            is_cheaper<-which(x$diff==min(x$diff))
          }
          
          is_cl<-is_cheaper[x$cl[is_cheaper]==max(x$cl[is_cheaper])]
          is_max_sample<-is_cl[x$sample_size[is_cl]==min(x$sample_size[is_cl])]
          is_max_em<-is_max_sample[x$cost[is_max_sample]==min(x$cost[is_max_sample])]
          is_max_ratio_list<-is_max_em[x$ratio_list[is_max_sample]==min(x$ratio_list[is_max_sample])]

          return(x[is_max_ratio_list,])
        }) %>% bind_rows()

write.csv(result, "simulation_sampling.csv")









