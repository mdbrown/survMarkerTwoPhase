#script to profile code to focus on what we shd code into C

library(profr)

load("exampleData.Rdata")


#####################
## cch 
#estimate accuracy measures using non-parametric estimates
#by setting estimation.method = "NP"
p.cch.np <- profr(survMTP.cch(time =survTime, 
            event = status, 
            marker = Y,
            weights = weights,
            subcoh = subcohort,
            data = cohortData_cch, 
            estimation.method = "NP",  
            predict.time = 2, 
            marker.cutpoint = 0))
plot(p.cch.np)
#estimate accuracy measures using semi-parametric estimates
survMTP.cch(time =survTime, 
            event = status, 
            marker = Y,
            weights = weights,
            subcoh = subcohort,
            data = cohortData_cch, 
            estimation.method = "SP",  
            predict.time = 2, 
            marker.cutpoint = 0)

###################
## ncc
#Nonparametric estimates 


p.ncc.np <- profr(survMTP.ncc(time = survTime, 
            event = status, 
            marker=Y, 
            subcoh=subcohort, 
            id = id, 
            data = cohortData_ncc,
            sets = RiskSets, 
            estimation.method = "NP", 
            predict.time = 2 , 
            marker.cutpoint = 0))

plot(p.ncc.np)

#Semiparametric estimates 
p.ncc.sp <- profr(survMTP.ncc(time = survTime, 
            event = status, 
            marker=Y, 
            subcoh=subcohort, 
            id = id, 
            data = cohortData_ncc,
            sets = RiskSets, 
            estimation.method = "SP", 
            predict.time = 2 , 
            marker.cutpoint = 0))


plot(p.ncc.sp)




### test out some of the C code
load("inputs.Rdata")

CondSurv_FUN_C(IPW, xk, dk, yk, t0,hhat.v )
  
