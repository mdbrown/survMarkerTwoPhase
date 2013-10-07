## Script to check the standard error estimates under CCH designs


#need the function to simulate data
source("../survAccuracyMeasuresPower/subroutines.R")

#get betas needed for AUC of .5, .6, .7, .8, .9

round(get.Betas(parameter = "AUC", a = .1, predict.time =2, cutoff = 0, 
                parval.H0 =.7, parval.Ha=.9,  f_x = dnorm ), 4)
  
AUCVec  = c(.5, .6, .7, .8, .9)
betaVec = c(0, .3247, 0.6759, 1.1220, 1.9387)

    
ddply(type1.estimates, "AUC", function(df)c("MeanAUCest" = mean(df$AUCest), 
                                            "MeanAUCSE" = mean(df$SE), 
                                            "EmpAUCSE"  = sd(df$AUCest), 
                                            "beta" = mean(df$beta),  
                                            "MeanbetaEST" = mean(df$betaest), 
                                            "MeanbetaSE" = mean(df$betaSE), 
                                            "EmpbetaSE" = sd(df$betaest)))

AUC MeanAUCest  MeanAUCSE   EmpAUCSE   beta  MeanbetaEST MeanbetaSE EmpbetaSE
1 0.5  0.4985593 0.06456151 0.04056454 0.0000 -0.001847552  0.2081279 0.1367906
2 0.6  0.6012137 0.06141551 0.03827648 0.3247  0.346383106  0.2055433 0.1348679
3 0.7  0.6943891 0.05880387 0.03520918 0.6759  0.691708183  0.2182586 0.1422142
4 0.8  0.7912670 0.04991397 0.02526562 1.1220  1.147318749  0.2532104 0.1389707
5 0.9  0.8870104 0.03300787 0.01355715 1.9387  1.964855911  0.3425785 0.1536204

ddply(type2.estimates, "AUC", function(df)c("MeanAUCest" = mean(df$AUCest), 
                                            "MeanAUCSE" = mean(df$SE), 
                                            "EmpAUCSE"  = sd(df$AUCest), 
                                            "beta" = mean(df$beta),  
                                            "MeanbetaEST" = mean(df$betaest), 
                                            "MeanbetaSE" = mean(df$betaSE), 
                                            "EmpbetaSE" = sd(df$betaest)))
AUC MeanAUCest  MeanAUCSE   EmpAUCSE   beta  MeanbetaEST MeanbetaSE EmpbetaSE
1 0.5  0.4974704 0.05044440 0.03847457 0.0000 -0.001571789  0.1609316 0.1287296
2 0.6  0.5991513 0.04956048 0.03638392 0.3247  0.339899843  0.1628435 0.1248197
3 0.7  0.6954977 0.04872714 0.03352307 0.6759  0.690574201  0.1711999 0.1347761
4 0.8  0.7931464 0.04302921 0.02532747 1.1220  1.140682976  0.1948035 0.1367206
5 0.9  0.8905602 0.02970482 0.01418483 1.9387  1.957841713  0.2609241 0.1652734


ddply(cohort.estimates, "AUC", function(df)c("MeanAUCest" = mean(df$AUCest), 
                                             "MeanAUCSE" = mean(df$SE), 
                                             "EmpAUCSE"  = sd(df$AUCest), 
                                             "beta" = mean(df$beta),  
                                             "MeanbetaEST" = mean(df$betaest), 
                                             "MeanbetaSE" = mean(df$betaSE), 
                                             "EmpbetaSE" = sd(df$betaest)))


  AUC MeanAUCest  MeanAUCSE   EmpAUCSE   beta  MeanbetaEST MeanbetaSE  EmpbetaSE
1 0.5  0.4980593 0.02546493 0.02610899 0.0000 -0.003161961 0.08171279 0.08499781
2 0.6  0.6002019 0.02473490 0.02500902 0.3247  0.331674320 0.08282186 0.08386807
3 0.7  0.6974502 0.02265945 0.02356091 0.6759  0.676506669 0.08618334 0.09093153
4 0.8  0.7983985 0.01801694 0.01862269 1.1220  1.132438389 0.09597059 0.09784758
5 0.9  0.8970353 0.01084510 0.01068358 1.9387  1.947708803 0.12473059 0.12548687


