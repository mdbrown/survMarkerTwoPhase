# just a bunch of random code to test ideas out as I write this package

data(SimData)


status <- eval(substitute(status), SimData)


#simulate some NCC data

tmpData <- simulateData.ncc(nn=500, 
                            beta = log(3), 
                            lam0 = .1, #baseline hazard 
                            cens.perc = .8, 
                            time.max = NULL, 
                            m.match = 3)
head(tmpData$data)
head(tmpData$Vi.k0.IND)
head(tmpData$Iik0)


cohortData <- data.frame(id = 1:nrow(tmpData$data), tmpData$data)
riskSets <- data.frame(tmpData[[2]][,-c(1:2)])


tmp <- survMTP.ncc(time = xi, 
         event = di, 
         marker=yi, 
         subcoh=vi, 
         id = id, 
         data = cohortData,
         risk.sets = riskSets, 
         estimation.method = "NP", 
         predict.time = 1 )
           
tmp.sp <- survMTP.ncc(time = xi, 
                   event = di, 
                   marker=yi, 
                   subcoh=vi, 
                   id = id, 
                   data = cohortData,
                   risk.sets = riskSets, 
                   estimation.method = "SP", 
                   predict.time = 1 )
tmp.sp
