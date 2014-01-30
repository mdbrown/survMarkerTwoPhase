## make some tables from the package validation simulations. 
setwd("inst/packageValidation/")

library(reshape2)
library(plyr)
library(ggplot2)
require(grid)

beta = log(3)
lam0 = 0.1
f_x = dnorm
cutoff = 0
predict.time = 2

#need functions to calculate the true parameter values
source("../../../../power/PowerSAM/subroutines.R")
source("../../../PowerSAM/subroutines.R")
#truevalues 
truevals <- data.frame(coef = beta, 
                       AUC = get.AUC.given.beta.intmethod(beta = beta, a = lam0, t=predict.time, f_x = f_x), 
                       TPR = TPF.fun(c=cutoff, beta = beta, a = lam0, t = predict.time, f_x = f_x), 
                       FPR = FPF.fun(c=cutoff, beta = beta, a = lam0, t = predict.time, f_x = f_x), 
                       PPV = PPV.fun(c=cutoff, beta = beta, a = lam0, t = predict.time, f_x = f_x), 
                       NPV = NPV.fun(c=cutoff, beta = beta, a = lam0, t = predict.time, f_x = f_x))

truevals.long <- melt(truevals)
####
#cch
####

##np
load("cch.np.sim.Rdata")

##sp
load("cch.sp.sim.Rdata")
cchData <- data.frame(matrix(nrow=1000*2, ncol = 7))
names(cchData) = c("coef", "AUC", "TPR", "FPR", "PPV", "NPV", "method")

cchData[1:1000,-c(7)] <- data.frame(t(matrix(cch.sp.sim, nrow = 6)))
cchData[-c(1:1000),-c(1,7)] <- data.frame(t(matrix(cch.np.sim, nrow = 5)))
cchData[,7] <- factor(rep(c("semiparametric", "nonparametric"), c(1000, 1000)))


cchData.long <- melt(cchData)

mdat <- ddply(cchData.long, .(variable, method), summarize, mean = mean(value))
tvlong <- melt(truevals)
tvlong$method = rep("truth", 6)
names(tvlong) = c("variable", "mean", "method"); tvlong <- tvlong[,c(1,3,2)]

mdat <- rbind(mdat, tvlong)

p <- ggplot(cchData.long, aes(x=value, fill = method))
p <- p + geom_density(alpha = .5) + facet_grid(~variable, scales="free_x")
p <- p + geom_vline(data = truevals.long, aes(xintercept=value), size = 1,
                    linetype = 2, colour = "purple")
p <- p+theme(text=element_text(size=25))
p <- p+theme(axis.text.x = element_text(size=13), legend.position= "top", legend.title = element_blank())
p <- p 
print(p)

"#F8766D" "#00BA38" "#619CFF"

data_table <- ggplot(mdat, aes(x = factor(variable), y = factor(method),
                              label = format(round(mean,3), nsmall = 1), colour = factor(method))) +
  geom_text(size = 5, hjust = 0) + theme_bw()  +
  theme(panel.grid.major = element_blank(), legend.position = "none",
       panel.border = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(), 
        plot.margin =unit(c(-0.5,1, 0, 0.5), "line")) + 
       xlab(NULL) + ylab(NULL) + 
       scale_color_manual(values = c("#F8766D" ,"#00BFC4", "purple")) +
       scale_y_discrete(labels = c("NP", "SP", "true"))



 Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2, 0.25), c("null", "null")))
 grid.show.layout(Layout)
 vplayout <- function(...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}

subplot <- function(x, y) viewport(layout.pos.row = x,
                                     layout.pos.col = y)
 mmplot <- function(a, b) {
  vplayout()
  print(a, vp = subplot(1, 1))
  print(b, vp = subplot(2, 1))
}

png(filename = "cch_validation.png", width = 15, height = 8, units = "in", res = 90)
 mmplot(p, data_table)
dev.off()

ggsave("cch_validation.png", width = 15, height = 8 )




####
#ncc
####

#estimates

##np

NPresults <- NULL
for( i in 1:10){
load(paste("ncc.np.sim_", i, ".Rdata", sep = ""))
NPresults <- rbind(NPresults, data.frame(ldply(ncc.np.sim[1,], I)) )

}
##sp
load("ncc.sp.sim.Rdata")
nccData <- data.frame(matrix(nrow=1000*2, ncol = 7))
names(nccData) = c("coef", "AUC", "FPR", "TPR", "NPV", "PPV", "method")

nccData[1:1000,-c(7)] <- data.frame(ldply(ncc.sp.sim[1,], I)) #estimates
#not done yet
nccData[-c(1:1000),-c(1,7)] <- NPresults

nccData[,7] <- factor(rep(c("semiparametric", "nonparametric"), c(1000, 1000)))

nccData.long <- melt(nccData)

ddply(nccData.long, .(variable, method), summarize, mean = mean(value), 
      se = sd(value))

p <- ggplot(nccData.long, aes(x=value, fill = method))
p <- p + geom_density(alpha = .6) + facet_grid(~variable, scales="free_x")

p <- p + geom_vline(data = truevals.long, aes(xintercept=value), size = 2, colour = "purple")
p <- p+theme(text=element_text(size=25))
p <- p+theme(axis.text.x = element_text(size=13))
print(p)

ggsave("ncc_validation.png", width = 15, height = 7 )



#se 

EmpiricalSE.long <- ddply(nccData.long, .(variable, method), summarize, mean = mean(value), 
                          se = sd(value))
NPresults <- NULL
for( i in 1:10){
  load(paste("ncc.np.sim_", i, ".Rdata", sep = ""))
  NPresults <- rbind(NPresults, data.frame(ldply(ncc.np.sim[2,], I)) )
  
}
##sp



load("ncc.sp.sim.Rdata")
nccData <- data.frame(matrix(nrow=1000*2, ncol = 7))
names(nccData) = c("coef", "AUC", "FPR", "TPR", "NPV", "PPV", "method")

nccData[1:1000,-c(7)] <- data.frame(ldply(ncc.sp.sim[2,], I)) #se
#not done yet
nccData[-c(1:1000),-c(1,7)] <- NPresults

nccData[,7] <- factor(rep(c("semiparametric", "nonparametric"), c(1000, 1000)))

nccData.long <- melt(nccData)

cbind( ddply(nccData.long, .(variable, method), summarize, meanSE = mean(value)), EmpiricalSE.long$se)


p <- ggplot(nccData.long, aes(x=value, fill = method))
p <- p + geom_density(alpha = 1) + facet_grid(method~variable, scales="free_x")

p <- p + geom_vline(data = EmpiricalSE.long, aes(xintercept=se), size = 1.5, colour = "purple", linetype=2)
p <- p+theme(text=element_text(size=25))
p <- p+theme(axis.text.x = element_text(size=13)) #+ ylim(c(0,100))
print(p)

ggsave("ncc_validation_SE.png", width = 15, height = 11 )




