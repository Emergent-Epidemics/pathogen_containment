#Compare extinction probs.
#SV Scarpino
#Jan 2023

#https://github.com/cmmid/ringbp
#Hellewell et al. 2020 https://www.thelancet.com/article/S2214-109X(20)30074-7/fulltext
#Hebert-dufresne et al. 2020 https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0393ks <- c(seq(from = 0.001, to = 1.01, length.out = 50))

###########
#libraries#
###########
library(ringbp)

#########
#Globals#
#########


######
#Data#
######


##########
#Analysis#
##########
R0 <- 2.5
exps <- c()
pb <- txtProgressBar(1, length(ks), style = 3)
for(i in 1:length(ks)){
  res.i <- ringbp::scenario_sim(n.sim = 100, num.initial.cases = 1,prop.asym=0,prop.ascertain = 0.95, cap_cases = 50, cap_max_days = 150, r0isolated = 0, r0community = 1.5, disp.com = ks[i], disp.iso = 1, delay_shape = 1.651524,delay_scale = 4.287786,k = 0, quarantine = FALSE)
  
  exp.i <- ringbp::extinct_prob(res.i,cap_cases = 50)
  
  exps <- c(exps, exp.i)
  setTxtProgressBar(pb, i)
}

solved <- solve_for_z(R0, ks)

##################
#Plots and models#
##################
y <- 1 - exps
plot(1-solved, exps, xlab = "Hebert-dufresne et al. (1-final size)", ylab = "Hellewell et al. (Extinction prob. 95% containment)", main = "R0 = 1.5, k = 0.001 - 1", pch = 16, bty = "n", xlim = c(0, 1), ylim = c(0, 1))
abline(0, 1, lwd = 3, col = "#4d4d4d", lty = 3)

summary(lm(y ~ solved))
        