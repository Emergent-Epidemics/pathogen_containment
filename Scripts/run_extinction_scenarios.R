#Calculate extinction probs
#SV Scarpino
#Jan 2023

#https://github.com/cmmid/ringbp
#Hellewell et al. 2020 https://www.thelancet.com/article/S2214-109X(20)30074-7/fulltext
#Hebert-dufresne et al. 2020 https://royalsocietypublishing.org/doi/10.1098/rsif.2020.0393

###########
#libraries#
###########
library(reticulate)
library(googlesheets4)
library(RColorBrewer)
library(lamW)
library(wesanderson)
library(ringbp)
library(ggplot2)

#########
#Globals#
#########
save_new <- FALSE
run_new <- FALSE
timestamp <- as.numeric(Sys.time())

######
#Data#
######


#######
#Model#
#######
source_python('compute_final_p.py')

##########
#Analysis#
##########
if(run_new == TRUE){
  k <- c(seq(from = 0.01, to = 2.5, length.out = 50))
  
  R0 <- seq(0.99, 5.5, length.out = 50)
  
  traced <- seq(0.001, 0.999, length.out = 100)
  
  stopped <- c()
  R0s <- c()
  ks <- c()
  c95 <- c()
  pb <- txtProgressBar(1, length(R0), style = 3)
  for(i in 1:length(R0)){
    for(j in 1:length(k)){
      res.ij <- ringbp::scenario_sim(n.sim = 100, num.initial.cases = 1, prop.asym=0, prop.ascertain = 0.95, cap_cases = 500, cap_max_days = 150, r0isolated = 0, r0community = R0[i], disp.com = k[j], disp.iso = 1, delay_shape = 1.651524, delay_scale = 4.287786, k = 0, quarantine = FALSE)
      contain.ij <- ringbp::extinct_prob(res.ij,cap_cases = 500)
      c95 <- c(c95, contain.ij)
      solved.ij <- solve_for_z(traced*R0[i], k[j])
      o0.01.ij <- which(solved.ij > 0.01)
      if(length(o0.01.ij) == 0){
        stopped.ij <- 0
      }else{
        if(length(which(solved.ij < 0.01)) == 0){
          stopped.ij <- 1
        }else{
          stopped.ij <- 1-traced[min(o0.01.ij)]
        }
      }
      stopped <- c(stopped, stopped.ij)
      R0s <- c(R0s, R0[i])
      ks <- c(ks, k[j])
    }
    setTxtProgressBar(pb, i)
  }
  
  
  c95[which(is.na(c95) == TRUE)] <- 1
  dat.plot <- data.frame(R0s, ks, stopped, c95)
}else{
  dat.plot <- read.csv("../Data/control_sims_hellewell_hebert_Jan1923.csv")
}

#######
#Plots#
#######
pal <- rev(wes_palette(name = "Zissou1", 10, type = "continuous"))
ggplot(dat.plot, aes(x = R0s, y = ks,  z = c95)) + geom_contour_filled() + scale_y_log10() + xlab("R0") + ylab("Overdispersion (log-scale)") + theme(legend.position = "right", legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffff75", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_fill_manual(values = pal, name = "Percent controlled")

ggplot(dat.plot, aes(x = (1-stopped), y = c95, color = ks)) + geom_point() + xlab("Hebert-Dufresne et al. (1 - final size)") + ylab("Hellewell et al. (Extinction prob. at 95% containment)") + theme(legend.position = c(0.8,0.2), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffff75", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_color_gradientn(colors = pal, name = "Overdispersion (k)") + geom_abline(slope = 1, intercept = 0)

###########
#Save Data#
###########
if(save_new == TRUE){
  filename <- paste0("../Output/", timestamp, "control_sims_hellewell_hebert.csv")
  write.csv(dat.plot, file = filename, row.names = FALSE)
}
