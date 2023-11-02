library(data.table)
library(dplyr)
library(gridExtra)
library(readxl)
library(zoo)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(gganimate)
library(lubridate)
library(tidyr)
library(stringr)
library(stringi)
library(cobs)

cutoff <- 5                           #Set the cutoff value for Minimum Allele Frequency
pathA <- "C:/Pat.A Table of Ct.csv"  
Patient_A <- "C:/Pat.A Stanford Analysis.csv"

#--- PATIENT A ---- 
  
  #--- CT PLOT ---- #
  ctA <- data.table(fread(pathA, sep=";", na="."))
  colnames(ctA) <- c("DAY","E","RDRP","N")
  ctA$DAY <- as.numeric(ctA$DAY)
  ctA$E <- as.numeric(ctA$E)
  ctA$RDRP <- as.numeric(ctA$RDRP)
  ctA$N <- as.numeric(ctA$N)
  
  plotctA <- ggplot(ctA)+geom_point(data = ctA[!is.na(N)],aes(x=DAY, y=N))+
  geom_path(data = ctA[!is.na(N)], aes(x=DAY, y=N, colour="N gene"),col="blue", size=0.5)+
  geom_hline(yintercept = 42, linetype=4, size=1)+
  geom_vline(xintercept = 1, linetype=2, col="red")+
  geom_vline(xintercept = 19, linetype=2, col="red")+
  geom_vline(xintercept = 38, linetype=2, col="red")+
  scale_y_continuous(breaks=c(1:42),trans = "reverse")+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60))+
  labs(x="Days", y="Ct value", col="Target")+
  theme_bw()
  
  ggsave(paste0("C:/Figure 1.jpeg"), plotctA, type = "cairo-png", scale=.5, width=60, height=30,  units="cm")
  
  
  #--- MUTATION ANALYSIS ----#
  
  heat_A <-  data.table(fread(Patient_A, sep=";", na="."))
  heat_A <- heat_A[, c("MUTATION", "GENE", "PERC", "LINEAGE", "TP", "DELTA")]
  hmmod <- heat_A %>%
    group_by(MUTATION) %>%
    count(TP==1 | TP==2 | TP==3, name="nTP")               #count number of time points in the report
  
  hmmod <- data.table(hmmod)
  hmmod <- heat_A[heat_A$MUTATION %in% hmmod[nTP==1]$MUTATION] 
  
  #generate missing time points
  insertp <- data.table(MUTATION=hmmod$MUTATION, GENE=hmmod$GENE, PERC=0, LINEAGE=NA, TP=NA, DELTA=NA) 
  
  heat_A <- data.table(setorder(rbind(heat_A, insertp, insertp, fill=T), MUTATION)) #insert missing time points
  heat_A <- fill(heat_A, TP, .direction="down")
  heat_A[is.na(LINEAGE)]$TP <- heat_A[is.na(LINEAGE)]$TP%%3+1 #name each missing time points (1, 2 or 3)
  
  heat_A <- heat_A %>%
    group_by(MUTATION) %>%
    mutate(dup=duplicated(TP))
  
  heat_A <- data.table(heat_A)
  heat_A[dup==T]$TP <- (heat_A[dup==T]$TP)%%3+1               #change value of duplicated time points
  heat_A <- data.table(setorder(heat_A[,-7], MUTATION, TP))
  heat_A <- fill(heat_A, LINEAGE, .direction="down")          
  
  delt <- heat_A %>% 
    group_by(MUTATION) %>% 
    summarise(maxp=max(PERC), dif=abs(PERC-maxp),.groups = "drop") #estimate delta value for each mutation in each time point
  
  delt <- data.table(delt)
  heat_A$DELTA <- delt$dif
  heat_Afilt <- heat_A[heat_A$MUTATION %in% delt[dif>=cutoff]$MUTATION] #filtered dataset of Pat.A with delta > cutoff
  
  heatmap_A <- ggplot(heat_A, aes(x=MUTATION, y=TP,fill=PERC))+ 
  geom_tile(color = "white", size=.5)+theme_bw()+scale_fill_gradient(high = "dodgerblue4", low = "lightcyan", labels=c(0, 25, 50,75,100))+
  scale_y_continuous(breaks = c(1,2,3), labels = c("day1", "day19", "day38"), trans = "reverse")+
  scale_x_discrete(position = "top")+
  labs(x="Mutation", y="Time Points", fill="Frequency (%)")+theme(legend.position="right", axis.text.y = element_text(face="bold"))+ 
  guides(x = guide_axis(angle = 90))
  
  ggsave(paste0("C:/Figure2A.jpeg"), heatmap_A, type = "cairo-png", scale=.5, width=40, height=20,  units="cm")
  
  heatmapA_5 <- ggplot(heat_Afilt, aes(x=MUTATION, y=TP,fill=PERC))+ 
  geom_tile(color = "white", size=.5)+theme_bw()+
  scale_fill_gradient(high = "dodgerblue4", low = "lightcyan", labels=c(0, 25, 50,75,100))+
  scale_y_continuous(breaks = c(1,2,3), labels = c("day1", "day19", "day38"), trans = "reverse")+
  scale_x_discrete(position = "top")+
  labs(x="Mutation", y="Time Points", fill="Frequency (%)")+
  theme(legend.position="right", axis.text.y = element_text(face="bold"))+ 
  guides(x = guide_axis(angle = 90))
  
  ggsave(paste0("C:/Figure2B.jpeg"), heatmap_A5, type = "cairo-png", scale=.5, width=40, height=20,  units="cm")
