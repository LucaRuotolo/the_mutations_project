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
pathB <- "C:/Pat.B Table of Ct.csv"  
Patient_B <- "C:/Pat.B Stanford Analysis.csv"

#--- PATIENT B ---- 
  
  #--- CT PLOT ---- #
  ctB <- data.table(fread(pathB, sep=";", na="."))
  colnames(ctB) <- c("DAY","E","RDRP","N")
  ctB$DAY <- as.numeric(ctB$DAY)
  ctB$E <- as.numeric(ctB$E)
  ctB$RDRP <- as.numeric(ctB$RDRP)
  ctB$N <- as.numeric(ctB$N)
  
  plotctB <- ggplot(ctB)+geom_point(data = ctB[!is.na(N)],aes(x=DAY, y=N))+
  geom_path(data = ctB[!is.na(N)], aes(x=DAY, y=N, colour="N gene"),col="blue", size=0.5)+
  geom_hline(yintercept = 42, linetype=4, size=1)+
  geom_vline(xintercept = 1, linetype=2, col="red")+
  geom_vline(xintercept = 19, linetype=2, col="red")+
  geom_vline(xintercept = 38, linetype=2, col="red")+
  scale_y_continuous(breaks=c(1:42),trans = "reverse")+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60))+
  labs(x="Days", y="Ct value", col="Target")+
  theme_bw()
  
  ggsave(paste0("C:/Figure 3.jpeg"), plotctB, type = "cairo-png", scale=.5, width=60, height=30,  units="cm")
  
  
  #--- MUTATION ANALYSIS ----#
  
  heat_B <-  data.table(fread(Patient_B, sep=";", na="."))
  heat_B <- heat_B[, c("MUTATION", "GENE", "PERC", "LINEAGE", "TP", "DELTA")]
  heat_B[GENE=="S"]$MUTATION <- paste("S:",heat_B[GENE=="S"]$MUTATION,"")
  heat_B$MUTATION <- str_remove_all(heat_B$MUTATION, " ")
  
  hzmod <- heat_B %>%
    group_by(MUTATION) %>%
    count(TP==1 | TP==2, name="nTP")            #count number of time points in report
  
  hzmod <- data.table(hzmod)
  hzmod <- heat_B[heat_B$MUTATION %in% hzmod[nTP==1]$MUTATION]
  
  insertp <- data.table(MUTATION=hzmod$MUTATION, GENE=hzmod$GENE, PERC=0, LINEAGE=NA, TP=NA, DELTA=NA)
  heat_B <- data.table(setorder(rbind(heat_B, insertp, fill=T), MUTATION))  #insert missing time points
  heat_B <- fill(heat_B, TP, .direction="down")
  
  heat_B[is.na(LINEAGE)]$TP <- heat_B[is.na(LINEAGE)]$TP%%2+1
  heat_B <- data.table(setorder(heat_B, MUTATION, TP))
  heat_B[is.na(LINEAGE)]$LINEAGE <- levels(factor(heat_B$LINEAGE))[heat_B[is.na(LINEAGE)]$TP] #fill missing time point record
  
  delt <- heat_B %>% 
    group_by(MUTATION) %>% 
    summarise(maxp=max(PERC), dif=abs(PERC-maxp),.groups = "drop")    #estimate delta value for each time mutation in each TP
  delt <- data.table(delt)    
  
  heat_B$DELTA <- delt$dif
  heat_Bfilt <- heat_B[heat_B$MUTATION %in% delt[dif>=cutoff]$MUTATION] #Pat.B dataset filtered with delta > cutoff
  
  heatmap_B <- ggplot(heat_B, aes(x=MUTATION, y=TP,fill=PERC))+ 
  geom_tile(color = "white", size=.5)+theme_bw()+
  scale_fill_gradient(high = "dodgerblue4", low = "lightcyan", labels=c(0, 25, 50,75,100))+
  scale_y_continuous(breaks = c(1,2), labels = c("day1", "day86"), trans = "reverse")+
  scale_x_discrete(position = "top")+
  labs(x="Mutation", y="Time Points", fill="Frequency (%)")+
  theme(legend.position="right", axis.text.y = element_text(face="bold"))+ 
  guides(x = guide_axis(angle = 90))
  
  ggsave(paste0("C:/Figure 4A.jpeg"), heatmap_B, type = "cairo-png", scale=.5, width=60, height=20,  units="cm")
  
  heatmap_B5 <- ggplot(heat_Bfilt, aes(x=MUTATION, y=TP,fill=PERC))+ 
  geom_tile(color = "white", size=.5)+
  theme_bw()+
  scale_fill_gradient(high = "dodgerblue4", low = "lightcyan", labels=c(0, 25, 50,75,100))+
  scale_y_continuous(breaks = c(1,2), labels = c("day1", "day86"), trans = "reverse")+
  scale_x_discrete(position = "top")+labs(x="Mutation", y="Time Points", fill="Frequency (%)")+
  theme(legend.position="right", axis.text.y = element_text(face="bold"))+ 
  guides(x = guide_axis(angle = 90))
  
  ggsave(paste0("C:/Figure 4B.jpeg"), heatmap_B5, type = "cairo-png", scale=.5, width=60, height=20,  units="cm")