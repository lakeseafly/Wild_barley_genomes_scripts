### This script is edited and adapted from Claire Mérot's script: https://github.com/clairemerot/angsd_pipeline/blob/master/01_scripts/Rscripts/Hobs_sliding.r

library(windowscanr)

argv <- commandArgs(T)
FILE <- argv[1]
WINDOW<- as.numeric(argv[2])
WINDOW_STEP<- as.numeric(argv[3])

print(paste("loading",FILE))
hwe<-read.table(FILE, header=T)

#output the same hwe matrix results but with estimated Hobs - big file - value are for each position
write.table(hwe, paste0(FILE,".Hobs"), row.names=F, quote=F, sep="\t")

#do sliding-window (time consuming)
print(paste("compute sliding window of size",WINDOW, "with step",WINDOW_STEP))
hwe_win <- winScan(x = hwe,groups = "Chromo", position = "Position",values = c("Hobs", "F"),win_size = WINDOW,win_step = WINDOW_STEP,funs = c("mean"))
head(hwe_win)

#output the sliding-windows file
write.table(hwe_win, paste0(FILE,".slidingwindows"), row.names=FALSE, quote=FALSE, sep="\t")

