library(hierfstat)
library(adegenet)
library(PopGenReport)
library(plotrix)
library(poppr)
library(tidyr)

setwd("~/snails/L_virgata/population-genomics/virgata_m3M3n3_R80_maf025/")
sumstats<-read.table("populations.sumstats.tsv", header=FALSE)
sumstats<-sumstats[!grepl("#",sumstats),] #Gets rid of the needless rows added by STACKS

#Split by population and keep in order of original popmap and what the populations output files are
sumstats$V5<-factor(sumstats$V5, levels=unique(sumstats$V5))
sumstats_list<-split(sumstats,sumstats$V5)

#calculate averages for each population
stats_table=data.frame()
for (i in sumstats_list){
  pop_data<-as.data.frame(i)
  pop<-pop_data[1,5]
  Ho<-round(mean(pop_data[,10]),digits=4)
  Ho_SE<-round(std.error(pop_data[,10]),digits=4)
  He<-round(mean(pop_data[,12]),digits=4)
  He_SE<-round(std.error(pop_data[,12]),digits=4)
  Pi<-round(mean(pop_data[,14]),digits=4)
  Pi_SE<-round(std.error(pop_data[,14]),digits=4)
  Fis<-round(mean(pop_data[,17]),digits=4)
  Fis_SE<-round(std.error(pop_data[,17]),digits=4)
  
  ##plot FIS values and write as pdf
  name_plot<-paste(pop,".FIS-plot.pdf",sep="")
  pdf(name_plot)
  plot(pop_data[,17], main = pop, xlab = "locus", ylab = "Fis")
  abline(h=Fis, col="red")
  dev.off()
  
  stats_table<-rbind(stats_table,list(pop,Ho,Ho_SE,He,He_SE,Pi,Pi_SE,Fis,Fis_SE))
}

##Allelic Richness##

##Function for calculating standard error on all columns
colSe <- function (x, na.rm = TRUE) {
  if (na.rm) {
    n <- colSums(!is.na(x))
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x*x, na.rm = na.rm) - (colMeans(x, na.rm = na.rm))^2
  return(sqrt(colVar/n))
}


data1<-read.genepop("populations.haps.gen")

data_fstat<-genind2hierfstat(data1)
data_fstat$pop
AR<-allelic.richness(data_fstat) ##Didn't use min.n because this one was the lowest, so default is OK.
AR$min.all
AR_mean<-colMeans(AR$Ar,na.rm=TRUE)
AR_SE<-colSe(AR$Ar,na.rm=TRUE)
AR_mean
AR_SE


###Make Table###
# Fix column names on stats_table first
colnames(stats_table) <- c("Site", "Ho", "Ho_SE", "He", "He_SE", "Pi", "Pi_SE", "Fis", "Fis_SE")

# Force AR names to match stats_table site order exactly
AR_df <- data.frame(
  Site  = stats_table$Site,   # use stats_table names directly
  AR    = round(AR_mean, 4),
  AR_SE = round(AR_SE, 4)
)

# Merge by Site
combined_table <- merge(stats_table, AR_df, by = "Site", sort = FALSE)

# Format all stats as "mean (SE)"
final_table <- data.frame(
  Site = combined_table$Site,
  `Ho`  = paste0(combined_table$Ho,  " (", combined_table$Ho_SE,  ")"),
  `He`  = paste0(combined_table$He,  " (", combined_table$He_SE,  ")"),
  `Pi`  = paste0(combined_table$Pi,  " (", combined_table$Pi_SE,  ")"),
  `Fis` = paste0(combined_table$Fis, " (", combined_table$Fis_SE, ")"),
  `AR`  = paste0(combined_table$AR,  " (", combined_table$AR_SE,  ")"),
  check.names = FALSE
)

print(final_table, row.names = FALSE)

# Optional: export
write.csv(final_table, "diversity_stats_table.csv", row.names = FALSE)
