library(poppr)
library(ade4)

setwd("~/snails/L_virgata/population-genomics/virgata_m3M3n3_R80_maf025/")
rm(list=ls())

##Read Data##
virgata.geneid = read.genepop("populations.snps.gen")
virgata.genclone<-as.genclone(virgata.geneid)
virgata.genclone$pop
strata(virgata.genclone)<-(as.data.frame(virgata.geneid$pop))
virgata.genclone$strata

##AMOVA ade4##
virgata.site.amova.pop = poppr.amova(virgata.genclone, ~virgata.geneid.pop, cutoff = 0.5, within=FALSE, method = "ade4")
virgata.site.amova.pop

#Randomization Tests for ade4
virgata.site.amova.pop.rtest<-randtest(virgata.site.amova.pop,nrepet = 1000)
virgata.site.amova.pop.rtest
plot(virgata.site.amova.pop.rtest)

##AMOVA PEGAS##
virgata.site.amova.pop.pegas = poppr.amova(virgata.genclone, ~virgata.geneid.pop, cutoff = 0.5, nperm = 1000, within=FALSE, method = "pegas")
virgata.site.amova.pop.pegas
