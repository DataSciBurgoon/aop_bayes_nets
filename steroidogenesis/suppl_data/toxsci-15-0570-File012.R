##----------------------------------------------------------------------------##
## R Script for Analyses comprising Tables and Figures for the manuscript
## "High-throughput screening of chemical effects on steroidogenesis in H295R
##  human adrenocortical carcinoma cells"
##
## Authors -- Agnes Karmaus & Matt Martin
## Emails  -- karmaus.agnes@epa.gov, martin.matt@epa.gov
##
## October 21, 2015
##----------------------------------------------------------------------------##

library(data.table)

setwd("set-working-directory-to-folder-containig-supplementary-files")

mtc.mtt <- fread(file.path("Supplemental_file_1.csv"))
mtc0 <- fread(file.path("Supplemental_file_2.csv"))
mtc <- fread(file.path("Supplemental_file_3.csv"))

cr.mtt <- fread(file.path("Supplemental_file_4.csv"))
cr0 <- fread(file.path("Supplemental_file_5.csv"))
cr <- fread(file.path("Supplemental_file_6.csv"))


################################################################################
################################################################################
### sTART: CALCULATE MANUSCRIPT TABLE 2 DATA

#-------------------------------------------------------------------------------
# Calculating selection sensitivity
#-------------------------------------------------------------------------------

mtc[ , hitcs := sum(hitc), by = list(spid)]
hmtc <- unique(mtc[ , list(spid, chid, hitcs)])

cr[ , hitcm := sum(hitc[hitc > 0]), by = list(spid)]
hcr <- unique(cr[ , list(spid, chid, hitcm)])

setkey(hmtc, spid)
setkey(hcr, spid)
hits.chem <- hmtc[hcr]

h <- hits.chem[, .N , by = list(hitcs, hitcm)]
setkeyv(h, c("hitcs", "hitcm"))

h[ hitcs >= 4 , sum(N[hitcm > 0]) / sum(N)]


#-------------------------------------------------------------------------------
# Correlation of maximum fold change
#-------------------------------------------------------------------------------

mtc[ , assay := sub("_dn", "", aenm)]
mtc[ , assay := sub("_up", "", assay)]
mtc[ , assay := sub("CEETOX_H295R_", "", assay)]
mtc[ , assay_spid := paste(assay, spid, sep = "_")]
mtc1 <- mtc[ , list(resps = max(abs(max_med))),
            by = list(assay_spid, assay, spid)]
mtc1 <- mtc1[!grep(c("DMSO|FOR|PRO"), spid),]

cr[ , assay := sub("_dn", "", aenm)]
cr[ , assay := sub("_up", "", assay)]
cr[ , assay := sub("CEETOX_H295R_", "", assay)]
cr[ , assay_spid := paste(assay, spid, sep = "_")]
cr1 <- cr[ , list(respm = max(abs(max_med))),
             by = list(assay_spid, assay, spid)]

setkeyv(mtc1, c("assay_spid","assay","spid"))
setkeyv(cr1, c("assay_spid","assay","spid"))
max.fc <- merge(mtc1, cr1)
horm.cor <- max.fc[ , list(horm.cor = summary(lm(respm ~ resps))$adj.r.squared),
                by = list(assay)]
horm.cor

#-------------------------------------------------------------------------------
# Calculating summary statistics (Z' values & SSMD values)
#-------------------------------------------------------------------------------

m0 <- unique(mtc0[ , list(acid, acnm, rval, wllt, wllq, apid, rowi, coli,
                          type = "sc")])
c0 <- unique(cr0[ , list(acid, acnm, rval, wllt, wllq, apid, rowi, coli,
                         type = "mc")])
dat0 <- rbind(m0,c0)

dat.qual <- dat0[ , list(
  zprm.p = 1 - ((3 * (
    mad(rval[wllt=="p"], na.rm = TRUE) + 
      mad(rval[wllt=="n"], na.rm = TRUE))) / 
      abs(
        median(rval[wllt == "p"], na.rm = TRUE) - 
          median(rval[wllt == "n"], na.rm = TRUE))),
  zprm.m = 1 - ((3 * (
    mad(rval[wllt == "m"], na.rm = TRUE) + 
      mad(rval[wllt == "n"], na.rm = TRUE))) / 
      abs(
        median(rval[wllt == "m"], na.rm = TRUE) - 
          median(rval[wllt == "n"], na.rm = TRUE))),
  ssmd.p = (median(rval[wllt == "p"], na.rm = TRUE) -
       median(rval[wllt == "n"], na.rm = TRUE)) / 
      sqrt(
        mad(rval[wllt == "p"], na.rm = TRUE)^2 +
          mad(rval[wllt == "n"], na.rm = TRUE)^2 ),
  ssmd.m = (median(rval[wllt == "m"], na.rm = TRUE) -
       median(rval[wllt == "n"], na.rm = TRUE)) / 
      sqrt(
        mad(rval[wllt == "m"], na.rm = TRUE)^2 +
          mad(rval[wllt == "n"], na.rm = TRUE)^2 )
  ), by = list(acid, acnm, apid, type)] 

dat.qual[!is.na(zprm.p) | !is.na(zprm.m),
         zprm := max(c(zprm.p, zprm.m), na.rm = TRUE),
         by = list(acid, acnm, apid, type)]
dat.qual[!is.na(ssmd.p) | !is.na(ssmd.m),
         ssmd := max(c(ssmd.p, ssmd.m), na.rm = TRUE),
         by = list(acid, acnm, apid, type)]

acqu <- dat.qual[ , list(zprm.p = round(median(zprm.p, na.rm = TRUE), 2),
                         zprm.m = round(median(zprm.m, na.rm = TRUE), 2),
                         ssmd.p = round(median(ssmd.p, na.rm = TRUE), 0),
                         ssmd.m = round(median(ssmd.m, na.rm = TRUE), 0)
),
by = list(acid, acnm)]
acqu[zprm.p < 0, zprm.p :=0 ]
acqu[zprm.m < 0, zprm.m :=0 ]
setkey(acqu,"acnm")
acqu

#-------------------------------------------------------------------------------
# Concentration-response BMAD and cutoff values
#-------------------------------------------------------------------------------

cr[ , assay := sub("_dn", "", aenm)]
cr[ , assay := sub("_up", "", assay)]
cr[ , assay := sub("CEETOX_H295R_", "", assay)]

cutoffs <- unique(cr[ , c("assay","bmad","coff"), with = FALSE])
cutoffs[ , coff_fold := 2^coff] ## coff's in log2 (as are all response measures)
cutoffs

### END: CALCULATE MANUSCRIPT TABLE 2 DATA
################################################################################
################################################################################
################################################################################
################################################################################
### START: CREATE MANUSCRIPT FIGURES

library(data.table)
library(gplots)
library(RColorBrewer)

mtc.mtt <- fread(file.path("Supp1_ToxCast_H295R_MTC_MTT_151021.csv"))
mtc0 <- fread(file.path("Supp2_ToxCast_H295R_MTC_RawNormData_151021.csv"))
mtc <- fread(file.path("Supp3_ToxCast_H295R_MTC_ActivityCalls_151021.csv"))

cr.mtt <- fread(file.path("Supp4_ToxCast_H295R_CR_MTT_151021.csv"))
cr0 <- fread(file.path("Supp5_ToxCast_H295R_CR_RawNormConcRespData_151021.csv"))
cr <- fread(file.path("Supp6_ToxCast_H295R_MC_ModelsHitCallsFlags_151021.csv"))


#-------------------------------------------------------------------------------
# Figure 3: Conazoles profiling heatmap
#-------------------------------------------------------------------------------

conazole <- cbind(c('triazole','triazole','triazole','triazole','triazole',
                    'triazole','imidazole','triazole','triazole','triazole',
                    'triazole','triazole','triazole','imidazole','imidazole',
                    'imidazole','imidazole','triazole','triazole','triazole',
                    'imidazole','imidazole','imidazole','imidazole','triazole',
                    'triazole','triazole','triazole','triazole','imidazole',
                    'triazole','triazole'),
                  c('102676-31-3','107534-96-3','112281-77-3','112809-51-5',
                    '114369-43-6','116255-48-2','117337-19-6','119446-68-3',
                    '120511-73-1','125116-23-6','125225-28-7','131983-72-7',
                    '133855-98-8','142459-58-3','23593-75-1','24169-02-6',
                    '35554-44-0','43121-43-3','55219-65-3','60207-90-1',
                    '65277-42-1','67485-29-4','67747-09-5','68694-11-1',
                    '76738-62-0','79983-71-4','83657-24-3','85509-19-9',
                    '86386-73-4','87674-68-8','88671-89-0','94361-06-5'))
colnames(conazole) <- c("chm_cat", "casrn")
conazole <- as.data.table(conazole)

cr[ , assay := sub("_dn", "", aenm)]
cr[ , assay := sub("_up", "", assay)]
cr[ , assay := sub("CEETOX_H295R_", "", assay)]

max.resp <- cr[hitc == 1 , list(assay, aenm, spid, chid, casn, chnm,
                                       modl_ga, max_med, modl_tp)]
max.resp[ , max_resp := max(modl_tp), by = list(chid, assay)]
max.resp[ , ac50 := round(10^modl_ga, digits = 2)]

setkey(conazole, casrn)
setkey(max.resp, casn)
con <- max.resp[conazole]
con <- con[!is.na(assay)]

con[ , mac50 := ac50[max_resp == modl_tp], by = list(chid, assay)]
con[grep("_dn", aenm), resp_tp := -1*max_resp]
con[grep("_up", aenm), resp_tp := max_resp]

con1 <- con[ , list(chid, chnm, chm_cat, assay, resp_tp, mac50)]

hmap.ac50 <- dcast(con1, chid + chnm + chm_cat ~ assay,
                   value.var = "mac50",
                   fun.aggregate = min,
                   fill = 0)
setcolorder(hmap.ac50, c(1,2,3,10,12,11,7,6,4,5,13,9,8))
hmap.mac50 <- data.matrix(hmap.ac50[ , !c("chid", "chnm", "chm_cat"), with = FALSE])
rownames(hmap.mac50) <- hmap.ac50[ , chid]


hmap.dat <- dcast(con1, chid + chnm + chm_cat ~ assay,
                   value.var = "resp_tp",
                   fun.aggregate = min,
                   fill = 0)
setcolorder(hmap.dat, c(1,2,3,10,12,11,7,6,4,5,13,9,8))
hmap.mat <- data.matrix(hmap.dat[ , !c("chid", "chnm", "chm_cat"), with = FALSE])
rownames(hmap.mat) <- hmap.dat[ , chid]

f <- factor(hmap.dat$chm_cat)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-7, -0.7, length = 100),
               seq(-0.69, 0.89, length = 100),
               seq(0.9, 5, length = 100))

heatmap.2(hmap.mat,
          cellnote = hmap.mac50,
          notecol = "white",
          density.info = "none",
          trace = "none",
          margins = c(7, 7),
          dendrogram = "row",
          key = TRUE,
          Colv = NA,
          Rowv = TRUE,
          col = my_palette,
          breaks = col_breaks,
          symbreaks = FALSE,
          RowSideColors = rev(brewer.pal(4,"Set1"))[f])

#-------------------------------------------------------------------------------
# Figure 4: All chemicals profiling k-means heatmap
#-------------------------------------------------------------------------------

max.resp[grep("_dn", aenm), resp_tp := -1*max_resp]
max.resp[grep("_up", aenm), resp_tp := max_resp]

hmap.all <- dcast(max.resp, spid + chid + chnm ~ assay,
                  value.var = "resp_tp",
                  fun.aggregate = max,
                  fill = 0)
hmap.mall <- data.matrix(hmap.all[ , !c("chid", "spid", "chnm"),
                                  with = FALSE])
rownames(hmap.mall) <- hmap.all[ , chid]

fit <- kmeans(hmap.mall, 5)

kmean.dat <- cbind(hmap.all, fit$cluster)
setnames(kmean.dat, "V2", "kcluster")
setcolorder(kmean.dat, c(1,2,3,14,10,12,11,7,6,4,5,13,9,8))
setkey(kmean.dat, "kcluster")

f <- factor(kmean.dat$kcluster)

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 299)
col_breaks = c(seq(-7, -0.7, length = 100),
               seq(-0.69, 0.89, length = 100),
               seq(0.9, 5, length = 100))

heatmap(data.matrix(kmean.dat)[ , 5:14],
        Colv = NA,
        Rowv = NA,
        scale = "none",
        labRow = NA,
        col = my_palette,
        breaks = col_breaks,
        margins = c(7,10),
        RowSideColors = rev(brewer.pal(5,"Set1"))[f])


### END: CREATE MANUSCRIPT FIGURES
################################################################################
################################################################################

