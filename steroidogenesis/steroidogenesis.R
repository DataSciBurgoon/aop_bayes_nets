library(gRain)
library(data.table)

#Building the Bayes Net
yn <- c("yes", "no")
cholesterol <- cptable(~ cholesterol, values=c(0.95, 0.05), levels=yn)
cyp11a1 <- cptable(~ cyp11a1, values=c(0.50, 0.50), levels=yn)
pregnenolone <- cptable(~ pregnenolone | cyp11a1+cholesterol, values=c(0.99, 0.01, 0.01, 0.01,
                                                                       0.01, 0.99, 0.99, 0.99), levels=yn)
cyp17a1_hydroxylase <- cptable(~ cyp17a1_hydroxylase, values=c(0.50, 0.50), levels=yn)
cyp17a1_lyase <- cptable(~ cyp17a1_lyase, values=c(.50, 0.50), levels=yn)
hsd3b1 <- cptable(~ hsd3b1, values=c(0.50, 0.50), levels=yn)
hsd17b3 <- cptable(~ hsd17b3, values=c(0.50, 0.50), levels=yn)
cyp19a1 <- cptable(~ cyp19a1, values=c(0.50, 0.50), levels=yn)
cyp21a1 <- cptable(~ cyp21a1, values=c(0.50, 0.50), levels=yn)
cyp11b1 <- cptable(~ cyp11b1, values=c(0.50, 0.50), levels=yn)
hydroxypregnenolone <- cptable(~ hydroxypregnenolone | cyp17a1_hydroxylase+pregnenolone, 
                               values=c(0.99, 0.01, 0.01, 0.99,
                                        0.01, 0.99, 0.01, 0.99), 
                               levels=yn)
dhea <- cptable(~ dhea | cyp17a1_lyase+hydroxypregnenolone,
                values = c(0.99, 0.01, 0.01, 0.99,
                           0.01, 0.99, 0.01, 0.99),
                levels=yn)
progesterone <- cptable(~ progesterone | hsd3b1+pregnenolone,
                        values = c(0.99, 0.01, 0.01, 0.99,
                                   0.01, 0.99, 0.01, 0.99),
                        levels=yn)
hydroxyprogesterone <- cptable(~hydroxyprogesterone | cyp17a1_hydroxylase+progesterone,
                               values = c(0.99, 0.01, 0.01, 0.99,
                                          0.01, 0.99, 0.01, 0.99),
                               levels=yn)
# androsteinedione <- cptable(~ androsteinedione | hsd3b1+dhea+cyp17a1_lyase+hydroxyprogesterone,
#                             values = c(0.99, 0.50, 0.50, 0.50, 0.50, 0.01, 0.01, 0.01, 0.50, 0.01, 0.01, 0.01, 0.50, 0.01, 0.01, 0.01,
#                                        0.01, 0.50, 0.50, 0.50, 0.50, 0.99, 0.99, 0.99, 0.50, 0.99, 0.99, 0.99, 0.50, 0.99, 0.99, 0.99),
#                             levels=yn)
androsteinedione1 <- cptable(~ androsteinedione1 | hsd3b1+dhea,
                            values = c(0.99, 0.01, 0.01, 0.99,
                                       0.01, 0.99, 0.01, 0.99),
                            levels=yn)
androsteinedione2 <- cptable(~ androsteinedione2 | cyp17a1_lyase+hydroxyprogesterone,
                             values = c(0.99, 0.01, 0.01, 0.99,
                                        0.01, 0.99, 0.01, 0.99),
                             levels=yn)
androsteinedione <- cptable(~ androsteinedione | androsteinedione1+androsteinedione2,
                            values = c(0.99, 0.01, 0.99, 0.01,
                                       0.99, 0.01, 0.01, 0.99),
                            levels=yn)
estrone <- cptable(~ estrone | cyp19a1+androsteinedione,
                   values = c(0.99, 0.01, 0.01, 0.99,
                              0.01, 0.99, 0.01, 0.99),
                   levels=yn)
testosterone <- cptable(~ testosterone | hsd17b3+androsteinedione,
                        values=c(0.99, 0.01, 0.01, 0.99,
                                 0.01, 0.99, 0.01, 0.99),
                        levels=yn)
estradiol <- cptable(~ estradiol | cyp19a1+testosterone,
                     values=c(0.99, 0.01, 0.01, 0.99,
                              0.01, 0.99, 0.01, 0.99),
                     levels=yn)
deoxycortisone <- cptable(~ deoxycortisone | cyp21a1+progesterone,
                          values=c(0.99, 0.01, 0.01, 0.99,
                                   0.01, 0.99, 0.01, 0.99),
                          levels=yn)
corticosterone <- cptable(~ corticosterone | cyp11b1+deoxycortisone,
                          values=c(0.99, 0.01, 0.01, 0.99,
                                   0.01, 0.99, 0.01, 0.99),
                          levels=yn)
deoxycortisol <- cptable(~ deoxycortisol | cyp21a1+hydroxyprogesterone,
                         values=c(0.99, 0.01, 0.01, 0.99,
                                  0.01, 0.99, 0.01, 0.99),
                         levels=yn)
cortisol <- cptable(~ cortisol | cyp11b1+deoxycortisol,
                    values=c(0.99, 0.01, 0.01, 0.99,
                             0.01, 0.99, 0.01, 0.99),
                    levels=yn)


plist <- compileCPT(list(cholesterol, cyp11a1, pregnenolone, cyp17a1_hydroxylase, cyp17a1_lyase, cyp21a1, cyp11b1, hsd3b1, hsd17b3, cyp19a1, hydroxypregnenolone,
                         dhea, progesterone, hydroxyprogesterone, androsteinedione1, androsteinedione2, androsteinedione, estrone, testosterone, estradiol,
                         deoxycortisone, corticosterone, deoxycortisol, cortisol))
bn_test <- grain(plist)
png("steroidogenesis_network.png", width=5000, height=5000, res=800, pointsize=20)
plot(bn_test)
dev.off()

#Test Queries
querygrain(setEvidence(bn_test, evidence=list(dhea="no", androsteinedione = "no", testosterone="no")))

querygrain(setEvidence(bn_test, evidence=list(cholesterol="yes", pregnenolone="yes")))


#Queries Using Agnes' data

format_results <- function(x){
  xdf <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T))
  rownames(xdf) <- names(x)
  colnames(xdf) <- c("yes", "no")
  return(xdf)
}

#Ketoconazole
ketoconazole_results <- querygrain(setEvidence(bn_test, evidence=list(hydroxyprogesterone="no",
                                              androsteinedione="no", testosterone="no",
                                              estrone="no", deoxycortisone="yes",
                                              cortisol="no", deoxycortisol="no")))
ketoconazole_results_out <- format_results(ketoconazole_results)
write.table(ketoconazole_results_out, file="ketoconazole_results.txt", sep="\t", row.names=TRUE, col.names=NA)


#Imazalil
imazalil_results <- querygrain(setEvidence(bn_test, evidence=list(progesterone="yes", hydroxyprogesterone="yes",
                                              cortisol="no", deoxycortisol="no", 
                                              androsteinedione="no", testosterone="no",
                                              estrone="no")))
imazalil_results_out <- format_results(imazalil_results)
write.table(imazalil_results_out, file="imazalil_results.txt", sep="\t", row.names=TRUE, col.names=NA)


#Fenbuconazole
fb_results <- querygrain(setEvidence(bn_test, evidence=list(progesterone="yes", hydroxyprogesterone="yes",
                                                            cortisol="no")))
fb_results_out <- format_results(fb_results)
write.table(fb_results_out, file="fenbuconazole_results.txt", sep="\t", row.names=TRUE, col.names=NA)


#Difenoconazole
dif_results <- querygrain(setEvidence(bn_test, evidence=list(hydroxyprogesterone="yes", deoxycorticosterone="no",
                                                             cortisol="no")))
dif_results_out <- format_results(dif_results)
write.table(dif_results_out, file="difenoconazole_results.txt", sep="\t", row.names=TRUE, col.names=NA)


#Metconazole
met_results <- querygrain(setEvidence(bn_test, evidence=list(hydroxyprogesterone="no", deoxycorticosterone="no")))
met_results_out <- format_results(met_results)
write.table(met_results_out, file="metconazole_results.txt", sep="\t", row.names=TRUE, col.names=NA)

format_results <- function(x){
  xdf <- data.frame(matrix(unlist(x), nrow=length(x), byrow=T))
  rownames(xdf) <- names(x)
  colnames(xdf) <- c("yes", "no")
  return(xdf)
}


format_evidence_list <- function(x){
  evidence_list <- NULL
  for(i in 1:nrow(x)){
    evidence_list[[as.character(x[i,2])]] <- as.character(x[i,3])
  }
  return(as.list(evidence_list))
}

steroidogenesis_bn <- function(evidence){
  evidence_list <- format_evidence_list(evidence)
  yn <- c("yes", "no")
  cholesterol <- cptable(~ cholesterol, values=c(0.95, 0.05), levels=yn)
  cyp11a1 <- cptable(~ cyp11a1, values=c(0.50, 0.50), levels=yn)
  pregnenolone <- cptable(~ pregnenolone | cyp11a1+cholesterol, values=c(0.99, 0.01, 0.01, 0.01,
                                                                         0.01, 0.99, 0.99, 0.99), levels=yn)
  cyp17a1_hydroxylase <- cptable(~ cyp17a1_hydroxylase, values=c(0.50, 0.50), levels=yn)
  cyp17a1_lyase <- cptable(~ cyp17a1_lyase, values=c(.50, 0.50), levels=yn)
  hsd3b1 <- cptable(~ hsd3b1, values=c(0.50, 0.50), levels=yn)
  hsd17b3 <- cptable(~ hsd17b3, values=c(0.50, 0.50), levels=yn)
  cyp19a1 <- cptable(~ cyp19a1, values=c(0.50, 0.50), levels=yn)
  cyp21a1 <- cptable(~ cyp21a1, values=c(0.50, 0.50), levels=yn)
  cyp11b1 <- cptable(~ cyp11b1, values=c(0.50, 0.50), levels=yn)
  hydroxypregnenolone <- cptable(~ hydroxypregnenolone | cyp17a1_hydroxylase+pregnenolone, 
                                 values=c(0.99, 0.01, 0.01, 0.99,
                                          0.01, 0.99, 0.01, 0.99), 
                                 levels=yn)
  dhea <- cptable(~ dhea | cyp17a1_lyase+hydroxypregnenolone,
                  values = c(0.99, 0.01, 0.01, 0.99,
                             0.01, 0.99, 0.01, 0.99),
                  levels=yn)
  progesterone <- cptable(~ progesterone | hsd3b1+pregnenolone,
                          values = c(0.99, 0.01, 0.01, 0.99,
                                     0.01, 0.99, 0.01, 0.99),
                          levels=yn)
  hydroxyprogesterone <- cptable(~hydroxyprogesterone | cyp17a1_hydroxylase+progesterone,
                                 values = c(0.99, 0.01, 0.01, 0.99,
                                            0.01, 0.99, 0.01, 0.99),
                                 levels=yn)
  # androsteinedione <- cptable(~ androsteinedione | hsd3b1+dhea+cyp17a1_lyase+hydroxyprogesterone,
  #                             values = c(0.99, 0.50, 0.50, 0.50, 0.50, 0.01, 0.01, 0.01, 0.50, 0.01, 0.01, 0.01, 0.50, 0.01, 0.01, 0.01,
  #                                        0.01, 0.50, 0.50, 0.50, 0.50, 0.99, 0.99, 0.99, 0.50, 0.99, 0.99, 0.99, 0.50, 0.99, 0.99, 0.99),
  #                             levels=yn)
  androsteinedione1 <- cptable(~ androsteinedione1 | hsd3b1+dhea,
                               values = c(0.99, 0.01, 0.01, 0.99,
                                          0.01, 0.99, 0.01, 0.99),
                               levels=yn)
  androsteinedione2 <- cptable(~ androsteinedione2 | cyp17a1_lyase+hydroxyprogesterone,
                               values = c(0.99, 0.01, 0.01, 0.99,
                                          0.01, 0.99, 0.01, 0.99),
                               levels=yn)
  androsteinedione <- cptable(~ androsteinedione | androsteinedione1+androsteinedione2,
                              values = c(0.99, 0.01, 0.99, 0.01,
                                         0.99, 0.01, 0.01, 0.99),
                              levels=yn)
  estrone <- cptable(~ estrone | cyp19a1+androsteinedione,
                     values = c(0.99, 0.01, 0.01, 0.99,
                                0.01, 0.99, 0.01, 0.99),
                     levels=yn)
  testosterone <- cptable(~ testosterone | hsd17b3+androsteinedione,
                          values=c(0.99, 0.01, 0.01, 0.99,
                                   0.01, 0.99, 0.01, 0.99),
                          levels=yn)
  estradiol <- cptable(~ estradiol | cyp19a1+testosterone,
                       values=c(0.99, 0.01, 0.01, 0.99,
                                0.01, 0.99, 0.01, 0.99),
                       levels=yn)
  deoxycortisone <- cptable(~ deoxycortisone | cyp21a1+progesterone,
                            values=c(0.99, 0.01, 0.01, 0.99,
                                     0.01, 0.99, 0.01, 0.99),
                            levels=yn)
  corticosterone <- cptable(~ corticosterone | cyp11b1+deoxycortisone,
                            values=c(0.99, 0.01, 0.01, 0.99,
                                     0.01, 0.99, 0.01, 0.99),
                            levels=yn)
  deoxycortisol <- cptable(~ deoxycortisol | cyp21a1+hydroxyprogesterone,
                           values=c(0.99, 0.01, 0.01, 0.99,
                                    0.01, 0.99, 0.01, 0.99),
                           levels=yn)
  cortisol <- cptable(~ cortisol | cyp11b1+deoxycortisol,
                      values=c(0.99, 0.01, 0.01, 0.99,
                               0.01, 0.99, 0.01, 0.99),
                      levels=yn)
  
  
  plist <- compileCPT(list(cholesterol, cyp11a1, pregnenolone, cyp17a1_hydroxylase, cyp17a1_lyase, cyp21a1, cyp11b1, hsd3b1, hsd17b3, cyp19a1, hydroxypregnenolone,
                           dhea, progesterone, hydroxyprogesterone, androsteinedione1, androsteinedione2, androsteinedione, estrone, testosterone, estradiol,
                           deoxycortisone, corticosterone, deoxycortisol, cortisol))
  bn_test <- grain(plist)
  bn_query <- querygrain(setEvidence(bn_test, evidence=evidence_list))
  bn_results <- format_results(bn_query)
  return(bn_results)
}

convert_data <- function(x){
  return(switch(x$aenm,
         CEETOX_H295R_11DCORT_dn = c("deoxycortisol", "no"),
         CEETOX_H295R_11DCORT_up = c("deoxycortisol", "yes"),
         CEETOX_H295R_OHPREG_dn = c("hydroxypregnenolone", "no"),
         CEETOX_H295R_OHPREG_up = c("hydroxypregnenolone", "yes"),
         CEETOX_H295R_PROG_dn = c("progesterone", "no"),
         CEETOX_H295R_PROG_up = c("prosterone", "yes"),
         CEETOX_H295R_OHPROG_dn = c("hydroxyprogesterone", "no"),
         CEETOX_H295R_OHPROG_up = c("hydroxyprogesterone", "yes"),
         CEETOX_H295R_ANDR_dn = c("androsteinedione", "no"),
         CEETOX_H295R_ANDR_up = c("androsteinedione", "yes"),
         CEETOX_H295R_TESTO_dn = c("testosterone", "no"),
         CEETOX_H295R_TESTO_up = c("testosterone", "yes"),
         CEETOX_H295R_CORTISOL_dn = c("cortisol", "no"),
         CEETOX_H295R_CORTISOL_up = c("cortisol", "yes"),
         CEETOX_H295R_DOC_dn = c("deoxycortisone", "no"),
         CEETOX_H295R_DOC_up = c("deoxycortisone", "yes"),
         CEETOX_H295R_ESTRADIOL_dn = c("estradiol", "no"),
         CEETOX_H295R_ESTRADIOL_up = c("estradiol", "yes"),
         CEETOX_H295R_ESTRONE_dn = c("estrone", "no"),
         CEETOX_H295R_ESTRONE_up = c("estrone", "yes")
         )
  )
}

mtc <- fread("suppl_data/toxsci-15-0570-File009.csv")
length(unique(mtc$chnm))
length(unique(mtc[hitc == 1, chnm]))

hits <- mtc[hitc == 1]
hits2 <- na.omit(hits)

extra_data <- NULL
extra_data$measure <- NULL
extra_data$yes_no <- NULL

for(i in 1:nrow(hits2)){
  hits_conv_info <- convert_data(hits2[i,])
  extra_data$measure <- c(extra_data$measure, hits_conv_info[1])
  extra_data$yes_no <- c(extra_data$yes_no, hits_conv_info[2])
}

hits3 <- cbind(hits2[, 4], extra_data$measure, extra_data$yes_no)
colnames(hits3) <- c("chemical_name", "measure", "yes_no")
  
# toxcast_data_split <- split(toxcast_data, f=toxcast_data$chemical_name)
toxcast_data_split <- split(hits3, f=hits3$chemical_name)
toxcast_steroidogenesis_results <- lapply(toxcast_data_split, steroidogenesis_bn)

length(toxcast_steroidogenesis_results)

save(toxcast_steroidogenesis_results, file="toxcast_steroidogenesis_results.RData")





