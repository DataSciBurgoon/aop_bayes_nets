library(gRain)

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
#Ketoconazole
querygrain(setEvidence(bn_test, evidence=list(hydroxyprogesterone="no",
                                              androsteinedione="no", testosterone="no",
                                              estrone="no", deoxycortisone="yes",
                                              cortisol="no", deoxycortisol="no")))

#Imazalil
querygrain(setEvidence(bn_test, evidence=list(progesterone="yes", hydroxyprogesterone="yes",
                                              cortisol="no", deoxycortisol="no", 
                                              androsteinedione="no", testosterone="no",
                                              estrone="no")))

