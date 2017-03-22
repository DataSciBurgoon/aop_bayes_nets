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
hydroxypregnenolone <- cptable(~ hydroxypregnenolone | cyp17a1_hydroxylase+pregnenolone, 
                               values=c(0.99, 0.01, 0.01, 0.01, 
                                        0.01, 0.99, 0.99, 0.99), 
                               levels=yn)
dhea <- cptable(~ dhea | cyp17a1_lyase+hydroxypregnenolone,
                values = c(0.99, 0.01, 0.01, 0.01, 
                           0.01, 0.99, 0.99, 0.99),
                levels=yn)
progesterone <- cptable(~ progesterone | hsd3b1+pregnenolone,
                        values = c(0.99, 0.01, 0.01, 0.01, 
                                   0.01, 0.99, 0.99, 0.99),
                        levels=yn)
hydroxyprogesterone <- cptable(~hydroxyprogesterone | cyp17a1_hydroxylase+progesterone,
                               values = c(0.99, 0.01, 0.01, 0.01, 
                                          0.01, 0.99, 0.99, 0.99),
                               levels=yn)
# androsteinedione <- cptable(~ androsteinedione | hsd3b1+dhea+cyp17a1_lyase+hydroxyprogesterone,
#                             values = c(0.99, 0.50, 0.50, 0.50, 0.50, 0.01, 0.01, 0.01, 0.50, 0.01, 0.01, 0.01, 0.50, 0.01, 0.01, 0.01,
#                                        0.01, 0.50, 0.50, 0.50, 0.50, 0.99, 0.99, 0.99, 0.50, 0.99, 0.99, 0.99, 0.50, 0.99, 0.99, 0.99),
#                             levels=yn)
androsteinedione <- cptable(~ androsteinedione | hsd3b1+dhea+cyp17a1_lyase+hydroxyprogesterone,
                            values = c(0.99, 0.75, 0.75, 0.75, 0.75, 0.01, 0.01, 0.01, 0.75, 0.01, 0.01, 0.01, 0.75, 0.01, 0.01, 0.01,
                                       0.01, 0.25, 0.25, 0.25, 0.25, 0.99, 0.99, 0.99, 0.25, 0.99, 0.99, 0.99, 0.25, 0.99, 0.99, 0.99),
                            levels=yn)
estrone <- cptable(~ estrone | cyp19a1+androsteinedione,
                   values = c(0.99, 0.01, 0.01, 0.01, 
                              0.01, 0.99, 0.99, 0.99),
                   levels=yn)
testosterone <- cptable(~ testosterone | hsd17b3+androsteinedione,
                        values=c(0.99, 0.01, 0.01, 0.01, 
                                 0.01, 0.99, 0.99, 0.99),
                        levels=yn)
estradiol <- cptable(~ estradiol | cyp19a1+testosterone,
                     values=c(0.99, 0.01, 0.01, 0.01, 
                              0.01, 0.99, 0.99, 0.99),
                     levels=yn)


plist <- compileCPT(list(cholesterol, cyp11a1, pregnenolone, cyp17a1_hydroxylase, cyp17a1_lyase, hsd3b1, hsd17b3, cyp19a1, hydroxypregnenolone,
                         dhea, progesterone, hydroxyprogesterone, androsteinedione, estrone, testosterone, estradiol))
bn_test <- grain(plist)
png("steroidogenesis_network.png", width=5000, height=5000, res=800, pointsize=20)
plot(bn_test)
dev.off()
querygrain(setEvidence(bn_test, evidence=list(cholesterol="yes", pregnenolone="yes", dhea="no",
                                              androsteinedione="no", hydroxypregnenolone="no",
                                              hydroxyprogesterone="no", testosterone="no")))

querygrain(setEvidence(bn_test, evidence=list(cholesterol="yes", pregnenolone="yes")))
