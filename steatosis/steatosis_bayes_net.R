library(gRain)

#Building the Bayes Net
yn <- c("yes", "no")
nrf2 <- cptable(~ nrf2, values=c(0.5, 0.5), levels=yn)
fxr <- cptable(~ fxr | nrf2, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
shp <- cptable(~ shp | fxr, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
lxr <- cptable(~ lxr | shp, values=c(0.05, 0.95, 0.95, 0.05), levels=yn)
ppar_alpha <- cptable(~ ppara | fxr+shp+lxr, 
                      values=c(0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.01, 0.99,
                               0.99, 0.01, 0.50, 0.50, 0.50, 0.50, 0.01, 0.99),
                      levels=yn)
hsd17b4 <- cptable(~ hsd17b4 | ppara, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
fatty_acid_beta_oxidation <- cptable(~ fatty_acid_beta_oxidation | hsd17b4,
                                     values = c(0.99, 0.01, 0.01, 0.99), levels=yn)
steatosis <-  cptable(~steatosis | cytosolic_fatty_acids, 
                      values=c(0.99, 0.01, 0.01, 0.99), levels=yn)
lrh1 <- cptable(~lrh1 | shp,
                values=c(0.05, 0.95, 0.95, 0.05), levels=yn)
pi3k <- cptable(~pi3k, values=c(0.50, 0.50), levels=yn)
mtorc2 <- cptable(~mtorc2 | ir, values=c(0.99, 0.01, 0.01, 0.99), levels=yn)
ir <- cptable(~ir, values=c(0.50, 0.50), levels=yn)
akt <- cptable(~akt | pi3k+mtorc2, values=c(0.95, 0.05, 0.05, 0.95,
                                            0.95, 0.05, 0.05, 0.95),
               levels=yn)
lfabp <- cptable(~lfabp | akt+pi3k, values=c(0.95, 0.05, 0.05, 0.05, 0.50, 0.50, 0.05, 0.95),
                 levels=yn)
pparg <- cptable(~pparg | lfabp, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
fas <- cptable(~fas | lrh1+lxr+pparg, values=c(0.95, 0.05, 0.75, 0.25, 0.75, 0.25, 0.50, 0.50,
                                               0.75, 0.25, 0.50, 0.50, 0.50, 0.50, 0.01, 0.99),
               levels=yn)
mtorc1 <- cptable(~mtorc1 | akt, values=c(0.95, 0.05, 0.05, 0.95), levels=yn)
apkc <- cptable(~apkc | pi3k, values=c(0.99, 0.01, 0.01, 0.99), levels=yn)
srebp1 <- cptable(~srebp1 | mtorc1 + apkc, values=c(0.95, 0.05, 0.05, 0.95,
                                                    0.95, 0.05, 0.05, 0.95),
                  levels=yn)
scd1 <- cptable(~scd1 | srebp1, values=c(1, 0, 0, 1), levels=yn)
lipogenesis <- cptable(~lipogenesis | scd1 + fas, values=c(0.95, 0.05, 0.05, 0.95,
                                                           0.95, 0.05, 0.05, 0.95),
                       levels=yn)
cytosolic_fatty_acids <- cptable(~cytosolic_fatty_acids | lipogenesis + lfabp + fatty_acid_beta_oxidation, 
                                 values=c(0.01, 0.99, 0.01, 0.99, 0.01, 0.99, 0.01, 0.99,
                                          0.99, 0.01, 0.66, 0.34, 0.66, 0.34, 0.50, 0.50),
                                 levels=yn)

plist <- compileCPT(list(nrf2, fxr, shp, lxr, ppar_alpha, hsd17b4, 
                         fatty_acid_beta_oxidation, steatosis, lrh1, pi3k, 
                         mtorc2, ir, akt, lfabp, pparg, fas, mtorc1, apkc, srebp1, scd1, lipogenesis,
                         cytosolic_fatty_acids))
bn_test <- grain(plist)
plot(bn_test)
querygrain(setEvidence(bn_test, evidence=list(hsd17b4="no")))
querygrain(setEvidence(bn_test, evidence=list(nrf2="no")))
querygrain(setEvidence(bn_test, evidence=list(nrf2="yes")))
querygrain(setEvidence(bn_test, evidence=list(hsd17b4="yes")))
querygrain(setEvidence(bn_test, evidence=list(lipogenesis="yes", hsd17b4="no")))
querygrain(setEvidence(bn_test, evidence=list(lipogenesis="yes", hsd17b4="yes")))
querygrain(setEvidence(bn_test, evidence=list(lipogenesis="yes", nrf2="no", akt="yes")))
