#1150 LTP ophthalmologist; 174 received 1 or more LTPs on the same eye during follow-up
#234 LTP optometrist; 84 received 1 or more LTPs on the same eye during follow-up

library(rstan)

opt_eyes <- 234
opt_return <- 84
oph_eyes <- 1150
oph_return <- 174

# THE MODEL.
eye_modelString = "
  data {
    int<lower=0> N; //number of items
    int y[N]; // y is an N-length vector of ints
  } parameters {
    real <lower=0, upper=1> theta;
  } model {
    theta ~ beta(1,1);
    y ~ bernoulli(theta);
  }
"

stanDso <- stan_model( model_code=eye_modelString )
N <- opt_eyes
z <- opt_return
y <- c( rep( 1, z), rep( 0, N-z))
dataList <- list( y = y , N = N )
#This is going to do the MCMC for Tokar's data. Note we're using a flat prior.
stanFit <- sampling( object = stanDso , data = dataList , chains = 3, iter = 5000 , warmup = 200 , thin = 1 )
stan_hist(stanFit)

N <- oph_eyes
z <- oph_return
y <- c( rep( 1, z), rep( 0, N-z))
dataList <- list( y = y , N = N )
#This is going to do the MCMC for Tokar's data. Note we're using a flat prior.
stanFitOph <- sampling( object = stanDso , data = dataList , chains = 3, iter = 5000 , warmup = 200 , thin = 1 )
stan_hist(stanFitOph)


opt_posterior <- as.data.frame(extract(stanFit)[[1]])
oph_posterior <- as.data.frame(extract(stanFitOph)[[1]])
colnames(opt_posterior) <- "theta"
colnames(oph_posterior) <- "theta"
combined_data <- rbind(opt_posterior, oph_posterior)
code <- c(rep("optometrist", nrow(opt_posterior)), rep("ophthamologist", nrow(oph_posterior)))
combined_data <- cbind(combined_data, code=code)
ggplot(combined_data, aes(theta, fill = code)) + geom_density(alpha = 0.2)


posterior_diff <- as.data.frame(opt_posterior - oph_posterior)
colnames(posterior_diff) <- "difference"
ggplot(posterior_diff, aes(difference, fill=1)) +
  geom_density(alpha=1)

quantile(posterior_diff$difference, probs=c(0.05, 1))
