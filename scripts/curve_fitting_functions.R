# Functions to perform bootstrap resampling and Bayesian concentration-response curve fitting methods
library(rjags)

#' Bootstrap resampling of concentration-response curves
#' outputs curve shape parameters to save bootstrap resampled curves
#'
#' @param jj integer | index of curve in dataframe
#' @param curve.info data.frame | dataframe to store single or mixture curve information
#' @param boot_num integer | number of bootstrap samples to run
#' @param type If "single", pulls curve information from single component data frame. If "mix", pulls curve information from mixture data frame
#' @param tc If tc = FALSE, data is pulled from the test dataset. If tc = TRUE, data is pulled from the legacy ToxCast chemical library
#' 
#' @return dataframe of bootstraped resampled curve parameters for all parameters of the input curve shape type.
#'  Returns NA if curve is inactive
curve_bootstrap <- function(jj, curve.info, boot_num=1000, type="single"){
  if (curve.info$hitc[jj] >= 0.9){
    if (type=="single"){
      m4id <- curve.info$m4id[jj]
    } else if (type == "mix"){
      m4id <- curve.info$m4id_mix[jj]
    }
    comp_dat <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==m4id)])
    bmad <- curve.info$bmad[jj]
    datchemval <- data.table(conc=comp_dat$conc, resp=comp_dat$resp, bmad=bmad)
    datchemstat = data.table(conc = sort(unique(datchemval$conc)))
    tempmat <- c()
    
    for (conc_val in datchemstat[, conc]){  
      number_measurements <- length(which(datchemval$conc==conc_val))
      resp_vect <- datchemval[conc == conc_val, resp] #gets the measured responses at this concentration
      
      if(number_measurements > 1){
        for (numpoint in 1:number_measurements){
          vect <- sample(resp_vect, boot_num, replace = TRUE) +
            rnorm(n = boot_num, mean = 0, sd = bmad)
          #rst(n=replicates, mu=0, sigma=bmad, nu=4)
          tempmat <- rbind(tempmat, vect)
        }
      }
      else{
        vect <- resp_vect + #This is necessary when only one dose, as sample will do 1:X if length of X == 1
          rnorm(n = boot_num, mean = 0, sd = bmad)
        #rst(n=replicates, mu=0, sigma=bmad, nu=4)
        tempmat <- rbind(tempmat, vect)
        
      }
    }
    datchemsample <- data.table(tempmat)
    # Fit each bootstrap sample using tcplFit
    fname <- paste0("fit", curve.info$modl[jj]) # requires each model function have name "fit____" where ____ is the model name
    # use do.call to call fit function; cnst has different inputs than others.
    if(fname != "fitpoly2"){
      temp2 <- sapply(datchemsample,
                      match.fun(fname),
                      conc = sort(datchemval$conc),
                      bidirectional = TRUE,
                      verbose = FALSE,
                      nofit = FALSE,
                      errfun = "dt4")
    }else{
      temp2 <- sapply(datchemsample,
                      match.fun(fname),
                      conc = sort(datchemval$conc),
                      bidirectional = TRUE,
                      verbose = FALSE,
                      nofit = FALSE,
                      errfun = "dt4",
                      biphasic = TRUE)
    }
    temp3 <- as.data.frame(t(temp2)[,temp2["pars",][[1]]])
    #unlist columns
    setDT(temp3)
    boot_params <- temp3[, lapply(.SD, unlist)]
    # save bootstrapped parameters, remove NA
    return(na.omit(boot_params))
  } else {
    return(NA)
  }
}

#' Bayesian method for fitting concentration-response curves
#' outputs posterior parameter distributions for the parameters of a given curve shape type
#'
#' @param jj integer | index of curve in dataframe
#' @param curve.info data.frame | dataframe to store single or mixture curve information
#' @param n_iter integer | number MCMC iterations to run
#' @param n_burn integer | number of MCMC initial iterations to run and delete as burnin
#' @param type If "single", pulls curve information from single component data frame. If "mix", pulls curve information from mixture data frame
#' @param poly2_correct If poly2_correct=TRUE, poly2 correction for concave up curve is applied
#' @param tc If tc = FALSE, data is pulled from the test dataset. If tc = TRUE, data is pulled from the legacy ToxCast chemical library
#'  
#' @return dataframe of bayesian posterior parameter distribution samples for all parameters of the input curve shape type
#'  Returns NA if input curve is inactive.
bayes_curve <- function(ii, curve.info, n_iter=10000, n_burn=10000, type="single", poly2_correct=FALSE){
    if (curve.info$hitc[ii]<0.9){
      bayesparams <- NA
    } else {
      if (type=="single"){
        m4id_curve <- curve.info$m4id[ii]
      } else if (type=="mix"){
        m4id_curve <- curve.info$m4id_mix[ii]
      }
      # # define the data
      if (poly2_correct){
        single_dat <- curve.info
        comp_dat <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==m4id_curve)])
      } else{
        single_dat <- mc5[m4id==m4id_curve,]
        comp_dat <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==m4id_curve)])
      }
      # medians
      rmds <- tapply(comp_dat$resp, comp_dat$conc, median)
      # For each curve type: (1) define bayesian likelihood model and (2) define prior as uniform with same bounds as ToxCast
      if (single_dat$modl=="poly1"){
        # define the model
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- a*conc[i]
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       a ~ dnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "sigma")
      } else if (single_dat$modl=="poly2"){
        # data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        # model_string <- textConnection("model{
        #                                # Likelihood
        #                                for(i in 1:n){
        #                                 resp[i] ~ dt(mu[i],(1/(sigma^2)),k)
        #                                 mu[i] <- b1*conc[i] + b2*(conc[i]*conc[i])
        #                                }
        #                                # Priors
        #                                k <- 4
        #                                sigma ~ dt(0,(2.5^(-2)),1)
        #                                b1 ~ dnorm(0,1)
        #                                b2 ~ dnorm(0,1)
        #                                }")
        # # define the initial values
        # inits <- list(b1=single_dat$a/single_dat$b, b2=single_dat$a/(single_dat$b^2), sigma=exp(single_dat$er))
        # params <- c("b1", "b2", "sigma")
        
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- b1*conc[i] + b2*(conc[i]*conc[i])
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       b1 ~ dnorm(0,1)
                                       b2 ~ dnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(b1=single_dat$a/single_dat$b, b2=single_dat$a/(single_dat$b^2), tau=(1/exp(single_dat$er))^2)
        params <- c("b1", "b2", "sigma")
      } else if (single_dat$modl=="pow"){
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- a*(conc[i]^p)
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       a ~ dlnorm(0,1)
                                       p ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "p", "sigma")
      } else if (single_dat$modl=="hill"){
        # for hill, concentrations are internally converted to log in ToxCast
        # here, keep it not in log scale
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- tp/(1+(ga/conc[i])^p)
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       ga ~ dlnorm(1,1)
                                       tp ~ dlnorm(0,1)
                                       p ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "p", "sigma")
      } else if (single_dat$modl=="gnls"){
        # # for gnls, concentrations are internally converted to log in ToxCast
        # # here, we just keep it in non-log space for simplicity
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- tp/((1+(ga/conc[i])^p)*(1+(conc[i]/la)^q))
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       ga ~ dlnorm(1,1)
                                       la ~ dlnorm(3,1)
                                       tp ~ dlnorm(0,1)
                                       p ~ dlnorm(0,1)
                                       q ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, p=single_dat$p, q=single_dat$q, la=single_dat$la, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "p", "q", "la", "sigma")
      } else if (single_dat$modl=="exp2"){
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- a*(exp(conc[i]/b)-1)
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       a ~ dlnorm(0,1)
                                       b ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, b=single_dat$b, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "b", "sigma")
      } else if (single_dat$modl=="exp3"){
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- a*(exp((conc[i]/b)^p)-1)
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       a ~ dlnorm(0,1)
                                       b ~ dlnorm(0,1)
                                       p ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, b=single_dat$b, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "b", "p", "sigma")
      } else if (single_dat$modl=="exp4"){
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- tp*(1-2^(-conc[i]/ga))
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       ga ~ dlnorm(1,1)
                                       tp ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "sigma")
      } else if (single_dat$modl=="exp5"){
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat))        
        model_string <- textConnection("model{
                                       # Likelihood
                                       for(i in 1:n){
                                        resp[i] ~ dt(mu[i],tau,k)
                                        mu[i] <- tp*(1-2^(-(conc[i]/ga)^p))
                                       }
                                       # Priors
                                       k <- 4
                                       tau ~ dgamma(0.1,0.1)
                                       sigma <- 1/sqrt(tau)
                                       ga ~ dlnorm(1,1)
                                       tp ~ dlnorm(0,1)
                                       p ~ dlnorm(0,1)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "p", "sigma")
      }
      
      # define model
      model <- jags.model(model_string, data=data, inits=inits, n.chains=1, quiet=TRUE)
      # model <- jags.model(model_string, data=data, inits=inits, n.chains=3, quiet=TRUE)
      # using one chain, but can use more to evaluate convergence
      
      # run burn-in
      update(model, n_burn, progress.bar="none")
      
      # run samples
      samples <- coda.samples(model, variable.names=params, n.iter=n_iter, progress.bar="none")
      
      # summary(samples)
      # plot(samples, trace=TRUE,density=FALSE)
      # gelman.diag(samples)
      
      # convert poly2 back to a and b parameters
      if (single_dat$modl=="poly2"){
        poly2_temp <- data.frame(samples[[1]])
        bayesparams <- data.frame(a=(poly2_temp$b1*poly2_temp$b1/poly2_temp$b2),
                                  b=(poly2_temp$b1/poly2_temp$b2),
                                  sigma=poly2_temp$sigma)
      } else{
        # convert chain to data frame
        bayesparams <- data.frame(samples[[1]])     
      }
    }
  
  return(bayesparams)
}
