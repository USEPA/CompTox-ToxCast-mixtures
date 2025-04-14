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
curve_bootstrap <- function(jj, curve.info, boot_num=1000, type="single",tc=FALSE){
  if (curve.info$hitc[jj] >= 0.9){
    if (type=="single"){
      m4id <- curve.info$m4id[jj]
    } else if (type == "mix"){
      m4id <- curve.info$m4id_mix[jj]
    }
    if (tc==FALSE){
      comp_dat <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==m4id)])
    } else if (tc==TRUE){
      # TC chemical library
      comp_dat <- mc3_tcsubset %>% filter(m3id %in% mc4_agg_tcsubset$m3id[which(mc4_agg_tcsubset$m4id==m4id)])
    }
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
bayes_curve <- function(ii, curve.info, n_iter=10000, n_burn=2000, type="single", poly2_correct=FALSE, tc=FALSE){
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
        if (tc){
          comp_dat <- mc3_tcsubset %>% filter(m3id %in% mc4_agg_tcsubset$m3id[which(mc4_agg_tcsubset$m4id==m4id_curve)])
        } else if (!tc){
          comp_dat <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==m4id_curve)])
        }
      } else{
        if (tc){
          single_dat <- mc5_attg_endpts[m4id==m4id_curve,]
          comp_dat <- mc3_tcsubset %>% filter(m3id %in% mc4_agg_tcsubset$m3id[which(mc4_agg_tcsubset$m4id==m4id_curve)])
        } else if (!tc){
          single_dat <- mc5[m4id==m4id_curve,]
          comp_dat <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==m4id_curve)])
        }
      }
      # medians
      rmds <- tapply(comp_dat$resp, comp_dat$conc, median)
      # For each curve type: (1) define bayesian likelihood model and (2) define prior as uniform with same bounds as ToxCast
      if (single_dat$modl=="poly1"){
        # define the model
        a0 <- max(abs(rmds))/max(comp_dat$conc)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat), a_low=-1e8*a0, a_upp=1e8*a0)
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
                                       a ~ dunif(a_low, a_upp)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "sigma")
      } else if (single_dat$modl=="poly2"){
        a0 <- max(abs(rmds))
        b0 <- max(comp_dat$conc)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat), b1_low=-1e8*a0,
                     b1_upp=1e8*a0, b2_low=-1e8*b0, b2_upp=1e8*b0)
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
                                       b1 ~ dunif(b1_low, b1_upp)
                                       b2 ~ dunif(b2_low, b2_upp)
                                       }")
        # define the initial values
        inits <- list(b1=single_dat$a/single_dat$b, b2=single_dat$a/(single_dat$b^2), tau=(1/exp(single_dat$er))^2)
        params <- c("b1", "b2", "sigma")
      } else if (single_dat$modl=="pow"){
        a0 <- max(abs(rmds))
        b0 <- max(comp_dat$conc)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat),
                     a_low=-1e8*a0, a_upp=1e8*a0)
        # p bound set 0.3-20 in ToxCast
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
                                       a ~ dunif(a_low, a_upp)
                                       p ~ dunif(0.3, 20)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "p", "sigma")
      } else if (single_dat$modl=="hill"){
        tp_val <- 1.2*max(abs(min(comp_dat$resp)),abs(max(comp_dat$resp)))
        ga_low <- min(comp_dat$conc)/10
        ga_upp <- max(comp_dat$conc)*(10^0.5)
        # p bounds set in TOxCast of 0.3 - 8
        # for hill, concentrations are internally converted to log in ToxCast
        # try keeping it not in log scale because I am lazy
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat),
                     tp_val=tp_val,ga_low=ga_low, ga_upp=ga_upp)
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
                                       ga ~ dunif(ga_low, ga_upp)
                                       tp ~ dunif(-tp_val, tp_val)
                                       p ~ dunif(0.3, 8)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "p", "sigma")
      } else if (single_dat$modl=="gnls"){
        tp_val <- 1.2*max(abs(min(comp_dat$resp)),abs(max(comp_dat$resp)))
        ga_low <- min(comp_dat$conc)/10
        #ga_upp <- max(comp_dat$conc)*(10^0.5)
        ga_upp <- max(comp_dat$conc) # change from ToxCast bounds
        #la_low <- min(comp_dat$conc)/10
        la_low <- min(comp_dat$conc) # change from ToxCast bounds
        #la_upp <- max(comp_dat$conc)*(10^2)
        la_upp <- max(comp_dat$conc)*(50) # change from ToxCast bounds
        # p and q bounds set in TOxCast of 0.3 - 8
        # for gnls, concentrations are internally converted to log in ToxCast
        # try keeping it not in log scale because I am lazy
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat),
                     tp_val=tp_val,ga_low=ga_low, ga_upp=ga_upp, la_low=la_low, la_upp=la_upp)
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
                                       la ~ dunif(la_low, la_upp)
                                       tp ~ dunif(-tp_val, tp_val)
                                       p ~ dunif(0.3, 8)
                                       q ~ dunif(0.3, 8)
                                       ga ~ dunif(ga_low, ga_upp)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, p=single_dat$p, q=single_dat$q, la=single_dat$la, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "p", "q", "la", "sigma")
      } else if (single_dat$modl=="exp2"){
        a0 <- max(abs(rmds))
        b0 <- max(comp_dat$conc)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat), a_low=-1e8*a0,
                     a_upp=1e8*a0, b_low=1e-2*b0, b_upp=1e8*b0)
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
                                       a ~ dunif(a_low, a_upp)
                                       b ~ dunif(b_low, b_upp)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, b=single_dat$b, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "b", "sigma")
      } else if (single_dat$modl=="exp3"){
        a0 <- max(abs(rmds))
        b0 <- max(comp_dat$conc)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat), a_low=-1e8*a0,
                     a_upp=1e8*a0, b_low=1e-2*b0, b_upp=1e8*b0)
        # p bounds are set to 0.3-8 in toxcast
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
                                       a ~ dunif(a_low,a_upp)
                                       b ~ dunif(b_low,b_upp)
                                       p ~ dunif(0.3, 8)
                                       }")
        # define the initial values
        inits <- list(a=single_dat$a, b=single_dat$b, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("a", "b", "p", "sigma")
      } else if (single_dat$modl=="exp4"){
        a0 <- max(abs(rmds))
        ga_low <- min(comp_dat$conc)/10
        ga_upp <- max(comp_dat$conc)*(10^0.5)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat), tp_low= -1.2*a0,
                     tp_upp=1.2*a0, ga_low=ga_low, ga_upp=ga_upp)
        # p bounds are set to 0.3-8 in toxcast
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
                                       ga ~ dunif(ga_low, ga_upp)
                                       tp ~  dunif(tp_low, tp_upp)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "sigma")
      } else if (single_dat$modl=="exp5"){
        a0 <- max(abs(rmds))
        ga_low <- min(comp_dat$conc)/10
        ga_upp <- max(comp_dat$conc)*(10^0.5)
        data <- list(conc=comp_dat$conc, resp=comp_dat$resp, n=nrow(comp_dat), tp_low= -1.2*a0,
                     tp_upp=1.2*a0, ga_low=ga_low, ga_upp=ga_upp)
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
                                       ga ~ dunif(ga_low, ga_upp)
                                       tp ~  dunif(tp_low, tp_upp)
                                       p ~ dunif(0.3, 8)
                                       }")
        # define the initial values
        inits <- list(ga=single_dat$ga, tp=single_dat$tp, p=single_dat$p, tau=(1/exp(single_dat$er))^2)
        params <- c("ga", "tp", "p", "sigma")
      }
      
      # define model
      model <- jags.model(model_string, data=data, inits=inits, n.chains=1, quiet=TRUE)
      # using one chain, but can use more to evaluate convergence
      
      # run burn-in
      update(model, n_burn, progress.bar="none")
      
      # run samples
      samples <- coda.samples(model, variable.names=params, n.iter=n_iter, progress.bar="none")
      
      #summary(samples)
      
      #plot(samples)
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
