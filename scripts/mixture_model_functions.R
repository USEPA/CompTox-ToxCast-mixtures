library(rootSolve)

# Functions to compute mixture models

#' Compute concentration-response curve from selected ToxCast model type
#'
#' @param ss integer | index of sample number, used if model is evaluated with sampled parameters/curves
#' @param dat data.frame | data frame containing mc5 level ToxCast information for curve (modl,conc_max,conc_min)
#' @param params data.frame | dataframe containing parameter values for curve type
#' @param XX vector | concentration values for which to calculate curve
#' @param sample logical | if sample=TRUE, parameters are provided as a sampled sets from bootsrap of Bayesian fitting.
#' If sample=FALSE, only one set of curve parameters are provided for each component
#' 
#' @return vector | response values calculated for curve across input concentration values for given curve type
toxcast_model <- function(ss=NA, dat, params, XX, sample=FALSE){
  if (sample){
    params <- params[ss,]
  }
  if (dat$hitc < 0.9){
    mu <- 0
  } else if (dat$modl=="cnst"){
    mu <- tcplfit2::cnst(c(params$er),XX)
  } else if (dat$modl=="poly1"){
    mu <- tcplfit2::poly1(c(params$a,params$er),XX)
  } else if (dat$modl=="poly2"){
    mu <- tcplfit2::poly2(c(params$a,params$b,params$er),XX)
  } 
  else if (dat$modl=="pow"){
    mu <- tcplfit2::pow(c(params$a,params$p,params$er),XX)
  } else if (dat$modl=="hill"){
    mu <- tcplfit2::hillfn(c(params$tp,params$ga,params$p,params$er),XX)
  } else if (dat$modl=="gnls"){
    mu <- tcplfit2::gnls(c(params$tp,params$ga,params$p,params$la,params$q,params$er),XX)
  } else if (dat$modl=="exp2"){
    mu <- tcplfit2::exp2(c(params$a,params$b,params$er),XX)
  } else if (dat$modl=="exp3"){
    mu <- tcplfit2::exp3(c(params$a,params$b,params$p,params$er),XX)
  } else if (dat$modl=="exp4"){
    mu <- tcplfit2::exp4(c(params$tp,params$ga,params$er),XX)
  } else if (dat$modl=="exp5"){
    mu <- tcplfit2::exp5(c(params$tp,params$ga,params$p,params$er),XX)
  }
  return(mu)
}

#' Compute inverse of concentration response curve for a given response value
#'
#' @param y vector | response values at which to calculate the concentration for the curve
#' @param params data.frame | dataframe containing parameter values for curve type
#' @param single_dat data.frame | data frame containing mc5 level ToxCast information for curve (modl,conc_max,conc_min)
#' 
#' @return vector | concentration values of curve for given response values
single_curve_inv <- function(y, params, single_dat){
  if (single_dat$modl %in% c('poly1', 'poly2','pow', 'exp2', 'exp3')){
    ACy_calc <- sapply(y, tcplfit2::acy, list(p=params$p, a=params$a, b=params$b), type=single_dat$modl)
  } else if (single_dat$modl %in% c('hill', 'exp4', 'exp5')){
    ACy_calc <- sapply(y, tcplfit2::acy, list(tp=params$tp, ga=params$ga, p=params$p, la=params$la, q=params$q), type=single_dat$modl)
  } else if (single_dat$modl == 'gnls'){
    ACy_calc <- sapply(y,gnls_inv,params=params,single_dat=single_dat)
  }
  return(ACy_calc)
}

#' Compute gainloss inverse of concentration response curve for a given response value
#'
#' @param y vector | response value at which to calculate the concentration for the curve
#' @param params data.frame | dataframe containing parameter values for gainloss curve
#' @param single_dat data.frame | data frame containing mc5 level ToxCast information for curve (modl,conc_max,conc_min)
#' 
#' @return vector | concentration value of curve for given response value
gnls_inv <- function(y,params,single_dat){
  #### copied from tcplfit2
  toploc = try(uniroot(tcplfit2::gnlsderivobj, c(params$ga,params$la), tp=params$tp,
                       ga=params$ga, p=params$p, la=params$la, q=params$q, tol=1e-8)$root,silent=TRUE)
  #If toploc fails, set topval to tp, set toploc to NA
  if(is(toploc,"try-error")){
    # warning("toploc could not be found numerically")
    topval = params$tp
    toploc = NA_real_
  } else {
    #get actual top, as opposed to theoretical tp.
    topval = tcplfit2::gnls(c(params$tp, params$ga, params$p, params$la, params$q), toploc)
  }
  #If y >= top, don't try to solve.
  if(abs(y) > abs(topval)) {
    # warning("y is greater than gnls top in function acy, returning NA")
    ACy_calc <- NA
  }
  #solve for acy
  output = try(uniroot(tcplfit2::acgnlsobj, c(1e-8, toploc), y=y, tp=params$tp,
                       ga=params$ga, p=params$p, la=params$la, q=params$q, tol=1e-8)$root,silent=TRUE)
  if(is(output,"try-error")){
    f1 = function(x){
      (params$tp/((1+((params$ga/x)^params$p))*(1+((x/params$la)^params$q)))) - y}
    
    ans <- uniroot.all(f1,lower=single_dat$conc_min,upper=single_dat$conc_max)
    if(length(ans)==0){
      ACy_calc <- NA
    } else {
      ACy_calc <- min(ans,na.rm=TRUE)
    }
  } else {
    ACy_calc <- output
  }
}

#' Compute top of a concentration-response curve
#'
#' @param dat data.frame | data frame containing mc5 level ToxCast information for curve (modl,conc_max,conc_min)
#' @param params data.frame | dataframe containing parameter values for curve type
#' 
#' @return numeric | top response value of concentration-response curve
single_curve_top <- function(dat,params){
  mu <- toxcast_model(dat=dat, params=params, XX=seq(dat$conc_min, dat$conc_max, length.out=500))
  topval <- mu[which.max(abs(mu))]    
  return(topval)
}

#' Compute minimum of a concentration-response curve
#'
#' @param dat data.frame | data frame containing mc5 level ToxCast information for curve (modl,conc_max,conc_min)
#' @param params data.frame | dataframe containing parameter values for curve type
#' 
#' @return numeric | minimum response value of concentration-response curve
single_curve_min <- function(dat,params){
  # could change X to only cover data range, but use a small concentration to get calculatable minimum
  mu <- toxcast_model(dat=dat, params=params, XX=seq(1E-6, dat$conc_max, length.out=500))
  resp_min <- mu[which.min(abs(mu))]
  return(resp_min)
}


#' Compute Concentration Addition mixture model for one given response value
#' 
#' @param ss integer | index of sample number for bootstrap resamples / bayesian posterior parameter samples
#' @param YY numeric | given non-normalized response value to calculate mixture CA model at
#' @param comp1 data frame | component 1 concentration-response curve information from ToxCast mc level 5
#' @param comp2 data frame | component 2 concentration-response curve information from ToxCast mc level 5
#' @param params1 data frame | component 1 curve parameters from ToxCast, bootstrap, or Bayesian curve fit
#' @param params2 data frame | component 2 curve parameters from ToxCast, bootstrap, or Bayesian curve fit
#' @param cfrac1 numeric | mixture concentration fraction of component 1
#' @param cfrac2 numeric | mixture concentration fraction of component 1
#' @param condition1 logical | component 1 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param condition2 logical | component 2 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param sample logical | if sample=TRUE, parameters are provided as a sampled sets from bootsrap of Bayesian fitting. If sample=FALSE, only one set of curve parameters are provided for each component
#'
#' @return numeric vector| mixture concentration computed from CA model for one given response YY
ca_point_conc <- function(ss=NA,YY,comp1,comp2,params1,params2,cfrac1,cfrac2,
                          condition1,condition2,sample=FALSE){
  # YY is a non-normalized response-value
  # return concentration for that response
  if (sample==FALSE){
    comp1_top <- single_curve_top(dat=comp1,params=params1)
    # comp1_top <- comp1$top
    comp1_min <- single_curve_min(dat=comp1,params=params1)
    comp2_top <- single_curve_top(dat=comp2,params=params2)
    # comp2_top <- comp2$top
    comp2_min <- single_curve_min(dat=comp2,params=params2)
  } else if (sample==TRUE){
    if (condition1){
      params1 <- params1[ss,]
      comp1_top <- single_curve_top(dat=comp1,params=params1)
      comp1_min <- single_curve_min(dat=comp1,params=params1)
    } else {
      params1 <- NA
      comp1_top <- NA
      comp1_min <- NA
      
    }
    if (condition2){
      params2 <- params2[ss,]
      comp2_top <- single_curve_top(dat=comp2,params=params2)
      comp2_min <- single_curve_min(dat=comp2,params=params2)
    } else {
      params2 <- NA
      comp2_top <- NA
      comp2_min <- NA
    }
  }
  tu_fact <- 0.7
  # eliminate any active curves with negative tops
  if ((comp1_top<0 & condition1) | (comp2_top<0 & condition2)){
    ACy_mix <- rep(NA,length(YY))
  } else if (!condition1 & !condition2){
    ACy_mix <- rep(NA,length(YY))
  } else if (!condition1 & condition2){
    if (YY>comp2_top){
      # Add option for certain curve shapes to extend beyond concentration range to compute a value at the cutoff for bootstrapped/bayesian parameters
      if (comp2$modl %in% c("poly1","pow","exp2","exp3")){
        ACy_mix <- 1/(cfrac2/single_curve_inv(y=YY, single_dat=comp2, params=params2))
      } else if(comp2$modl == 'poly2') {
        if(!(params2$a<0 & params2$b<0)){
          ACy_mix <- 1/(cfrac2/single_curve_inv(y=YY, single_dat=comp2, params=params2))
        } else{
          ACy_mix <- 10000
        }
      } else {
        ACy_mix <- 10000
      }
    } else if (YY<comp2_min){
      ACy_mix <- -10000
    } else if (YY>=comp2_min & YY<=comp2_top){
      ACy_mix <- 1/(cfrac2/single_curve_inv(y=YY, single_dat=comp2, params=params2))
    }
  } else if (condition1 & !condition2){
    if (YY>comp1_top){
      # Add option for certain curve shapes to extend beyond concentration range to compute a value at the cutoff for bootstrapped/bayesian parameters
      if (comp1$modl %in% c("poly1","pow","exp2","exp3")){
        ACy_mix <- 1/(cfrac1/single_curve_inv(y=YY, single_dat=comp1, params=params1))
      } else if(comp1$modl == 'poly2') {
        if(!(params1$a<0 & params1$b<0)){
          ACy_mix <- 1/(cfrac1/single_curve_inv(y=YY, single_dat=comp1, params=params1))
        } else{
          ACy_mix <- 10000
        }
      } else {
        ACy_mix <- 10000
      }
    } else if (YY<comp1_min){
      ACy_mix <- -10000
    } else if (YY>=comp1_min & YY<=comp1_top){
      ACy_mix <- 1/(cfrac1/single_curve_inv(y=YY, single_dat=comp1, params=params1))
    }
  } else if (condition1 & condition2){
    # component 1 curve  evaluated only when valid
    # component 2 curve evaluated only when valid
    if (YY>=comp1_min & YY<=comp1_top & YY>=comp2_min & YY<=comp2_top){
      comp1_curve <- single_curve_inv(y=YY, single_dat=comp1, params=params1)
      comp2_curve <- single_curve_inv(y=YY, single_dat=comp2, params=params2)
      ACy_mix <- 1/((cfrac1/comp1_curve) + (cfrac2/comp2_curve))
    } else if (YY>comp1_top & YY<=comp2_top & YY>=comp1_min & YY>=comp2_min){
      ec1_70 <- single_curve_inv(y=tu_fact*comp1_top, single_dat=comp1, params=params1)
      ec2_70 <- single_curve_inv(y=tu_fact*comp1_top, single_dat=comp2, params=params2)
      comp2_curve <- single_curve_inv(y=YY, single_dat=comp2, params=params2)
      ACy_mix <- (comp2_curve/cfrac2)*(1-(cfrac1*(1/((cfrac1/ec1_70) + (cfrac2/ec2_70)))/ec1_70))
    } else if (YY>comp2_top & YY<=comp1_top & YY>=comp1_min & YY>=comp2_min){
      ec1_70 <- single_curve_inv(y=tu_fact*comp2_top, single_dat=comp1, params=params1)
      ec2_70 <- single_curve_inv(y=tu_fact*comp2_top, single_dat=comp2, params=params2) 
      comp1_curve <- single_curve_inv(y=YY, single_dat=comp1, params=params1)
      ACy_mix <- (comp1_curve/cfrac1)*(1-(cfrac2*(1/((cfrac1/ec1_70) + (cfrac2/ec2_70)))/ec2_70))
    } else if (YY<comp1_min | YY<comp2_min){
      ACy_mix <- -10000
    } else if (YY>comp1_top & YY>comp2_top){
      ACy_mix <- 10000
    }
  }
  
  return(ACy_mix)
}

#' Compute Concentration Addition mixture model for multiple given response values
#' 
#' @param ss integer | index of sample number for bootstrap resamples / bayesian posterior parameter samples
#' @param YY numeric vector| given non-normalized response values at which to calculate mixture CA model
#' @param comp1 data frame | component 1 concentration-response curve information from ToxCast mc level 5
#' @param comp2 data frame | component 2 concentration-response curve information from ToxCast mc level 5
#' @param params1 data frame | component 1 curve parameters from ToxCast, bootstrap, or Bayesian curve fit
#' @param params2 data frame | component 2 curve parameters from ToxCast, bootstrap, or Bayesian curve fit
#' @param cfrac1 numeric | mixture concentration fraction of component 1
#' @param cfrac2 numeric | mixture concentration fraction of component 1
#' @param condition1 logical | component 1 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param condition2 logical | component 2 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param sample logical | if sample=TRUE, parameters are provided as a sampled sets from bootsrap of Bayesian fitting. If sample=FALSE, only one set of curve parameters are provided for each component
#'
#' @return numeric vector| mixture concentrations computed from CA model for each input response YY
ca_full_conc <- function(ss=NA,YY,comp1,comp2,params1,params2,cfrac1,cfrac2,
                         condition1,condition2,sample=FALSE){
  # return concentration for that response
  if (!sample){
    comp1_top <- single_curve_top(dat=comp1,params=params1)
    comp1_min <- single_curve_min(dat=comp1,params=params1)
    comp2_top <- single_curve_top(dat=comp2,params=params2)
    comp2_min <- single_curve_min(dat=comp2,params=params2)
  } else if (sample){
    if (condition1){
      params1 <- params1[ss,]
      comp1_top <- single_curve_top(dat=comp1,params=params1)
      comp1_min <- single_curve_min(dat=comp1,params=params1)
    } else {
      params1 <- NA
      comp1_top <- NA
      comp1_min <- NA
      
    }
    if (condition2){
      params2 <- params2[ss,]
      comp2_top <- single_curve_top(dat=comp2,params=params2)
      comp2_min <- single_curve_min(dat=comp2,params=params2)
    } else {
      params2 <- NA
      comp2_top <- NA
      comp2_min <- NA
    }
  }
  tu_fact <- 0.7
  # eliminate any active curves with negative tops
  if ((comp1_top<0 & condition1) | (comp2_top<0 & condition2)){
    ACy_mix <- rep(NA,length(YY))
  } else if (!condition1 & !condition2){
    ACy_mix <- rep(NA,length(YY))
  } else if (!condition1 & condition2){
    ACy_mix_temp <- ifelse(YY>comp2_top,10000,-10000)
    y_indices <- which(YY>=comp2_min & YY<=comp2_top)
    comp2_curve <- c(rep(NA,length(YY)))
    comp2_curve[y_indices] <- single_curve_inv(y=YY[y_indices], single_dat=comp2, params=params2)
    ACy_mix <- ifelse((YY>=comp2_min & YY<=comp2_top),1/(cfrac2/comp2_curve),
                      ACy_mix_temp)
  } else if (condition1 & !condition2){
    ACy_mix_temp <- ifelse(YY>comp1_top,10000,-10000)
    y_indices <- which(YY>=comp1_min & YY<=comp1_top)
    comp1_curve <- c(rep(NA,length(YY)))
    comp1_curve[y_indices] <- single_curve_inv(y=YY[y_indices], single_dat=comp1, params=params1)
    ACy_mix <- ifelse((YY>=comp1_min & YY<=comp1_top),1/(cfrac1/comp1_curve),
                      ACy_mix_temp)
  } else if (condition1 & condition2){
    y_indices1 <- which(YY>=comp1_min & YY<=comp1_top)
    comp1_curve <- c(rep(NA,length(YY)))
    comp1_curve[y_indices1] <- single_curve_inv(y=YY[y_indices1], single_dat=comp1, params=params1)
    y_indices2 <- which(YY>=comp2_min & YY<=comp2_top)
    comp2_curve <- c(rep(NA,length(YY)))
    comp2_curve[y_indices2] <- single_curve_inv(y=YY[y_indices2], single_dat=comp2, params=params2)
    if (comp1_top<comp2_top){
      ec1_70 <- single_curve_inv(y=tu_fact*comp1_top, single_dat=comp1, params=params1)
      ec2_70 <- single_curve_inv(y=comp1_top*tu_fact, single_dat=comp2, params=params2)
      ACy_mix_temp <- ifelse((YY<comp1_min | YY<comp2_min),-10000,10000)
      ACy_mix_temp2 <- ifelse((YY>=comp1_min & YY>=comp2_min & YY<tu_fact*comp1_top),
                              (1/((cfrac1/comp1_curve) + (cfrac2/comp2_curve))),ACy_mix_temp)
      ACy_mix <- ifelse((YY>=tu_fact*comp1_top & YY<=comp2_top),
                        (comp2_curve/cfrac2)*(1-(cfrac1*(1/((cfrac1/ec1_70) + (cfrac2/ec2_70)))/ec1_70)),ACy_mix_temp2)
    }else if(comp2_top<comp1_top){
      ec1_70 <- single_curve_inv(y=tu_fact*comp2_top, single_dat=comp1, params=params1)
      ec2_70 <- single_curve_inv(y=comp2_top*tu_fact, single_dat=comp2, params=params2)
      ACy_mix_temp <- ifelse((YY<comp1_min | YY<comp2_min),-10000,10000)
      ACy_mix_temp2 <- ifelse((YY>=comp1_min & YY>=comp2_min & YY<tu_fact*comp2_top),
                              1/((cfrac1/comp1_curve) + (cfrac2/comp2_curve)),ACy_mix_temp)
      ACy_mix <- ifelse((YY>=tu_fact*comp2_top & YY<=comp1_top),
                        (comp1_curve/cfrac1)*(1-(cfrac2*(1/((cfrac1/ec1_70) + (cfrac2/ec2_70)))/ec2_70)),ACy_mix_temp2)
    }
  }
  return(ACy_mix)
  
}

#' Compute Independent Action mixture model for given concentration values
#' 
#' @param ss integer | index of sample number for bootstrap resamples / bayesian posterior parameter samples
#' @param val numeric | given concentration value to calculate mixture IA model at
#' @param comp1 data frame | component 1 concentration-response curve information from ToxCast mc level 5
#' @param comp2 data frame | component 2 concentration-response curve information from ToxCast mc level 5
#' @param params1 data frame | component 1 curve parameters from ToxCast, bootstrap, or Bayesian curve fit
#' @param params2 data frame | component 2 curve parameters from ToxCast, bootstrap, or Bayesian curve fit
#' @param cfrac1 numeric | mixture concentration fraction of component 1
#' @param cfrac2 numeric | mixture concentration fraction of component 1
#' @param condition1 logical | component 1 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param condition2 logical | component 2 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param norm_factor vector | normalization factor for endpoint
#' @param sample logical | if sample=TRUE, parameters are provided as a sampled sets from bootsrap of Bayesian fitting. If sample=FALSE, only one set of curve parameters are provided for each component
#'
#' @return numeric vector| mixture normalized response computed from IA model for a given concentration val
ia_full_resp <- function(ss=NA,val,comp1,comp2,params1,params2,cfrac1,cfrac2,
                         condition1,condition2,norm_factor,sample=FALSE){
  if (sample==TRUE){
    if (condition1){
      params1 <- params1[ss,]
    } else {
      params1 <- NA
    }
    if (condition2){
      params2 <- params2[ss,]
    } else {
      params2 <- NA
    }
  }
  # val is a concentration value
  curve1 <- toxcast_model(dat=comp1,params=params1,XX=val*cfrac1)
  curve2 <- toxcast_model(dat=comp2,params=params2,XX=val*cfrac2)
  # if the direction of one curve is negative, make return NA
  if((max(curve1,na.rm=T)<=0 & condition1) | (max(curve2,na.rm=T)<=0 & condition2)){
    ia_mix <- rep(NA,length(val))
  } else{
    ia_mix <- 1-((1-curve1/norm_factor)*(1-curve2/norm_factor))
  }
  # remove poly2 that go concave up negative
  if (comp1$modl=='poly2' & condition1){
    if (params1$a>0 & params1$b<0){
      ia_mix <- rep(NA,length(val))
    }
  }
  if (comp2$modl=='poly2' & condition2){
    if (params2$a>0 & params2$b<0){
      ia_mix <- rep(NA,length(val))
    }    
  }
  return(ia_mix)
}

#' Compute Concentration Addition mixture model Bayesian prediction interval samples for one given response value
#' 
#' @param ss integer | index of sample number for bootstrap resamples / bayesian posterior parameter samples
#' @param YY numeric | given non-normalized response value to calculate mixture CA model at
#' @param comp1 data frame | component 1 concentration-response curve information from ToxCast mc level 5
#' @param comp2 data frame | component 2 concentration-response curve information from ToxCast mc level 5
#' @param params1_samp data frame | component 1 curve parameters from Bayesian curve fit
#' @param params2_samp data frame | component 2 curve parameters from Bayesian curve fit
#' @param cfrac1 numeric | mixture concentration fraction of component 1
#' @param cfrac2 numeric | mixture concentration fraction of component 1
#' @param condition1 logical | component 1 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param condition2 logical | component 2 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#'
#' @return vector| prediction interval samples for mixture concentration computed from CA model for a given response YY
ca_point_pred <- function(ss=NA,YY,comp1,comp2,params1_samp,params2_samp,cfrac1,cfrac2,
                          condition1,condition2){
  if (condition1){
    params1 <- params1_samp[ss,]
    err_samp_1 <- rt(1,df=4)*params1$sigma
    comp1_top <- single_curve_top(dat=comp1,params=params1)
    comp1_min <- single_curve_min(dat=comp1,params=params1)
    y_err_1 <- YY - err_samp_1
  } else {
    params1 <- NA
    comp1_top <- NA
    comp1_min <- NA
  }
  if (condition2){
    params2 <- params2_samp[ss,]
    err_samp_2 <- rt(1,df=4)*params2$sigma
    comp2_top <- single_curve_top(dat=comp2,params=params2)
    comp2_min <- single_curve_min(dat=comp2,params=params2)
    y_err_2 <- YY - err_samp_2
  } else {
    params2 <- NA
    comp2_top <- NA
    comp2_min <- NA
  }
  tu_fact <- 0.7
  # eliminate any active curves with negative tops
  if ((comp1_top<0 & condition1) | (comp2_top<0 & condition2)){
    ACy_mix <- rep(NA,length(YY))
  } else if (!condition1 & !condition2){
    ACy_mix <- rep(NA,length(YY))
  } else if (!condition1 & condition2){
    if (y_err_2>comp2_top){
      if (comp2$modl %in% c("poly1","pow","exp2","exp3")){
        ACy_mix <- 1/(cfrac2/single_curve_inv(y=y_err_2, single_dat=comp2, params=params2))
      }else if(comp2$modl == 'poly2') {
        if(!(params2$a<0 & params2$b<0)){
          ACy_mix <- 1/(cfrac2/single_curve_inv(y=y_err_2, single_dat=comp2, params=params2))
        } else{
          ACy_mix <- 10000
        }
      }else {
        ACy_mix <- 10000
      }
    } else if (y_err_2<comp2_min){
      ACy_mix <- -10000
    } else if (y_err_2>=comp2_min & y_err_2<=comp2_top){
      ACy_mix <- 1/(cfrac2/single_curve_inv(y=y_err_2, single_dat=comp2, params=params2))
    }
  } else if (condition1 & !condition2){
    if (y_err_1>comp1_top){
      if (comp1$modl %in% c("poly1","pow","exp2","exp3")){
        ACy_mix <- 1/(cfrac1/single_curve_inv(y=y_err_1, single_dat=comp1, params=params1))
      }else if(comp1$modl == 'poly2') {
        if(!(params1$a<0 & params1$b<0)){
          ACy_mix <- ACy_mix <- 1/(cfrac1/single_curve_inv(y=y_err_1, single_dat=comp1, params=params1))
        } else{
          ACy_mix <- 10000
        }
      }else{
        ACy_mix <- 10000
      }
    } else if (y_err_1<comp1_min){
      ACy_mix <- -10000
    } else if (y_err_1>=comp1_min & y_err_1<=comp1_top){
      ACy_mix <- 1/(cfrac1/single_curve_inv(y=y_err_1, single_dat=comp1, params=params1))
    }
  } else if (condition1 & condition2){
    # component 1 curve  evaluated only when valid
    # component 2 curve evaluated only when valid
    if (y_err_1>=comp1_min & y_err_1<=comp1_top & y_err_2>=comp2_min & y_err_2<=comp2_top){
      comp1_curve <- single_curve_inv(y=y_err_1, single_dat=comp1, params=params1)
      comp2_curve <- single_curve_inv(y=y_err_2, single_dat=comp2, params=params2)
      ACy_mix <- 1/((cfrac1/comp1_curve) + (cfrac2/comp2_curve))
    } else if (y_err_1>comp1_top & y_err_2<=comp2_top & y_err_1>=comp1_min & y_err_2>=comp2_min){
      if (comp2_top<(tu_fact*comp1_top)){
        ACy_mix <- NA
      } else{
        ec1_70 <- single_curve_inv(y=tu_fact*comp1_top, single_dat=comp1, params=params1)
        ec2_70 <- single_curve_inv(y=tu_fact*comp1_top, single_dat=comp2, params=params2)
        comp2_curve <- single_curve_inv(y=y_err_2, single_dat=comp2, params=params2)
        ACy_mix <- (comp2_curve/cfrac2)*(1-(cfrac1*(1/((cfrac1/ec1_70) + (cfrac2/ec2_70)))/ec1_70))        
      }
    } else if (y_err_2>comp2_top & y_err_1<=comp1_top & y_err_1>=comp1_min & y_err_2>=comp2_min){
      if(comp1_top<(tu_fact*comp2_top)){
        ACy_mix <- NA
      } else {
        ec1_70 <- single_curve_inv(y=tu_fact*comp2_top, single_dat=comp1, params=params1)
        ec2_70 <- single_curve_inv(y=tu_fact*comp2_top, single_dat=comp2, params=params2) 
        comp1_curve <- single_curve_inv(y=y_err_1, single_dat=comp1, params=params1)
        ACy_mix <- (comp1_curve/cfrac1)*(1-(cfrac2*(1/((cfrac1/ec1_70) + (cfrac2/ec2_70)))/ec2_70))        
      }
    } else if (y_err_1<comp1_min | y_err_2<comp2_min){
      ACy_mix <- -10000
    } else if (y_err_1>comp1_top & y_err_2>comp2_top){
      ACy_mix <- 10000
    }
  }
  return(ACy_mix)
}

#' Compute Independent Action mixture model prediction interval samples for given concentration values
#' 
#' @param ss integer | index of sample number for bayesian posterior parameter samples
#' @param comp1 data frame | component 1 concentration-response curve information from ToxCast mc level 5
#' @param comp2 data frame | component 2 concentration-response curve information from ToxCast mc level 5
#' @param params1 data frame | component 1 curve parameters from Bayesian curve fit
#' @param params2 data frame | component 2 curve parameters from Bayesian curve fit
#' @param cfrac1 numeric | mixture concentration fraction of component 1
#' @param cfrac2 numeric | mixture concentration fraction of component 1
#' @param norm_factor vector | normalization factor for endpoint
#' @param x_vals vector | given concentration values to calculate mixture IA model at
#' @param condition1 logical | component 1 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#' @param condition2 logical | component 2 hitcall active/inactive evaluation, TRUE if hitc >= 0.9
#'
#' @return numeric vector| mixture normalized response computed from IA model for a given concentration val
ia_full_pred <- function(ss=NA,comp1,comp2,cfrac1,cfrac2,params1,params2,norm_factor,x_vals,condition1,condition2){
  if (condition1){
    # remove poly2 that go concave up negative
    if(comp1$modl=='poly2'){
      if (params1[ss,]$a>0 & params1[ss,]$b<0){
        curve1_err <- rep(NA,length(x_vals))
      } else{
        curve1 <- toxcast_model(dat=comp1,params=params1[ss,],XX=x_vals*cfrac1)
        if (max(curve1,na.rm=T)<=0){
          curve1_err <- rep(NA,length(x_vals))
        } else {
          curve1_err <- curve1 + rt(1,df=4)*params1[ss,]$sigma
        }        
      }
    } else {
      curve1 <- toxcast_model(dat=comp1,params=params1[ss,],XX=x_vals*cfrac1)
      if (max(curve1,na.rm=T)<=0){
        curve1_err <- rep(NA,length(x_vals))
      } else {
        curve1_err <- curve1 + rt(1,df=4)*params1[ss,]$sigma
      }
    }
  } else {
    curve1_err <- 0
  }
  if (condition2){
    # remove poly2 that go concave up negative
    if(comp2$modl=='poly2'){
      if (params2[ss,]$a>0 & params2[ss,]$b<0){
        curve2_err <- rep(NA,length(x_vals))
      } else{
        curve2 <- toxcast_model(dat=comp2,params=params2[ss,],XX=x_vals*cfrac2)
        if (max(curve2,na.rm=T)<=0){
          curve2_err <- rep(NA,length(x_vals))
        } else {
          curve2_err <- curve2 + rt(1,df=4)*params2[ss,]$sigma
        }        
      }
    } else {
      curve2 <- toxcast_model(dat=comp2,params=params2[ss,],XX=x_vals*cfrac2)
      if (max(curve2,na.rm=T)<=0){
        curve2_err <- rep(NA,length(x_vals))
      } else {
        curve2_err <- curve2 + rt(1,df=4)*params2[ss,]$sigma
      }
    }
  } else {
    curve2_err <- 0
  }  
  ia_mix <- 1-((1-curve1_err/norm_factor)*(1-curve2_err/norm_factor))
  if (comp1$hitcall<0.9 & comp2$hitc<0.9){
    return(rep(NA,length(x_vals)))
  } else {
    return(ia_mix)
  }
}

#' Calculate prediction interval samples for a single curve
#'
#' @param ss integer | index of sample number for sampled parameters
#' @param curve_dat data.frame | data frame containing mc5 level ToxCast information for curve (modl,conc_max,conc_min)
#' @param params data.frame | dataframe containing parameter values for curve type
#' @param XX vector | concentration values for which to calculate curve
#' 
#' @return vector | predicted response values calculated for curve across input concentration values used to calculate a prediction interval
single_curve_pred <- function(ss,curve_dat,params_sample,XX){
  single_curve <- toxcast_model(dat=curve_dat,params=params_sample[ss,],XX=XX)
  if(max(single_curve,na.rm=T)<=0){
    return(rep(NA,length(XX)))
  } else {
    return(single_curve + rt(1,df=4)*params_sample[ss,]$sigma)
  }
}

#' Compute fraction of modeled ACC values that are accurate or conservative
#' based on the difference in log10(ACC) values between the observed and modeled ACC point estimates
#' Threshold is set to +/- 0.5
#'
#' @param results vector| contains metrics for each analyzed mixture of observed ACC - modeled ACC difference for a given mixture model
#' @param conserv logical | if TRUE, returns fraction of analyzed mixtures with a conservative model estimate. If FALSE, returns fraction of analyzed mixtures with an accurate model estimate
#' 
#' @return numeric | fraction of total analyzed mixtures with a conservative or accurate ACC model estimate
point_fit_metric <- function(results,conserv=T){
  if(conserv){
    return(length(results[results >= -0.5 & !is.na(results)])/length(results[!is.na(results)]))
  }
  else{
    return(length(results[results <= 0.5 & results >= -0.5 & !is.na(results)])/length(results[!is.na(results)]))
  }
}

#' Calculate 95% interval around ACC from boostrapped or Bayesian ACC samples
#'
#' @param results vector| contains bootstrap or Bayesian ACC samples
#' 
#' @return lower and upper 95% boostrap confidence or Bayesian credible interval bounds
accint_quantiles <- function(results){
  # omit + infinity placeholder values (ACC=10000) from interval calculation
  # results_crop <- ifelse(results==10000,NA, results)
  int_95 <- quantile(results,probs=c(0.025,0.975),na.rm=TRUE)
  return(int_95)
}

#' Compute fraction of Residual Sum of Squares (RSS) ratio metrics of
#' RSS_modeled:RSS_observed that fall within a threshold
#' Threshold is set to <=10
#'
#' @param results vector| contains metrics for each analyzed mixture of RSS_modeled:RSS_observed for a given mixture model
#' 
#' @return numeric | fraction of total analyzed mixtures with a full curve model that represents the observed data as 
#' well within 10-fold or better than the experimentally-derived mixture curve
rss_fit_metric <- function(results){
  return(length(results[results<=10 & !is.na(results)])/length(results[!is.na(results)]))
}


#' Compute fraction of 95% prediction intervals of the analyzed mixtures that cover
#' the observed data points for a set threshold 
#' Threshold is set to 80% coverage
#'
#' @param results vector| contains metrics for each analyzed mixture of prediction interval coverage for a given mixture model
#' 
#' @return numeric | fraction of total analyzed mixtures with a full curve prediction interval that captures the observed data with at least 80% coverage
pred_fit_metric <- function(results){
  return(length(results[results>=0.8 & !is.na(results)])/length(results[!is.na(results)]))
}

#' Compute fraction of the modeled 95% bootstrapped confidence or Bayesian credible intervals
#' of the analyzed mixtures that overlap the observed mixture confidence/credible intervals
#' at the observed concentrations for a set threshold. 
#' Threshold is set to 100% overlap for the observed concentration values
#'
#' @param results vector| contains metrics for each analyzed mixture of confidence or credible interval overlap for a given mixture model
#' 
#' @return numeric | fraction of total analyzed mixtures with a full curve confidence or credible interval that overlaps with the
#' experimentally-derived interval for at least 80% of the observed concentration values
interval_fit_metric <- function(results){
  return(sum(results==1,na.rm=TRUE)/length(results[!is.na(results)]))
  # return(length(results[results>=0.8 & !is.na(results)])/length(results[!is.na(results)]))
}

# Calculate bootstrap confidence intervals or Bayesian credible intervals model overlap with mixture curve intervals
#' Calculate interval overlap of 95% bootstrap confidence intervals or 95% Bayesian credible intervals with the corresponding 
#' experimentally-derived mixture 95% intervals at all the unique observed mixture concentration values
#'
#' @param curve_samples data.frame | upper and lower 95% intervals of the model for a given test mixture at set concentration or response values
#' @param mix_curves data.frame | upper and lower 95% intervals of the experimentally-derived observed mixture curve at set concentration values
#' @param eval_concs numeric vector | all unique observed concentration points of the test mixture for evaluating the overlap
#' @param norm_facotr numeric | response normalization factor for the given endpoint
#' @param model_type ('ca','ia','single') | name of model for which the input model interval was evaluated
#' 
#' @return numeric | percent overlap of model confidence/credible 95% interval and 
#' experimentally observed and derived mixture confidence/credible 95% interval out of all concentration evaluation points
interval_overlap <- function(curve_samples,mix_curves,eval_concs,norm_factor,model_type){
  mix_low <- approx(x=mix_curves$conc, y=mix_curves$X2.5./norm_factor, 
                    xout=eval_concs, rule=1,ties = list("ordered", mean))[[2]] 
  mix_hi <-  approx(x=mix_curves$conc,y=mix_curves$X97.5./norm_factor, 
                    xout=eval_concs, rule=1,ties = list("ordered", mean))[[2]]
  if(model_type=='ca'){
    ca_info_1 <- subset(curve_samples,(X2.5.>=0) & (X2.5.<=1000) &
                          !(is.na(X2.5.)) & !(is.na(X97.5.)))
    ca_info_2 <- subset(curve_samples,(X97.5.>=0) & (X97.5.<=1000) &
                          !(is.na(X2.5.)) & !(is.na(X97.5.)))
    ca_hi <- approx(x=ca_info_1$X2.5.,y=ca_info_1$resp/norm_factor,
                    xout=eval_concs, rule=1,ties = list("ordered", mean))[[2]] 
    ca_low <-  approx(x=ca_info_2$X97.5.,y=ca_info_2$resp/norm_factor,
                      xout=eval_concs, rule=2,ties = list("ordered", mean))[[2]]   
    overlap <- ifelse(mix_hi<ca_low | mix_low>ca_hi,0,1)
  } else if(model_type=='ia'){
    ia_low <- approx(x=curve_samples$conc,y=curve_samples$X2.5., 
                     xout=eval_concs, rule=2,ties = list("ordered", mean))[[2]] 
    ia_hi <-  approx(x=curve_samples$conc,y=curve_samples$X97.5., 
                     xout=eval_concs, rule=2,ties = list("ordered", mean))[[2]]
    overlap <- ifelse(mix_hi<ia_low | mix_low>ia_hi,0,1)
  } else if(model_type=='single'){
    sing_low <- approx(x=curve_samples$conc,y=curve_samples$X2.5./norm_factor, 
                       xout=eval_concs, rule=2,ties = list("ordered", mean))[[2]] 
    sing_hi <-  approx(x=curve_samples$conc,y=curve_samples$X97.5./norm_factor, 
                       xout=eval_concs, rule=2,ties = list("ordered", mean))[[2]]
    overlap <- ifelse(mix_hi<sing_low | mix_low>sing_hi,0,1)
  }
  return(overlap)
}

#' Compute Concentration Addition mixture model for mixtures with more than 2 components at one given response value
#' CA model is not extrapolated and assumes the response value exists on the curve
#' 
#' @param YY numeric | given non-normalized response value at which to calculate mixture CA model
#' @param comps data frame | tcpl mc level 5 curve information for each component in mixture
#' @param cfracs numeric vector | concentration fractions of each component
#'
#' @return numeric vector| mixture concentration computed from CA model for one given response YY and one given set of concentration fractions of the components that make up the mixture
ca_point_conc_multi <- function(cfracs,YY,comps){
  comp_curve <- rep(NA,nrow(comps))
  if (sum(comps$hitcall < 0.9)==nrow(comps)){
    ACy_mix <- NA
  } else {
    for (qq in 1:nrow(comps)){
      # comp_top <- single_curve_top(dat=comps[qq,],params=comps[qq,])
      # comp_min <- single_curve_min(dat=comps[qq,],params=comps[qq,])
      if (comps$hitcall[qq] < 0.9){
        comp_curve[qq] <- Inf
      } else {
        comp_curve[qq] <- single_curve_inv(y=YY, single_dat=comps[qq,], params=comps[qq,])
      }
    }
    ac_units <- cfracs/comp_curve
    ACy_mix <- 1/(sum(ac_units))   
  }
  return(ACy_mix)
}

