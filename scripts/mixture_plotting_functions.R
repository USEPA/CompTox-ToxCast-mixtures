#plotting function for single line CA and IA fits
plotlines.func <- function(index, mix_input, tc=FALSE){
  # ii <- order_mixresults[index]
  ii <- index # do not order
  if (mix.results.all$m4id_mix[ii] %in% mix_input$m4id_mix){
    make_plot <- TRUE
    norm_factor <- mix_input[which(mix_input$m4id_mix==mix.results.all$m4id_mix[ii]),]$edpt_top
    coff <- mix_input[which(mix_input$m4id_mix==mix.results.all$m4id_mix[ii]),]$coff
    mixdata <-   mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==mix.results.all$m4id_mix[ii])])
    mixinfo <- mc5[which(mc5$m4id==mix.results.all$m4id_mix[ii]),]
    # mixture concentration at top response
    x_vals <- seq(min(mixdata$conc),max(mixdata$conc),length.out=800)
    mix_model <- toxcast_model(dat=mixinfo,params=mixinfo,XX=x_vals)
    mix_top_loc <- x_vals[which.max(mix_model)]
    if (tc == FALSE){
      # plot single curve from minimum to maximum concentration where the single component curve exists, as determined by the tested concentrations of that curve
      sing_max <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id[ii]),]$conc_max
      sing_min <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id[ii]),]$conc_min
      sing_info <- subset(mix.results.all$potent_single[[ii]],conc>=sing_min & conc<=sing_max)
      # plot all of CA curve
      ca_info <- mix.results.all$ca_curve[[ii]]
      # plot IA up to maximum response of IA curve or concentration of maximum observed mixture response
      ia_top_loc <- mix.results.all$ia_curve[[ii]]$conc[which.max(mix.results.all$ia_curve[[ii]]$resp)]
      ia_info <- subset(mix.results.all$ia_curve[[ii]], conc<=max(mix_top_loc,ia_top_loc))
      plot_title <- "Modeled Response (Test Components)"
    } else {
      # plot single curve from minimum to maximum concentration where the single component curve exists, as determined by the tested concentrations of that curve
      sing_max <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id.tc[ii]),]$conc_max
      sing_min <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id.tc[ii]),]$conc_min
      sing_info <- subset(mix.results.all$potent_single.tc[[ii]],conc>=sing_min & conc<=sing_max)
      # plot all of CA curve
      ca_info <- mix.results.all$ca_curve.tc[[ii]]
      # plot IA up to maximum response of IA curve or concentration of maximum observed mixture response
      ia_top_loc <- mix.results.all$ia_curve.tc[[ii]]$conc[which.max(mix.results.all$ia_curve.tc[[ii]]$resp)]
      ia_info <- subset(mix.results.all$ia_curve.tc[[ii]], conc<=max(mix_top_loc,ia_top_loc))
      plot_title <- "Modeled Response (Legacy Components)"
    }
  } else {
    make_plot <- FALSE
  }
  if (make_plot){
    # set x-axis limits
    lower_limit <- min(mixdata$conc)
    upper_limit <- max(mixdata$conc)
    # remove NA entries and parts of curve outside x-limits
    ca_info <- subset(ca_info,((conc != -10000) & !(is.na(conc)) & (conc != 10000)))
    # change CA where conc == 0 to a small number to be able to use a log10 scale x-axis
    ca_info$conc <- ifelse(ca_info$conc<=0,1E-50,ca_info$conc)
    ia_info <- subset(ia_info,!is.na(resp))
    sing_info <- subset(sing_info,!is.na(resp))
    # remove NA entries and parts of curve outside x-limits
    # ca_info <- subset(ca_info,((conc != -10000) & !(is.na(conc)) & (conc != 10000)) & (conc>=lower_limit) & (conc<=(upper_limit+100)))
    # ia_info <- subset(ia_info,!is.na(resp) & (conc>=lower_limit) & (conc<=(upper_limit+100)))
    # sing_info <- subset(sing_info,!is.na(resp) & (conc>=lower_limit) & (conc<=(upper_limit+100)))
    # calculate mixture ToxCast curve to plot
    XX <- 10^seq(log10(min(mixdata$conc)),log10(max(mixdata$conc)),length.out=300)
    mix_plot <- data.frame(conc=XX,resp=toxcast_model(dat=mixinfo,params=mixinfo,XX=XX)/norm_factor)
    # set y-axis units
    ggplot() +
      # Plot Mixture ToxCast curve
      geom_line(data=mix_plot,aes(x=conc,y=resp,color="Mixture"),lty=1,lwd=0.5)+
      # Plot CA and IA full models
      geom_line(data=ca_info,aes(x=conc,y=resp/norm_factor,color="CA"),lty=1,lwd=0.5)+
      geom_line(data=ia_info,aes(x=conc,y=resp,color="IA"),lty=6,lwd=0.5)+
      # plot most potent single curve
      geom_line(data=sing_info,aes(x=conc,y=resp/norm_factor,color="Single"),lty=2,lwd=0.5)+
      # Plot cutoff
      geom_hline(yintercept=coff/norm_factor, linetype=2, color="gray34") +
      # Plot mixture data
      geom_point(data = mixdata, aes(x = mixdata$conc, y = mixdata$resp/norm_factor,color="Mixture"), shape=19, size=1) + 
      scale_x_log10(labels = label_number(drop0trailing = TRUE))+
      coord_cartesian(xlim=c(lower_limit,upper_limit))+
      theme_bw(base_size=6) +
      labs(x=expression(paste("Concentration [",mu, "M]")), y="Normalized Response",
           title=plot_title) + #, tag = paste("rmse[CA]: ", round(fitvals[[ii]]$rmse_ca, digits=2), '\n', "rmse[IA]: ", round(fitvals[[ii]]$rmse_ia, digits=2))) +
      scale_color_manual(values=c(CA="#ffb000",IA="#006299",Single="#be548f",Mixture="black"),
                         breaks=c("Mixture","CA","IA","Single"),
                         labels=c(CA="CA",IA="IA",Single="MP",Mixture="Observed Mixture"))+
      theme(legend.title=element_blank(),
            legend.background = element_rect(fill=NA),
            legend.margin = margin(0.1,0.1,0.1,0.1, unit='cm'),
            legend.text = element_text(size=6),
            legend.justification = c(0,1),
            legend.position = "inside",
            legend.position.inside = c(0.01,0.99),
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.spacing.y = unit(0, "cm"),
            legend.spacing.x = unit(0, "cm"),
            axis.text = element_text(color="black"),
            plot.margin=unit(c(0.1,0.6, 0.1,0.6), "cm"),
            plot.title=element_text(size=9, hjust=0.5, face="bold", margin=margin(0,0,1,0)))
  } else{
    ggplot() + theme_void()
  }
}



#plotting bootstrap intervals CA and IA fits
plotint.func <- function(index,tc=FALSE){
  # ii <- order_mixresults[index]
  ii <- index # do not order
  if(!tc){
    if (mix.results.all$m4id_mix[ii] %in% mix.info$m4id_mix){
      make_plot <- TRUE
      norm_factor <- mix.info[which(mix.info$m4id_mix==mix.results.all$m4id_mix[ii]),]$edpt_top
      coff <- mix.info[which(mix.info$m4id_mix==mix.results.all$m4id_mix[ii]),]$coff
      mixdata <-   mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==mix.results.all$m4id_mix[ii])])
      mixinfo <- mc5[which(mc5$m4id==mix.results.all$m4id_mix[ii]),]
      # mixture concentration at top response
      x_vals <- seq(min(mixdata$conc),max(mixdata$conc),length.out=800)
      mix_model <- toxcast_model(dat=mixinfo,params=mixinfo,XX=x_vals)
      mix_top_loc <- x_vals[which.max(mix_model)]
      # plot all of CA curve interval
      ca_info <- mix.results.all$ca_bootint[[ii]]
      # plot IA interval up to the IA model maximum
      ia_max_index <- min(which.max(mix.results.all$ia_bootint[[ii]]$X2.5.),
                          which.max(mix.results.all$ia_bootint[[ii]]$X97.5.))
      ia_info <- mix.results.all$ia_bootint[[ii]]
      # plot single curve interval from minimum to maximum concentration where the single component curve exists, as determined by the tested concentrations of that curve
      sing_max <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id[ii]),]$conc_max
      sing_min <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id[ii]),]$conc_min
      single_info <- subset(mix.results.all$single_bootint[[ii]],conc>=sing_min & conc<=sing_max)
      plot_title <- "Bootstrap Intervals (Test Components)"
    } else {
      make_plot <- FALSE
    }
    
  }
  if(tc){
    if (mix.results.all$m4id_mix[ii] %in% mix.info.tc$m4id_mix){
      make_plot <- TRUE
      norm_factor <- mix.info.tc[which(mix.info.tc$m4id_mix==mix.results.all$m4id_mix[ii]),]$edpt_top
      coff <- mix.info.tc[which(mix.info.tc$m4id_mix==mix.results.all$m4id_mix[ii]),]$coff
      mixdata <-   mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==mix.results.all$m4id_mix[ii])])
      mixinfo <- mc5[which(mc5$m4id==mix.results.all$m4id_mix[ii]),]
      # mixture concentration at top response
      x_vals <- seq(min(mixdata$conc),max(mixdata$conc),length.out=800)
      mix_model <- toxcast_model(dat=mixinfo,params=mixinfo,XX=x_vals)
      mix_top_loc <- x_vals[which.max(mix_model)]
      # plot all of CA interval
      ca_info <- mix.results.all$ca_bootint.tc[[ii]]
      # plot IA interval up to the IA maximum
      ia_max_index <- min(which.max(mix.results.all$ia_bootint.tc[[ii]]$X2.5.),
                          which.max(mix.results.all$ia_bootint.tc[[ii]]$X97.5.))
      # ia_info <- mix.results.all$ia_bootint.tc[[ii]][1:ia_max_index,]
      ia_info <- mix.results.all$ia_bootint.tc[[ii]]
      # plot single curve interval from minimum to maximum concentration where the single component curve exists, as determined by the tested concentrations of that curve
      sing_max <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id.tc[ii]),]$conc_max
      sing_min <- mc5[which(mc5$m4id==mix.results.all$potent_sing_id.tc[ii]),]$conc_min
      single_info <- subset(mix.results.all$single_bootint.tc[[ii]],conc>=sing_min & conc<=sing_max)
      plot_title <- "Bootstrap Intervals (Legacy Components)"
    } else {
      make_plot <- FALSE
    }
  }
  if(make_plot){
    # set x-axis limits
    lower_limit <- min(mixdata$conc)
    upper_limit <- max(mixdata$conc) #max(c(mixdata$conc,ia_info$conc,ca_info$X97.5.))
    # calculate mixture ToxCast curve to plot
    XX <- 10^seq(log10(min(mixdata$conc)),log10(max(mixdata$conc)),length.out=300)
    mix_plot <- data.frame(conc=XX,resp=toxcast_model(dat=mixinfo,params=mixinfo,XX=XX)/norm_factor)
    # remove NA entries and entries outside plot limits
    # ca_info_1 <- subset(ca_info,(X2.5.>0) & (X2.5.<=1000)&
    #                       !(is.na(X2.5.)) & !(is.na(X97.5.)))
    # ca_info_2 <- subset(ca_info,(X97.5.>0) & (X97.5.<=1000) &
    #                       !(is.na(X2.5.)) & !(is.na(X97.5.)))
    # ca_info <- subset(ca_info,(X2.5.>0) & (X97.5.<=1000) &
    #                     !(is.na(X2.5.)) & !(is.na(X97.5.)))
    # select part of interval that has an upper bound of + infinity (top not in evaluated range)
    # ca_info_extra <- subset(ca_info,(X97.5.>1000 | X97.5.==max(ca_info$X97.5[ca_info$X97.5<=1000],na.rm=T)))
    # remove NA entries and upper concentration entries outside plot limits, will cut off interval plot at top of upper 95% interval which is not an accurate limit but is a good visual limit
    ca_info <- subset(ca_info, !(is.na(X2.5.)) & !(is.na(X97.5.)) & (X97.5.<=1000))
    # ca_info <- subset(ca_info, !(is.na(X2.5.)) & !(is.na(X97.5.)))
    # change negative values at lower extreme to a very small number for log transformation plotting
    ca_info$X2.5. <- ifelse(ca_info$X2.5.<=0,1E-50,ca_info$X2.5.)
    ca_info$X97.5. <- ifelse(ca_info$X97.5.<=0,1E-50,ca_info$X97.5.)
    # ca_info$X97.5. <- ifelse(ca_info$X97.5>1000,max(mixdata$conc),ca_info$X97.5.)
    # remove negative values from CA lower curve
    # ia_info <- subset(ia_info,!is.na(X2.5.) & !is.na(X97.5.) & (conc>0) & (conc<=1000))
    max_ia_conc <- max(ia_info$conc[ia_max_index],mix_top_loc)
    ia_info <- subset(ia_info, conc<=max_ia_conc & !is.na(X2.5.) & !is.na(X97.5.) & (conc>0))
    
    mix_int <- subset(mix.results.all$mix_bootint[[ii]],
                      conc>=lower_limit & conc<=upper_limit)
    single_info <- subset(single_info, !(is.na(X2.5.)) & !(is.na(X97.5.)))

    p1 <- ggplot() +
      scale_x_log10(labels = label_number(drop0trailing = TRUE))+
      coord_cartesian(xlim=c(lower_limit,upper_limit))+
      # plot mixture interval
      geom_ribbon(data=mix_int,
                  aes(x=conc,ymin=X2.5./norm_factor,ymax=X97.5./norm_factor,fill="Mix_Int"),color="black",linewidth=0.25,alpha=0.4)+
      # Plot IA model
      geom_ribbon(data=ia_info,
                  aes(x=conc,ymin=X2.5.,ymax=X97.5.,fill="IA"),color="#006299", linewidth=0.25, alpha=0.8)+
      # Plot CA models
      geom_ribbon(data=ca_info,
                  aes(y=resp/norm_factor,xmin=X2.5.,xmax=X97.5.,fill="CA"),color="#ffb000",linewidth=0.25,alpha=0.5)+
      # # plot extra part of CA model
      # geom_ribbon(data=ca_info_extra,aes(y=resp/norm_factor,xmin=X2.5.,xmax=rep(max(ca_info$X97.5.),nrow(ca_info_extra)),fill="CA"),color="#ffb000",linewidth=0.25,alpha=0.5, show.legend=FALSE)+
      # geom_line(data=ca_info_1,aes(x=X2.5.,y=resp/norm_factor,color="CA"),lty=1,lwd=0.5)+
      # geom_line(data=ca_info_2,aes(x=X97.5.,y=resp/norm_factor,color="CA"),lty=1,lwd=0.5)+

      # geom_line(data=ia_info,aes(x=conc,y=X2.5.,color="IA"),lty=6,lwd=0.5)+
      # geom_line(data=ia_info,aes(x=conc,y=X97.5.,color="IA"),lty=6,lwd=0.5)+
      # Plot most potent single chemical model
      geom_ribbon(data=single_info,
                  aes(x=conc,ymin=X2.5./norm_factor,ymax=X97.5./norm_factor,fill="Single"),linewidth=0.25,color="#be548f",alpha=0.2)+
      # geom_line(data=single_info,aes(x=conc, y=X2.5./norm_factor,color="Single"),lty=2,lwd=0.5)+
      # geom_line(data=single_info,aes(x=conc, y=X97.5./norm_factor,color="Single"),lty=2,lwd=0.5)+

      # Plot cutoff
      geom_hline(yintercept=coff/norm_factor, linetype="dashed", color="gray34") +
      # Plot Mixture ToxCast curve
      geom_line(data=mix_plot,aes(x=conc,y=resp,color="Mixture"),lty=1,lwd=0.5)+
      # Plot mixture data
      geom_point(data = mixdata, aes(x = mixdata$conc, y = mixdata$resp/norm_factor, color="Mixture"), shape=19, size=1) + 
      theme_bw(base_size=6) +
      labs(x=expression(paste("Concentration [",mu, "M]")), y="Normalized Response",
           title=plot_title) + 
      scale_fill_manual(values=c(CA="#ffb000",IA="#006299",Single="#be548f",Mix_Int="gray40"),
                        breaks=c("Mix_Int","CA","IA","Single"),
                        labels=c(CA="CA 95% Confidence Interval",IA="IA 95% Confidence Interval",Single="MP 95% Confidence Interval",Mix_Int="Observed Mixture 95% Confidence Interval"))+
      scale_color_manual(values=c(Mixture="black"),breaks=c("Mixture"),labels=c(Mixture="Observed Mixture"))+
      theme(legend.title=element_blank(),
            legend.background = element_rect(fill=NA),
            legend.margin = margin(0.1,0.1,0.1,0.1, unit='cm'),
            legend.text = element_text(size=6),
            legend.justification = c(0,1),
            legend.position = "inside",
            legend.position.inside=c(0.01,0.99),
            legend.key.width = unit(0.5, "cm"),
            legend.key.height = unit(0.13, "cm"),
            legend.spacing.y = unit(-0.15, "cm"),
            # legend.spacing.x = unit(0, "cm"),
            axis.text = element_text(color="black"),
            plot.margin=unit(c(0.1,0.6, 0.1,0.6), "cm"),
            plot.title=element_text(size=9, hjust=0.5, face="bold", margin=margin(0,0,1,0)))
    
  } else{
    ggplot() + theme_void()
  }
}


plotbayesintia.func <- function(index,tc=FALSE){
  if (tc==FALSE){
    if (mix.results.all$m4id_mix[index] %in% mix.info$m4id_mix){
      make_plot <- TRUE
      plot_title <- "Bayesian IA Model (Test Components)"
      cred_max_index <- min(which.max(mix.results.all$ia_bayesint[[index]]$X2.5.),
                            which.max(mix.results.all$ia_bayesint[[index]]$X97.5.))
      pred_max_index <- min(which.max(mix.results.all$ia_predint[[index]]$X2.5.),
                            which.max(mix.results.all$ia_predint[[index]]$X97.5.))
      # plot_cred <- mix.results.all$ia_bayesint[[index]][1:cred_max_index,]
      # plot_pred <- mix.results.all$ia_predint[[index]][1:pred_max_index,]
      plot_cred <- mix.results.all$ia_bayesint[[index]]
      plot_pred <- mix.results.all$ia_predint[[index]]
      norm_factor <- mix.info[which(mix.info$m4id_mix==mix.results.all$m4id_mix[index]),]$edpt_top
    } else {
      make_plot <- FALSE
    }
  } else if (tc == TRUE){
    if (mix.results.all$m4id_mix[index] %in% mix.info.tc$m4id_mix){
      make_plot <- TRUE
      plot_title <- "Bayesian IA Model (Legacy Components)"
      cred_max_index <- min(which.max(mix.results.all$ia_bayesint.tc[[index]]$X2.5.),
                            which.max(mix.results.all$ia_bayesint.tc[[index]]$X97.5.))
      pred_max_index <- min(which.max(mix.results.all$ia_predint.tc[[index]]$X2.5.),
                            which.max(mix.results.all$ia_predint.tc[[index]]$X97.5.))
      # plot_cred <- mix.results.all$ia_bayesint.tc[[index]][1:cred_max_index,]
      # plot_pred <- mix.results.all$ia_predint.tc[[index]][1:pred_max_index,]
      plot_cred <- mix.results.all$ia_bayesint.tc[[index]]
      plot_pred <- mix.results.all$ia_predint.tc[[index]]
      norm_factor <- mix.info.tc[which(mix.info.tc$m4id_mix==mix.results.all$m4id_mix[index]),]$edpt_top
    } else {
      make_plot <- FALSE
    }
  }
  if (make_plot){
    mix_data <- mc3 %>% filter(m3id %in% mc4_agg$m3id[which(mc4_agg$m4id==mix.results.all$m4id_mix[index])])
    lower_limit <- min(mix_data$conc)
    upper_limit <- max(mix_data$conc) 
    mixinfo <- mc5[which(mc5$m4id==mix.results.all$m4id_mix[index]),]
    x_vals <- seq(min(mix_data$conc),max(mix_data$conc),length.out=800)
    mix_model <- toxcast_model(dat=mixinfo,params=mixinfo,XX=x_vals)
    mix_top_loc <- x_vals[which.max(mix_model)]
    coff <- mc5[m4id==mix.results.all$m4id_mix[index],]$coff
    # calculate mixture ToxCast curve to plot
    XX <- 10^seq(log10(min(mix_data$conc)),log10(max(mix_data$conc)),length.out=300)
    mix_plot <- data.frame(conc=XX,resp=toxcast_model(dat=mixinfo,params=mixinfo,XX=XX)/norm_factor)
    # limit credible and prediction interval plots to go up to mixture top
    max_plot_conc <- max(plot_cred$conc[cred_max_index],
                         plot_pred$conc[pred_max_index],
                         mix_top_loc)
    plot_cred <- subset(plot_cred,conc<=max_plot_conc)
    plot_pred <- subset(plot_pred,conc<=max_plot_conc)
    ggplot()+
      geom_ribbon(data=plot_pred,aes(x=conc,ymin=X2.5.,ymax=X97.5.,fill="pred"))+
      geom_ribbon(data=plot_cred, aes(x=conc,ymin=X2.5.,ymax=X97.5.,fill="cred"))+
      # Plot Mixture ToxCast curve
      geom_line(data=mix_plot,aes(x=conc,y=resp,color="Mixture"),lty=1,lwd=0.5)+
      # Plot mixture data
      geom_point(data=mix_data, aes(x = conc, y = resp/norm_factor,color="Mixture"), size=1)+
      scale_fill_manual(values=c(pred="#CDEDFF",cred="#0082CC"),labels=c(pred="95% Prediction Interval",cred="95% Credible Interval"))+
      scale_color_manual(values=c(Mixture="black"),breaks=c("Mixture"),labels=c(Mixture="Observed Mixture"))+
      coord_cartesian(xlim=c(lower_limit,upper_limit))+
      scale_x_log10(labels = label_number(drop0trailing = TRUE))+
      theme_bw(base_size=6) +
      geom_hline(yintercept=coff/norm_factor, linetype="dashed", color="gray34") +
      labs(x=expression(paste("Concentration [",mu, "M]")), y="Normalized Response",
           title=plot_title)+
      theme(legend.title=element_blank(),
            legend.background = element_rect(fill=NA),
            legend.text = element_text(size=6),
            legend.margin = margin(0.1,0.1,0.1,0.1, unit='cm'),
            legend.justification = c(0,1),
            legend.position = "inside",
            legend.position.inside=c(0.01,0.99),
            legend.key.width = unit(0.6, "cm"),
            legend.key.height = unit(0.2, "cm"),
            legend.spacing.y = unit(-0.15, "cm"),
            axis.text = element_text(color="black"),
            plot.margin=unit(c(0.1,0.6, 0.1,0.6), "cm"),
            plot.title=element_text(size=9,hjust=0.5,face="bold",margin=margin(0,0,1,0)))
  } else{
    ggplot() + theme_void()
  }
}


#plotting function for single chemicals and data
plotdata.func <- function(index,tc_lib=FALSE){
  ii<- index
  mixinfo <- mc5[m4id == mix.results.all$m4id_mix[ii],]
  if (tc_lib == TRUE){
    aeid_mix <- mc5[m4id == mix.results.all$m4id_mix[ii],]$aeid
    dtxsid <- filter(mixtures_key, DTXSID.x == mc5[m4id == mix.results.all$m4id_mix[ii],]$dsstox_substance_id)$DTXSID.y
    m4id_1 <- filter(mc5, dsstox_substance_id == dtxsid[1] & aeid == aeid_mix & (spid %in% spids_tc))$m4id
    m4id_2 <- filter(mc5, dsstox_substance_id == dtxsid[2] & aeid == aeid_mix & (spid %in% spids_tc))$m4id
    component1 <- mc5 %>% filter(m4id == m4id_1)
    comp1_data <- (mc3 %>% filter(spid == component1$spid & aeid == component1$aeid))
    component2 <- mc5 %>% filter(m4id == m4id_2)
    comp2_data <- (mc3 %>% filter(spid == component2$spid & aeid == component2$aeid))
    plot_title <- "Legacy Single Component + Observed Mixture Data"
  } else {
    aeid_mix <- mc5[m4id == mix.results.all$m4id_mix[ii],]$aeid
    dtxsid <- filter(mixtures_key, DTXSID.x == mc5[m4id == mix.results.all$m4id_mix[ii],]$dsstox_substance_id)$DTXSID.y
    m4id_1 <- filter(mc5, dsstox_substance_id == dtxsid[1] & aeid == aeid_mix & (spid %in% spids))$m4id
    m4id_2 <- filter(mc5, dsstox_substance_id == dtxsid[2] & aeid == aeid_mix & (spid %in% spids))$m4id
    component1 <- mc5 %>% filter(m4id == m4id_1)
    comp1_data <- (mc3 %>% filter(spid == component1$spid & aeid == component1$aeid))
    component2 <- mc5 %>% filter(m4id == m4id_2)
    comp2_data <- (mc3 %>% filter(spid == component2$spid & aeid == component2$aeid))
    plot_title <- "Test Single Component + Observed Mixture Data"
  }
  # gather data
  df1 <- data.frame(resp=comp1_data$resp, conc=comp1_data$conc, category="comp1")
  df2 <- data.frame(resp=comp2_data$resp, conc=comp2_data$conc, category="comp2")
  mixdata <- (mc3 %>% filter(spid == mc5$spid[which(mc5$m4id==mix.results.all$m4id_mix[ii])] & aeid == component1$aeid))
  df3 <- data.frame(resp=mixdata$resp, conc=mixdata$conc, category="mix")
  data_pts <- rbind(df1,df2,df3)
  # setup curves to plot
  plot_colors <- c("Component 1" = "#D55E00", "Component 2" = "#009E73", "Mixture" = "black", '#D55E00'='#D55E00', '#009E73'='#009E73', 'black'='black' )
  curve1 <- data.frame(x=seq(min(comp1_data$conc),max(comp1_data$conc),length.out=300),
                       resp=toxcast_model(dat=component1,params=component1,
                                          XX=seq(min(comp1_data$conc),max(comp1_data$conc),length.out=300)))
  curve2 <- data.frame(x=seq(min(comp2_data$conc),max(comp2_data$conc),length.out=300),
                       resp=toxcast_model(dat=component2,params=component2,
                                          XX=seq(min(comp2_data$conc),max(comp2_data$conc),length.out=300)))
  curvemix <- data.frame(x=seq(min(mixdata$conc),max(mixdata$conc),length.out=300),
                         resp=toxcast_model(dat=mixinfo,params=mixinfo,
                                            XX=seq(min(mixdata$conc),max(mixdata$conc),length.out=300)))
  # make plot
  ggplot() +
    geom_point(data=data_pts, aes(x = conc, y = resp,
                                  color=category, shape=category,group=category), size=1) + 
    geom_point(data=data_pts, aes(x = conc, y = resp,
                                  color=category, shape=category,group=category), size=1) +
    geom_point(data=data_pts, aes(x = conc, y = resp,
                                  color=category, shape=category,group=category), size=1) +
    # geom_line(data = curve1, aes(x=x,y=resp),color="#D55E00", lwd=0.5, lty=5) +
    # geom_line(data = curve2, aes(x=x,y=resp),color="#009E73", lwd=0.5, lty=5) +
    # geom_line(data = curvemix, aes(x=x,y=resp),color="black", lwd=0.5, lty=1) +
    geom_line(data = curve1, aes(x=x,y=resp,color="comp1"), lwd=0.5, lty=5) +
    geom_line(data = curve2, aes(x=x,y=resp,color="comp2"), lwd=0.5, lty=5) +
    geom_line(data = curvemix, aes(x=x,y=resp,color="mix"), lwd=0.5, lty=1) +
    geom_hline(yintercept=component1$coff, linetype="dashed", color="gray34") +
    scale_x_log10(labels = label_number(drop0trailing = TRUE),
                  limits= c(min(comp1_data$conc,comp2_data$conc,mixdata$conc),
                            max(comp1_data$conc,comp2_data$conc,mixdata$conc)))+
    theme_bw(base_size=6) +
    labs(x=expression(paste("Concentration [",mu, "M]")), y="Response [Log2 Fold Induction]",
         title=plot_title)+
    # scale_color_manual(values=c("#D55E00","#009E73","black"),
    #                    labels=c(paste0(component1$chnm," (",component1$modl,", hitc=",signif(component1$hitc,3),")"),
    #                             paste0(component2$chnm," (",component2$modl,", hitc=",signif(component2$hitc,3),")"),
    #                             paste0("Mixture (",mixinfo$modl,", hitc=",signif(mixinfo$hitc,3),")")),
    #                    breaks=c("comp1", "comp2", "mix"))+
    scale_color_manual(values=c(comp1="#D55E00",comp2="#009E73",mix="black"),
                       labels=c(comp1=paste0(component1$chnm," (",component1$modl,", hitc=",signif(component1$hitc,3),")"),
                                comp2=paste0(component2$chnm," (",component2$modl,", hitc=",signif(component2$hitc,3),")"),
                                mix=paste0("Observed Mixture (",mixinfo$modl,", hitc=",signif(mixinfo$hitc,3),")")),
                       breaks=c("mix", "comp1", "comp2"))+
    scale_shape_manual(values=c(comp1=15,comp2=17,mix=19),
                       labels=c(comp1=paste0(component1$chnm," (",component1$modl,", hitc=",signif(component1$hitc,3),")"),
                                comp2=paste0(component2$chnm," (",component2$modl,", hitc=",signif(component2$hitc,3),")"),
                                mix=paste0("Observed Mixture (",mixinfo$modl,", hitc=",signif(mixinfo$hitc,3),")")),
                       breaks=c("mix", "comp1", "comp2"))+
    theme(legend.title=element_blank(),
          legend.background = element_rect(fill=NA),
          # legend.box.background = element_rect(color="black", linewidth=0.25),
          legend.margin = margin(0.1,0.1,0.1,0.1, unit='cm'),
          legend.text = element_text(size=6),
          legend.justification = c(0,1),
          legend.position = "inside",
          legend.position.inside = c(0.01,0.99),
          legend.key.width = unit(0.5, "cm"),
          legend.key.height = unit(0.3, "cm"),
          legend.spacing.y = unit(0, 'cm'),
          axis.text = element_text(color="black"),
          plot.margin=unit(c(0.1,0.6, 0.1,0.6), "cm"),
          plot.tag.position = c(0.98, 0.3),
          plot.title=element_text(size=9, hjust=0.5, face="bold", margin=margin(0,0,1,0))) 
}
