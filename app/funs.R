make_data <- function(ddata, cdata, bw){
  doses <- ddata$Dose
  doses[is.na(doses)] <- 0
  
  datehour_start <- lubridate::dmy_hm(paste(ddata$Date, ddata$Hour_start))
  datehour_end <- lubridate::dmy_hm(paste(ddata$Date, ddata$Hour_end))
  
  index_start <- is.na(datehour_start) & is.na(datehour_end) & doses != 0
  datehour_start[index_start] <- lubridate::dmy_hm(paste(ddata$Date[index_start], "00:00"))
  
  index_end <- !is.na(datehour_start) & is.na(datehour_end) & doses != 0
  datehour_end[index_end] <- datehour_start[index_end] + lubridate::hours(3)
  
  dur_infusion <- as.double.difftime(datehour_end - datehour_start, units = "hours")
  
  timecol <- NULL
  if(all(is.na(datehour_start))){
    timecol <- rep(0, length(datehour_start))
  }
  
  datehour_sample <- lubridate::dmy_hm(paste(cdata$Date, cdata$Hour_sample))
  
  ans <- adm_rows(
    time = timecol,
    cmt = 1,
    evid = 4,
    amt = doses,
    .datehour = datehour_start,
    DAY = ddata$DAY
  ) %>%
    dplyr::arrange(DAY) %>% 
    dplyr::mutate(rate = amt / dur_infusion, .after = "amt") %>%
    obs_rows(
      cmt = 1, 
      mdv = as.integer(cdata$Exclude), 
      DV = cdata$Concentration,
      DAY = cdata$DAY,
      .datehour = datehour_sample
    ) %>% 
    dplyr::mutate(BW = bw) %>% 
    filter(!is.na(.datehour))

  if(nrow(ans) > 0){
    day_of_last_adm <- max(ans$DAY[ans$evid==4 & ans$amt > 0])
    time_of_last_adm <- max(ans$time[ans$evid == 4 & ans$amt > 0])
    
    days_of_nextadm <- ddata$DAY[ddata$DAY > day_of_last_adm]
    time_of_next_adm <- time_of_last_adm + 24 * seq_along(days_of_nextadm)
    
    ans <- ans %>% 
      adm_rows(time = time_of_next_adm, 
               evid = 4L, cmt = 1, 
               amt = 0, tinf = 3, 
               DAY = days_of_nextadm) %>% 
      obs_rows(time = tail(time_of_next_adm, 1) + 24, 
               evid = 2, cmt = 1, 
               mdv = 1, DV = NA_real_, 
               DAY = tail(days_of_nextadm, 1)
      )
  }
  
  ans
}

make_post <- function(ddata, cdata, est, auc){
  
  concdata <- isolate(cdata)
  dosdata <- isolate(ddata)
  dosdata$Dose[is.na(dosdata$Dose)] <- 0
  MWbusulfan <- 246.304
  
  day_first_dose <- min(dosdata$DAY[dosdata$Dose > 0])
  day_with_explo <- unique(concdata$DAY[!is.na(concdata$Concentration)])
  
  # Getting days with doses in the protocol
  
  tab <- dosdata %>% 
    dplyr::select(DAY, Dose) %>% 
    dplyr::filter(!(Dose == 0 & DAY < day_first_dose)) 
  
  # Getting CL at day 0->4 with extrapolation of estimated ETAs.
  day_eta <- tibble::tibble(
    DAY = seq(0L, 4L), 
    ETAnum = seq(3L, by = 2, length.out = 5),
    ETAval = unname(get_eta(est, ETAnum))
  ) %>% 
    dplyr::mutate(ETAval = ifelse(ETAval == 0, NA, ETAval)) %>% 
    tidyr::fill(ETAval, .direction = "downup")
  
  idataeta <- get_eta(est)
  idataeta[paste0("ETA", day_eta$ETAnum)] <- day_eta$ETAval
  
  tabCL <- est$model %>% 
    data_set(
      dplyr::bind_cols(ID = 1, time = 0, DAY = 0:4, as.list(idataeta))
    ) %>% 
    zero_re() %>% 
    mrgsim(end = -1, carry_out = "DAY", output = "df", Req = "CL") %>% 
    dplyr::select(DAY, CL)
  
  tab$CL_estimated <- NA_real_
  tab$CL_estimated[tab$DAY %in% day_with_explo] <- tabCL$CL[tabCL$DAY %in% day_with_explo]
  
  tab$ETAIOVCL_used <- day_eta$ETAval[day_eta$DAY %in% tab$DAY]
  tab$ETAIOVCLnum <- day_eta$ETAnum[day_eta$DAY %in% tab$DAY]
  tab$CL_used <- tabCL$CL[tabCL$DAY %in% tab$DAY]
  
  # Calculating AUC and remaining dose to administer
  tab$AUC_exposed <- MWbusulfan * tab$Dose / tab$CL_used
  
  remaining_days <- dosdata$DAY[dosdata$DAY > max(day_with_explo) & dosdata$Dose == 0]
  
  remaining_AUC <- auc - sum(tab$AUC_exposed)
  
  tab$Dose_proposed <- NA_real_
  tab$Dose_proposed[tab$DAY %in% remaining_days] <- (remaining_AUC / length(remaining_days)) * mean(tab$CL_used[tab$DAY %in% remaining_days]) / MWbusulfan
  tab$Dose_proposed[tab$Dose_proposed < 0] <- 0
  
  tab$AUC_targeted <- MWbusulfan * tab$Dose_proposed / tab$CL_used
  tab$AUC_cumulative <- cumsum(tab$AUC_exposed + ifelse(is.na(tab$AUC_targeted), 0, tab$AUC_targeted))
  tab
}

make_plot <- function(est, post){
  # Updating ETAs IOV
  est2 <- est
  final_eta <- get_eta(est2)
  final_eta[post$ETAIOVCLnum] <- post$ETAIOVCL_used
  est2$final_eta[[1]] <- final_eta
  est2
  
  # Updating doses to come
  data2 <- get_data(est)
  data2$amt[data2$evid == 4 & data2$DAY %in% post$DAY[!is.na(post$Dose_proposed)]] <- post$Dose_proposed[!is.na(post$Dose_proposed)]
  data2$rate[data2$rate==0 & data2$amt > 0] <- data2$amt[data2$rate==0 & data2$amt > 0] / 3
  print(data2)
  # Make plot 
  est2 %>% 
    augment(data = data2, nocb = FALSE, end = max(data2$time)) %>% 
    plot()
}

apriori <- function(bw, malign, gsta1, auc, ndays){
  
  stopifnot(sapply(list(bw, malign, gsta1, auc, ndays), length) == 1)
  stopifnot(is.numeric(bw), is.numeric(auc), is.numeric(ndays))
  stopifnot(is.logical(malign), is.logical(gsta1))
  stopifnot(ndays %in% c(4,5))
  
  if(is.na(gsta1)){
    cl01 <- 6.16 * (bw/25)^0.821 * 0.810^malign
    cl234 <- cl01 * 0.942
  } else {
    cl01 <- 6.38 * (bw/25)^0.786 * 0.896^malign * 0.894^gsta1 
    cl234 <- cl01 * 0.951
  }
  
  daily_cl <- c(rep(cl01, ndays-3), rep(cl234, 3))
  daily_auc <- rep(auc/ndays, ndays)
  
  ans <- data.frame(
    DAY = seq.int(to = 4, by = 1, length.out = ndays), 
    CL = signif(daily_cl, 3), 
    AUC = as.integer(daily_auc), 
    DOSE = as.integer(daily_cl * daily_auc / 246.304)
  )
  names(ans) <- c("DAY", "CL (L/h)", "AUC (micromolar.h)", "DOSE (mg)")
  ans
}
