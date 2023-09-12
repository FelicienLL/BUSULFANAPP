
# Models and parameters

base_model <- mread("sim/base_busulfan_simulations.cpp")

param_BASE <- list(
  param = c(
    TVCL = 5.8517,
    TVV = 18.787,
    BW_CL = 0.83354,
    BW_V = 0.92716,
    DAY234_CLV = 0.9424,
    SEX_V = 1,
    MALIGN_V = 1,
    MALIGN_CL = 1,
    GSTA1_CL = 1
  ), 
  omega_iiv = bmat(0.06780630, 0.03581010, 0.02656750), 
  omega_iov = bmat(0.01981420, 0.01322990, 0.01915770), 
  sigma = bmat(0.00370852, 0, 0)
)

param_NOGENET <- list(
  param = c(
    TVCL = 6.381000 ,
    TVV = 20.190000 ,
    BW_CL = 0.821000 ,
    BW_V = 0.907300 ,
    DAY234_CLV = 0.941500 ,
    SEX_V = 0.923800 ,
    MALIGN_V = 0.917200 ,
    MALIGN_CL = 0.809700 ,
    GSTA1_CL = 1
  ), 
  omega_iiv = bmat(0.04977540, 0.02607890, 0.02044450), 
  omega_iov = bmat(0.01772870, 0.01080470, 0.01831950), 
  sigma = bmat(0.00355669, 0, 0)
)

param_GENET <- list(
  param = c(
    TVCL = 6.15685,
    TVV = 19.3446,
    BW_CL = 0.785861,
    BW_V = 0.875741,
    DAY234_CLV = 0.950817,
    SEX_V = 0.917867,
    MALIGN_V = 1,
    MALIGN_CL = 0.895993,
    GSTA1_CL = 0.894297
  ), 
  omega_iiv = bmat(0.04977540, 0.02607890, 0.02044450), 
  omega_iov = bmat(0.01772870, 0.01080470, 0.01831950), 
  sigma = bmat(0.00355669, 0, 0)
)

update_mrg <- function(x, parlist){
  x %>% 
    param(parlist$param) %>% 
    omat(
      IIV = parlist$omega_iiv,
      IOV_DAY0 = parlist$omega_iov,
      IOV_DAY1 = parlist$omega_iov,
      IOV_DAY2 = parlist$omega_iov,
      IOV_DAY3 = parlist$omega_iov,
      IOV_DAY4 = parlist$omega_iov
    ) %>%
    smat(parlist$sigma) 
}

# ---- General values ----
MWbusulfan <- 246.304
protocol_times <- c(3.00, 3.25, 3.50, 4.00, 5.00, 6.00, 8.00)
# ---- Unit functions ----
empty_data <- function(day, amt, sampletimes = c(3.25, 5, 8), covlist = list()){
  nid <- unique(sapply(covlist, length))
  cov <- bind_cols(ID = seq_len(nid), covlist, DAY = day, OCCA = day)
  bind_rows(
    data.frame(ID = cov$ID, time = 0, evid = 1, cmt = 1, amt = amt, DV = NA_real_, mdv = 1, rate  = amt / 3), 
    data.frame(ID = cov$ID,           evid = 0, cmt = 1, amt = 0,   DV = NA_real_, mdv = 0, rate = 0) %>% 
      expand_grid(time = sampletimes)
  ) %>% 
    left_join(cov, by = "ID") %>% 
    arrange(ID, time, evid)
}

pred_CL <- function(mod, dat){
  sim <- mod %>% 
    zero_re() %>% 
    idata_set(dat) %>% 
    mrgsim_df(end = 0)
  sim$CL
}

curr_equation <- function(AGE, BSA, ...){
  unname(ifelse(AGE < 1, 80, ifelse(AGE > 7, 130, 120)) * BSA)
}

init_list <- function(n = 1){
  l <- lapply(1:24, function(x) rep(NA_real_, n))
  names(l) <- unlist(lapply(c("dose", "actual_cl", "actual_AUC", "apri_cl", "apost_cl", "pred_AUC"), paste, 1:4, sep = "_"))
  l
}

#' Treatment
#'
#' @param target_AUC cumulative AUC in micromolar.min
#' @param method1,method2,method3,method4 either "curr_equation", "nogenet_equation", "genet_equation", "TDM_bayes_1", "TDM_bayes_2", "TDM_bayes_3"
#' @param patient a data.frame with columns "BW" (body weight in kg), "BSA" (body-surface area in mg/m2), "AGE" (age in years), "DIAG" (0 for malignant disease, 1 otherwise), "GSTA1" (0 for WT, 1 for HT and HM) and "SEX" (0 for male, 1 for female)
#' @param sampletimes a vector of numeric, time after start of infusion 
#' @param mrgmodel mrgsolve model. Parameters will be updated
#' @param parlist_nocov,parlist_nogenet,parlist_genet lists of population parameters, default values are provided in the script
#' 
treatment <- function(target_AUC = 16000, 
                      method1, 
                      method2, 
                      method3, 
                      method4,
                      patient = my_patients, 
                      sampletimes = c(3.25, 5, 8), 
                      mrgmodel = base_model, 
                      parlist_nocov = param_BASE, 
                      parlist_nogenet = param_NOGENET,
                      parlist_genet = param_GENET 
){
  
  sim_mrgmodel <- update_mrg(mrgmodel, parlist = parlist_genet)
  est_mrgmodel <- update_mrg(mrgmodel, parlist = parlist_nocov)
  
  nogenet_model <-  update_mrg(mrgmodel, parlist = parlist_nogenet)
  genet_model <- sim_mrgmodel
  
  seed <- as.integer(1000000*runif(1))
  ttt <- init_list(n = length(patient[[1]]))
  
  # --- DAY 1 ---------------------------------------------
  current_day <- 1
  
  # --- Calculation of Dose 1 -----------------------------
  if(method1 == "curr_equation"){
    ttt$dose_1 <- curr_equation(BSA = patient$BSA, AGE = patient$AGE)
  }
  
  if(method1 == "genet_equation"){
    ttt$apri_cl_1 <- pred_CL(mod = genet_model, dat = c(patient, list(DAY = current_day)))
    ttt$dose_1 <- ttt$apri_cl_1 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method1 == "nogenet_equation"){
    ttt$apri_cl_1 <- pred_CL(mod = nogenet_model, dat = c(patient, list(DAY = current_day)))
    ttt$dose_1 <- ttt$apri_cl_1 * target_AUC / (4.0 * MWbusulfan)
  }
  
  # --- Administration of Day 1 ---------------------------
  data1 <- empty_data(day = current_day, amt = ttt$dose_1, covlist = patient, sampletimes = sampletimes)
  set.seed(seed)
  sim1 <- mrgsim(sim_mrgmodel, data1)
  ttt$actual_cl_1 <- unique(sim1$CL)
  ttt$actual_AUC_1 <- MWbusulfan * ttt$dose_1 / ttt$actual_cl_1 
  
  
  
  # --- DAY 2 ---------------------------------------------
  current_day <- 2
  
  # --- Calculation of Dose 2 -----------------------------
  if(method2 == "curr_equation"){
    ttt$dose_2 <- curr_equation(BSA = patient$BSA, AGE = patient$AGE)
  }
  
  if(method2 == "genet_equation"){
    ttt$apri_cl_2 <- pred_CL(mod = genet_model, dat = c(patient, list(DAY = current_day)))
    ttt$dose_2 <- ttt$apri_cl_2 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method2 == "nogenet_equation"){
    ttt$apri_cl_2 <- pred_CL(mod = nogenet_model, dat = c(patient, list(DAY = current_day)))
    ttt$dose_2 <- ttt$apri_cl_2 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method2 == "TDM_bayes_1"){
    est1 <- mapbayest(est_mrgmodel, mutate(data1, DV = sim1$DV), hessian = FALSE, select_eta = c(1,2,5,6), verbose = FALSE)
    ttt$apost_cl_1 <- get_param(est1)[["CL"]]
    ttt$pred_AUC_1 <- MWbusulfan * ttt$dose_1 / ttt$apost_cl_1 
    ttt$apri_cl_2 <- ttt$apost_cl_1 * sim_mrgmodel$DAY234_CLV
    ttt$dose_2 <- sapply(((target_AUC - ttt$pred_AUC_1) / 3) * ttt$apri_cl_2 / MWbusulfan, max, 0) #take the decrease of CL into account
  }
  
  # --- Administration of Day 2 --------------------------
  data2 <- empty_data(day = current_day, amt = ttt$dose_2, covlist = patient, sampletimes = sampletimes)
  set.seed(seed)
  sim2 <- mrgsim(sim_mrgmodel, data2)
  ttt$actual_cl_2 <- unique(sim2$CL)
  ttt$actual_AUC_2 <- MWbusulfan * ttt$dose_2 / ttt$actual_cl_2 
  
  
  # --- DAY 3 ---------------------------------------------
  current_day <- 3
  
  # --- Calculation of Dose 3 -----------------------------
  if(method3 == "curr_equation"){
    ttt$dose_3 <- curr_equation(BSA = patient$BSA, AGE = patient$AGE)
  }
  
  if(method3 == "genet_equation"){
    ttt$apri_cl_3 <- pred_CL(mod = genet_model, dat = c(patient, list(DAY = current_day))) 
    ttt$dose_3 <- ttt$apri_cl_3 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method3 == "nogenet_equation"){
    ttt$apri_cl_3 <- pred_CL(mod = nogenet_model, dat = c(patient, list(DAY = current_day))) 
    ttt$dose_3 <- ttt$apri_cl_3 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method3 == "TDM_bayes_1"){
    # No TDM at Day 2, the day before
    # We assume a priori CL and a priori exposure
    ttt$pred_AUC_2 <- MWbusulfan * ttt$dose_2 / ttt$apri_cl_2
    ttt$apri_cl_3 <- ttt$apri_cl_2 
    ttt$dose_3 <- ttt$dose_2
  }
  
  if(method3 == "TDM_bayes_2"){
    est2 <- mapbayest(est_mrgmodel, mutate(data2, DV = sim2$DV), hessian = FALSE, select_eta = c(1,2,7,8), verbose = FALSE)
    ttt$apost_cl_2 <- get_param(est2)[["CL"]]
    ttt$pred_AUC_2 <- MWbusulfan * ttt$dose_2 / ttt$apost_cl_2 
    ttt$apri_cl_3 <- ttt$apost_cl_2
    ttt$dose_3 <- sapply(((target_AUC - (ttt$pred_AUC_1 + ttt$pred_AUC_2)) / 2) * ttt$apri_cl_3 / MWbusulfan, max, 0)
  }
  
  
  # --- Administration of Day 3 -----------------------------
  data3 <- empty_data(day = current_day, amt = ttt$dose_3, covlist = patient, sampletimes = sampletimes)
  set.seed(seed)
  sim3 <- mrgsim(sim_mrgmodel, data3)
  ttt$actual_cl_3 <- unique(sim3$CL)
  ttt$actual_AUC_3 <- MWbusulfan * ttt$dose_3 / ttt$actual_cl_3 
  
  # --- DAY 4 ----------------------------------------------
  current_day <- 4
  
  # --- Calculation of Dose 4 ------------------------------
  if(method4 == "curr_equation"){
    ttt$dose_4 <- curr_equation(BSA = patient$BSA, AGE = patient$AGE)
  }
  
  if(method4 == "genet_equation"){
    ttt$apri_cl_4 <- pred_CL(mod = genet_model, dat = c(patient, list(DAY = current_day))) 
    ttt$dose_4 <- ttt$apri_cl_4 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method4 == "nogenet_equation"){
    ttt$apri_cl_4 <- pred_CL(mod = nogenet_model, dat = c(patient, list(DAY = current_day))) 
    ttt$dose_4 <- ttt$apri_cl_4 * target_AUC / (4.0 * MWbusulfan)
  }
  
  if(method4 == "TDM_bayes_1"){
    # No TDM at Day 3
    # We assume a priori CL (apri_cl_3) and a priori exposure
    ttt$pred_AUC_3 <- MWbusulfan * ttt$dose_3 / ttt$apri_cl_3
    ttt$apri_cl_4 <- ttt$apri_cl_3 
    ttt$dose_4 <- ttt$dose_2
  }
  
  if(method4 == "TDM_bayes_2"){
    # No TDM at Day 3
    # We assume a priori CL (apri_cl_3) and a priori exposure
    ttt$pred_AUC_3 <- MWbusulfan * ttt$dose_3 / ttt$apri_cl_3 
    ttt$pred_cl_4 <- ttt$apri_cl_3
    ttt$dose_4 <- ttt$dose_3
  }
  
  if(method4 == "TDM_bayes_3"){
    est3 <- mapbayest(est_mrgmodel, mutate(data3, DV = sim3$DV), hessian = FALSE, select_eta = c(1,2,9,10), verbose = FALSE)
    ttt$apost_cl_3 <- get_param(est3)[["CL"]]
    ttt$pred_AUC_3 <- MWbusulfan * ttt$dose_3 / ttt$apost_cl_3 
    ttt$apri_cl_4 <- ttt$apost_cl_3
    ttt$dose_4 <- sapply((target_AUC - (ttt$pred_AUC_1 + ttt$pred_AUC_2 + ttt$pred_AUC_3)) * ttt$apri_cl_4 / MWbusulfan, max, 0)
  }
  
  # --- Administration of Day 4 -----------------------------
  data4 <- empty_data(day = 4, amt = ttt$dose_4, covlist = patient, sampletimes = sampletimes)
  set.seed(seed)
  sim4 <- mrgsim(sim_mrgmodel, data4)
  ttt$actual_cl_4 <- unique(sim4$CL)
  ttt$actual_AUC_4 <- MWbusulfan * ttt$dose_4 / ttt$actual_cl_4 
  
  ttt$cumAUC <- ttt$actual_AUC_1 + ttt$actual_AUC_2 + ttt$actual_AUC_3 + ttt$actual_AUC_4
  
  ttt
}