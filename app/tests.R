dose_data <- data.frame(
  DAY = 0:4,
  Date = paste0(11:15, "/10/2023"),
  Dose = c(0, 100, 100, 100, 0),
  Hour_start = rep(NA_character_, 5),
  Hour_end = rep(NA_character_, 5)
)

conc_data <- data.frame(
  DAY = rep(c(1L,3L), each = 3),
  Date = rep(c("12/10/2023", "14/10/2023"), each = 3),
  Hour_sample = c("3:00" ,rep(NA_character_, 5)),
  Concentration = c(3000, rep(NA, 5)),
  Exclude = rep(FALSE, 6)
)

apri <- apriori(bw = 45, malign = TRUE, gsta1 = TRUE, auc = 18000, ndays = 5)
mooodel <- mread("app/busulfan_shiny.cpp")
daaata <- make_data(dose_data, conc_data, bw = 25)
eeeest <- mapbayest(mooodel, daaata)
new_post <- make_post(dose_data, conc_data, eeeest, auc = 16000)
plooot <- make_plot(est = eeeest, post = new_post)

rmarkdown::render("app/busulfan_report.Rmd", params = list(
  all_inputs = list(BW = 26, gsta1 = "", ID = "123456789"),
  re_apriori = apri,
  concentrationData = conc_data, 
  dosingData = dose_data,
  nmtranData = daaata,
  re_post = new_post,
  re_plot = plooot
))

mooodel %>% 
  adm_rows(time = 0, rate = 85/3, amt = 85, cmt = 1) %>% 
  obs_rows(
    time = c(180, 195, 210, 240, 300, 360, 480)/60, 
    DV = c(5530, 5240, 4350, 3890, 2710, 2370, 1170), 
    cmt = 1
  ) %>% 
  mapbayest()


mooodel %>% 
  adm_rows(time = 0, rate = 85/3, amt = 85, cmt = 1, 18.35) %>% 
  obs_rows(
    time = c(195, 300, 480)/60, 
    DV = c(5240, 2710, 1170), 
    cmt = 1
  ) %>% 
  mapbayest()

dose_data2 <- data.frame(
  DAY = 0:4,
  Date = paste0(11:15, "/10/2023"),
  Dose = c(0, 85, 0, 0, 0),
  Hour_start = rep(NA_character_, 5),
  Hour_end = rep(NA_character_, 5)
)

conc_data2 <- data.frame(
  DAY = rep(c(1L,3L), each = 3),
  Date = rep(c("12/10/2023", "14/10/2023"), each = 3),
  Hour_sample = c("3:25" , "5h00", "8h00", rep(NA_character_, 3)),
  Concentration = c(c(5240, 2710, 1170), rep(NA, 3)),
  Exclude = rep(FALSE, 6)
)

data2 <-  make_data(dose_data2, conc_data2, bw = 18.35)
eeeest2 <- mapbayest(mooodel, data2)
new_post <- make_post(dose_data2, conc_data2, eeeest2, auc = 16000)
# should find AUC around 6500

