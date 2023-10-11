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
