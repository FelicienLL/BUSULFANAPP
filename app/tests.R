dose_data <- data.frame(
  DAY = 0:4,
  Date = lubridate::as_date(seq(from = Sys.Date(), by = 1, length.out = 5)),
  Dose = c(0, 100, 100, 100, 0),
  Hour_start = rep(NA_character_, 5),
  Hour_end = rep(NA_character_, 5)
)

conc_data <- data.frame(
  DAY = rep(c(1L,3L), each = 3),
  Date = rep(c(Sys.Date() + 1, Sys.Date() + 3), each = 3),
  Hour_sample = c("3:00" ,rep(NA_character_, 5)),
  Concentration = c(3000, rep(NA, 5)),
  Exclude = rep(FALSE, 6)
)

mooodel <- mread("app/busulfan_shiny.cpp")
daaata <- make_data(dose_data, conc_data, bw = 25)
eeeest <- mapbayest(mooodel, daaata)
new_post <- make_post(dose_data, conc_data, eeeest, auc = 16000)

