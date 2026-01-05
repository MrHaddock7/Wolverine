###This script is used to automatically knit together all reports for the zoos.

library(rmarkdown)
library(tibble)
library(dplyr)
library(purrr)
library(stringr)

zoos <- tibble(zoo = c(
  "Skansen",
  "Helsinki",
  "Gaia",
  "Herberstein",
  "Osnabrück",
  "Duisburg",
  "Salzburg",
  "Nordens ark",
  "Lycksele",
  "Vildriket",
  "Kolmården",
  "Borås",
  "Ähtäri"
)) %>%
  pull(zoo)

zoos_pdf <- tibble(zoo = c(
  "Skansen",
  "Helsinki",
  "Gaia",
  "Herberstein",
  "Osnabrück",
  "Duisburg",
  "Salzburg",
  "Nordens ark",
  "Lycksele",
  "Vildriket",
  "Kolmarden",
  "Boras",
  "Ahtari"
)) %>%
  pull(zoo)


zoos

reports <- tibble(
  input = "C:/Users/Lovisa/Documents/Wolverine/scripts/report_generation/testing_pdf.Rmd",
  output_file = str_glue("{zoos_pdf}.pdf"),
  params = map(zoos, ~ list(zoo = .))
)

pwalk(reports, render)
