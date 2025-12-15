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
  "Ústí nad Labem",
  "Nordens ark",
  "Lycksele",
  "Vildriket",
  "Kolmården",
  "Borås",
  "Ähtäri"
)) %>%
  pull(zoo)

zoos

reports <- tibble(
  input = "C:/Users/Lovisa/Documents/Wolverine/scripts/report_generation/testing_pdf.Rmd",
  output_file = str_glue("{zoos}.pdf"),
  params = map(zoos, ~ list(zoo = .))
)

pwalk(reports, render)
