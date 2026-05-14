#' run retrospective analysis removing most recent 10 years of data except for the WCGBTS index
#' 
#' @param species character, name of species to run retrospective for
#' @param dir character, directory of original model files
#' @param newdir character, directory to write new model files to
#' @param inputs list, model input list as read by r4ss::SS_read, if NULL will be read from dir
#' @param WCGBTS_fleet numeric, fleet number of the WCGBTS index

index_only_retro <- function(
  species,
  dir = NULL,
  newdir = NULL,
  inputs = NULL,
  WCGBTS_fleet = NULL
) {
  if (is.null(inputs)) {
    inputs <- r4ss::SS_read(dir)
  }
  endyr <- inputs$dat$endyr
  last_data_year <- endyr - 10
  dat <- inputs$dat
  newdat <- inputs$dat
  
  # remove all data from the past 10 years except for the data from the WCGBTS 
  newdat$lencomp <- dat$lencomp |>
    dplyr::filter(year <= last_data_year | fleet == WCGBTS_fleet)
  newdat$agecomp <- dat$agecomp |>
    dplyr::filter(year <= last_data_year | fleet == WCGBTS_fleet)
  newdat$CPUE <- dat$CPUE |>
    dplyr::filter(year <= last_data_year | index == WCGBTS_fleet)

  # replace the data in the model input with the new data
  inputs$dat <- newdat

  # recude run display detail
  inputs$start$run_display_detail <- 0

  # write the new data to a new directory
  inputs$dir <- newdir
  r4ss::SS_write(inputs, dir = inputs$dir, overwrite = TRUE)

  return(inputs$dir)
}