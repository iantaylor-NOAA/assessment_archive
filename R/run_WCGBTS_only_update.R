#' update model inputs with new WCGBTS data
#'
#' @param common_name common name of species to pull data for, e.g. "Dover sole"
#' @param dir_model directory with model input files
#' @param wcgbts_fleet fleet number or name for WCGBTS data (if NULL, will try to auto-detect based on fleet names containing "WCGBT")
#' @return list of input files that can be written by r4ss::SS_write()
#' @example
#' inputs <- update_assessment_with_wcgbts()

update_assessment_with_wcgbts <- function(
  common_name = "Dover sole",
  # TODO: make these inputs dynamic based on species name
  dir_model = "G:/Shared drives/NMFS NWC FRAM Assessment Archive/Archives/DoverSole/DoverSole_2021/2_base_model",
  wcgbts_fleet = NULL
) {
  # read model inputs
  inputs <- r4ss::SS_read(dir_model)

  # get fleet name for WCGBTS
  if (is.null(wcgbts_fleet)) {
    wcgbts_fleet <- which(grepl("WCGBT", inputs$dat$fleetnames))
    if (length(wcgbts_fleet) == 0) {
      cli::cli_abort(
        "No fleet mating 'WCGBT', please input a fleet name or number `wcgbts_fleet`"
      )
    }
  }
  cli::cli_alert_info("Updating fleet {wcgbts_fleet} with new WCGBTS data")

  if (is.character(wcgbts_fleet)) {
    wcgbts_fleet_input <- wcgbts_fleet
    wcgbts_fleet <- which(inputs$dat$fleetnames == wcgbts_fleet)
    if (length(wcgbts_fleet) == 0) {
      cli::cli_abort(
        "No fleet mating {wcgbts_fleet_input}. Fleetnames: {inputs$dat$fleetnames}"
      )
    }
    cli::cli_alert_info(
      "{inputs$dat$fleetnames[wcgbts_fleet]} is fleet number {wcgbts_fleet}"
    )
  }

  # check month of existing data
  wcgbts_month <- inputs$dat$CPUE |>
    dplyr::filter(index == wcgbts_fleet) |>
    dplyr::pull(month) |>
    as.numeric() |>
    unique()

  if (length(wcgbts_month) > 1) {
    cli::cli_abort(
      "Multiple months found for fleet {wcgbts_fleet}: {wcgbts_month}, so code adjustment is needed to handle this case."
    )
  }
  cli::cli_alert_info(
    "Month of existing data for WCGBTS index is {wcgbts_month}"
  )

  # get new survey index data
  index <- get_est_by_area_csv(common_name)
  index <- index |>
    dplyr::filter(area == "Coastwide") |>
    dplyr::select(year, est, se) |>
    dplyr::mutate(index = wcgbts_fleet, month = wcgbts_month) |>
    dplyr::rename(year = year, obs = est, se_log = se) |>
    dplyr::select(year, month, index, obs, se_log)

  # get new survey bio data
  # commands are all from https://pfmc-assessments.github.io/nwfscSurvey/articles/nwfscSurvey.html
  catch <- nwfscSurvey::pull_catch(common_name = common_name)
  bio <- nwfscSurvey::pull_bio(common_name = common_name)

  # TODO: make strata dynamic based on species
  strata <- nwfscSurvey::create_strata(
    names = paste(
      sep = "_",
      rep(times = 4, c("55m-183m", "183m-549m", "549m-900m", "900m-1280m")),
      rep(each = 4, c("32-34.5", "34.5-42.0", "42.0-46.0", "46.0-49"))
    ),
    depths_shallow = rep(times = 4, x = c(55, 183, 549, 900)),
    depths_deep = rep(times = 4, x = c(183, 549, 900, 1280)),
    lats_south = rep(each = 4, x = c(32, 34.5, 42, 46.0)),
    lats_north = rep(each = 4, x = c(34.5, 42.0, 46.0, 49))
  )

  # get length data bins
  length_bins <- inputs$dat$lbin_vector
  age_data_bins <- inputs$dat$agebin_vector

  # get ageing error matrix for WCGBTS ages
  wcgbts_ageerr <- inputs$dat$agecomp |>
    dplyr::filter(fleet == wcgbts_fleet) |>
    dplyr::pull(ageerr) |>
    unique()
  if (length(wcgbts_ageerr) > 1) {
    cli::cli_abort(
      "Multiple age error matrices found for fleet {wcgbts_fleet}: {wcgbts_ageerr}, so code adjustment is needed to handle this case."
    )
  }

  # figure out marginal vs conditional age comps
  fleets_with_marginal_ages <- inputs$dat$agecomp |>
    dplyr::filter(Lbin_lo == -1) |>
    dplyr::pull(fleet) |>
    unique()
  fleets_with_caal <- inputs$dat$agecomp |>
    dplyr::filter(Lbin_lo > 0) |>
    dplyr::pull(fleet) |>
    unique()

  length_comps <- nwfscSurvey::get_expanded_comps(
    bio_data = bio,
    catch_data = catch,
    comp_bins = length_bins,
    strata = strata,
    comp_column_name = "length_cm",
    output = "full_expansion_ss3_format",
    two_sex_comps = TRUE,
    input_n_method = "stewart_hamel"
  )$sexed |>
    dplyr::rename(part = partition, Nsamp = input_n) |>
    dplyr::mutate(month = wcgbts_month, fleet = wcgbts_fleet)

  age_comps <- nwfscSurvey::get_expanded_comps(
    bio_data = bio,
    catch_data = catch,
    comp_bins = age_data_bins,
    strata = strata,
    comp_column_name = "age",
    output = "full_expansion_ss3_format",
    two_sex_comps = TRUE,
    input_n_method = "stewart_hamel"
  )$sexed |>
    dplyr::rename(part = partition, Nsamp = input_n) |>
    dplyr::mutate(
      ageerr = wcgbts_ageerr,
      month = wcgbts_month,
      fleet = wcgbts_fleet
    )
  # exclude marginal data using negative fleet number if needed
  if (-wcgbts_fleet %in% fleets_with_marginal_ages) {
    age_comps <- age_comps |> dplyr::mutate(fleet = -1 * wcgbts_fleet)
    cli::cli_alert_info(
      "Marginal age comps for WCGBTS are excluded from the likelihood, so they will be marked with fleet number {-1 * wcgbts_fleet}."
    )
  }

  caal <- nwfscSurvey::get_raw_caal(
    data = bio,
    len_bins = length_bins,
    age_bins = age_data_bins
  ) |>
    dplyr::rename(part = partition, Nsamp = input_n) |>
    dplyr::mutate(
      ageerr = wcgbts_ageerr,
      month = wcgbts_month,
      fleet = wcgbts_fleet
    )
  # exclude caal data using negative fleet number if needed
  if (-wcgbts_fleet %in% fleets_with_caal) {
    caal <- caal |> dplyr::mutate(fleet = -wcgbts_fleet)
    cli::cli_alert_info(
      "Conditional age-at-length comps for WCGBTS are excluded from the likelihood, so they will be marked with fleet number {-1 * wcgbts_fleet}."
    )
  }

  # copy inputs
  inputs2 <- inputs

  # update index data
  inputs2$dat$CPUE <- inputs$dat$CPUE |>
    dplyr::filter_out(index == wcgbts_fleet) |>
    dplyr::bind_rows(index)

  # update length comp data
  inputs2$dat$lencomp <- inputs$dat$lencomp |>
    dplyr::mutate(month = as.numeric(month)) |>
    dplyr::filter_out(abs(fleet) == wcgbts_fleet) |>
    dplyr::bind_rows(length_comps |> dplyr::mutate(month = as.numeric(month)))

  # update age comp data
  inputs2$dat$agecomp <- inputs$dat$agecomp |>
    dplyr::mutate(month = as.numeric(month)) |>
    dplyr::filter_out(abs(fleet) == wcgbts_fleet) |>
    dplyr::bind_rows(age_comps, caal)

  # get catches from GEMM
  gemm <- nwfscSurvey::pull_gemm(common_name = common_name)

  # TODO: consider adding final year of catch as recent average
  
  # TODO: check for discards in the model and deal with those somehow

  # simplistic assumption of catch allocation among fleets constant over time
  # TODO: replace this with a system to assign GEMM sectors to fleet
  # and/or use pacfin to allow spatial separation of catch
  catch_prop <- inputs$dat$catch |>
    dplyr::filter(year > inputs$dat$endyr - 10) |> # use recent 10 years for average proportion
    dplyr::group_by(fleet) |>
    dplyr::summarise(catch = sum(catch)) |>
    dplyr::mutate(prop = catch / sum(catch)) |>
    dplyr::select(fleet, prop)
  gemm_catch_by_year <- gemm |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      catch_total = sum(total_discard_with_mort_rates_applied_and_landings_mt)
    )
  gemm_catch_by_year_and_fleet <- gemm_catch_by_year |>
    dplyr::cross_join(catch_prop) |>
    dplyr::mutate(catch = catch_total * prop) |>
    dplyr::mutate(seas = 1, catch_se = min(inputs$dat$catch$catch_se)) |>
    dplyr::select(year, seas, fleet, catch, catch_se)

  # check for similarity between fleet-specific catches in overlapping years
  # between GEMM and model inputs
  catch_comparison <- inputs$dat$catch |>
    dplyr::filter(year > inputs$dat$endyr - 10) |>
    dplyr::inner_join(
      gemm_catch_by_year_and_fleet,
      by = c("year", "seas", "fleet"),
      suffix = c("_model", "_gemm_avg_split")
    ) |>
    dplyr::mutate(
      catch_model = round(catch_model, 0),
      catch_gemm_avg_split = round(catch_gemm_avg_split, 0),
      catch_ratio = catch_model / catch_gemm_avg_split,
      catch_ratio_diff = round(catch_ratio - 1.0, 2)
    ) |>
    dplyr::select(
      year,
      fleet,
      catch_model,
      catch_gemm_avg_split,
      catch_ratio_diff
    ) |>
    dplyr::arrange(year, fleet)

  cli::cli_alert_info(
    "Comparison of model catch and GEMM catch in recent years:"
  )
  print(catch_comparison |> knitr::kable())

  # append catch data for years from endyr + 1 to the final year of the gemm catch data
  if (any(inputs$dat$catch$year > inputs$dat$endyr)) {
    cli::cli_alert_warning(
      "Model inputs already contain catch data for years after endyr, adjust code as needed."
    )
  }
  inputs2$dat$catch <- inputs$dat$catch |>
    dplyr::bind_rows(
      gemm_catch_by_year_and_fleet |>
        dplyr::filter(year > inputs$dat$endyr)
    )
  # update end year
  inputs2$dat$endyr <- max(inputs2$dat$catch$year)

  # change starter file to stop using the .par file (if not already set this way)
  inputs2$start$init_values_src <- 0
  # change output detail
  inputs2$start$run_display_detail <- 0

  return(inputs2)
}

if (FALSE) {
  # run model
  r4ss:::SS_write(inputs2, "test", overwrite = TRUE)
  r4ss:::run(
    "test",
    extras = "-nohess",
    show_in_console = TRUE,
    skipfinished = FALSE
  )

  # compare results
  mod1 <- r4ss::SS_output(dir_model)
  mod2 <- r4ss::SS_output("test")
  r4ss::SSplotComparisons(
    SSsummarize(
      list(mod1, mod2),
    ),
    legendlabels = c(
      glue::glue("Original (end year {inputs$dat$endyr})"),
      glue::glue("Updated with new WCGBTS data (end year {inputs2$dat$endyr})")
    ),
    indexPlotEach = TRUE,
    endyrvec = c(inputs$dat$endyr, inputs2$dat$endyr)
  )
}
# lessons learned:
# - recent catch by sector is easily available from GEMM, but not set up for spatial divisions (e.g. CA vs OR+WA in the Dover assessment)
#
