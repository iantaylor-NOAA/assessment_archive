# script to run retrospective analysis removing most recent 10 years of data except for the WCGBTS index

# switches for which things to run if sourcing this script
make_tables <- FALSE
modify_model_inputs <- FALSE
run_models <- FALSE
load_results_text <- FALSE
load_results_rda <- TRUE
make_plots_tables <- TRUE


if (make_tables) {
  cli::cli_alert_info("sourcing functions and loading model inputs")
  # run script which makes the info2_recent table with info on recent assessments
  source("R/load_assessments.R")
  # source function to modify model input files for the retrospective
  source("R/remove_recent_data.R")

  cli::cli_alert_info("making info3_wcgbts table")
  # identify the WCGBTS fleet number for each model
  for (i in 1:nrow(info2_recent)) {
    WCGBTS_fleet <- which(grepl("WCGBT", inputs_recent[[i]]$dat$fleetnames))
    if (length(WCGBTS_fleet) == 1) {
      info2_recent$WCGBTS_number[i] <- WCGBTS_fleet
      info2_recent$WCGBTS_name[i] <- inputs_recent[[i]]$dat$fleetnames[
        WCGBTS_fleet
      ]
    } else {
      WCGBTS_fleet <- which(grepl(
        "coastwide_NWFSC", # this is a subset of the name used by Canary
        inputs_recent[[i]]$dat$fleetnames
      ))
      if (length(WCGBTS_fleet) == 1) {
        info2_recent$WCGBTS_number[i] <- WCGBTS_fleet
        info2_recent$WCGBTS_name[i] <- inputs_recent[[i]]$dat$fleetnames[
          WCGBTS_fleet
        ]
      } else {
        info2_recent$WCGBTS_number[i] <- NA
        info2_recent$WCGBTS_name[i] <- NA
      }
    }
  }

  # black, copper, and quillback all don't use the WCGBTS, so skip them
  info3_wcgbts <- info2_recent |>
    dplyr::filter_out(
      Species %in% c("Black rockfish", "Copper rockfish", "Quillback rockfish")
    )

  # directory where retro will go
  info3_wcgbts$wcgbts_only_retro_dir <- file.path(
    "data-raw/wcgbts_only_retro",
    info3_wcgbts$Species
  )

  info3_wcgbts$Species |> data.frame()
  #   `info3_wcgbts$Species`
  #   <chr>
  # 1 Dover sole
  # 2 Pacific spiny dogfish
  # 3 Canary rockfish
  # 4 Petrale sole
  # 5 Yellowtail rockfish
  # 6 Rougheye/blackspotted rockfish
  # 7 Sablefish
}

if (modify_model_inputs) {
  cli::cli_alert_info("modifying input files")
  dir.create(
    "data-raw/wcgbts_only_retro",
    recursive = TRUE,
    showWarnings = FALSE
  )

  # modify models to remove 10 years from all data sources other than the WCGBTS
  for (i in 1:nrow(info3_wcgbts)) {
    cli::cli_inform("Writing input files for {info3_wcgbts$Species[i]}...")
    dir <- remove_recent_data(
      species = info3_wcgbts$Species[i],
      dir = info3_wcgbts$Model.archive.directory[i],
      newdir = info3_wcgbts$wcgbts_only_retro_dir[i],
      WCGBTS_fleet = info3_wcgbts$WCGBTS_number[i]
    )
  }
}

if (run_models) {
  cli::cli_alert_info("running models...")

  for (i in 1:nrow(info3_wcgbts)) {
    cli::cli_inform(
      "Running index-only retrospective for {info3_wcgbts$Species[i]}..."
    )
    r4ss::run(
      dir = info3_wcgbts$wcgbts_only_retro_dir[i],
      show_in_console = TRUE
    )
  }

  # problems with a few models:
  # * Sablefish needs change to Q param section for Recruitment_Index
  # * Dover has problem reading ss.par
  # * Canary has some selectivity blocks in starting after the removed years of data, causing bad estimates of selectivity (which impacts the OFLs). Need to modify block design and remove parameters for the blocks that are being removed.

  # fix sablefish
  inputs <- r4ss::SS_read("data-raw/wcgbts_only_retro/Sablefish")
  inputs$start$init_values_src <- 0 # don't use par because the number of parameters has changed
  inputs$ctl$Q_options <- inputs$ctl$Q_options |> dplyr::filter_out(fleet == 11)
  inputs$ctl$Q_parms <- inputs$ctl$Q_parms |>
    dplyr::filter_out(grepl("Recruitment_Index", rownames(inputs$ctl$Q_parms)))
  r4ss::SS_write(
    inputs,
    dir = "data-raw/wcgbts_only_retro/Sablefish",
    overwrite = TRUE
  )
  r4ss::run(
    dir = "data-raw/wcgbts_only_retro/Sablefish",
    show_in_console = TRUE
  )

  # fix dover (the ss.par file from the 2021 SS3 version included extra unused parameters for implementation error)
  inputs <- r4ss::SS_read("data-raw/wcgbts_only_retro/Dover sole")
  inputs$start$init_values_src <- 0 # don't use par
  r4ss::SS_write(
    inputs,
    dir = "data-raw/wcgbts_only_retro/Dover sole",
    overwrite = TRUE
  )
  r4ss::run(
    dir = "data-raw/wcgbts_only_retro/Dover sole",
    show_in_console = TRUE
  )

  # fix canary
  inputs <- r4ss::SS_read("data-raw/wcgbts_only_retro/Canary rockfish")
  # figure out which fleets have selectivity blocks
  inputs$ctl$size_selex_parms |>
    dplyr::filter(Block > 0) |>
    rownames() |>
    # extract fleet number, e.g. 29 from the parentheses at the end of the rowname
    # which is in the format "SizeSel_PFemOff_4_29_coastwide_Tri_early(29)"
    stringr::str_extract("\\(\\d+\\)") |>
    stringr::str_replace("\\)", "") |>
    stringr::str_replace("\\(", "") |>
    as.numeric() |>
    unique()
  # [1] 1 2 4 5 7 8 9

  # which blocks end after endyr - 10?
  maxblockyr <- inputs$dat$endyr - 10
  # which block patterns include blocks which start after maxblockyr?
  inputs$ctl$Block_Design
  # [[1]]
  # [1] 2000 2010 2011 2022

  # [[2]]
  # [1] 2000 2019 2020 2022

  # [[3]]
  # [1] 2004 2016 2017 2022

  # [[4]]
  # [1] 2004 2014 2015 2022

  # [[5]]
  # [1] 2006 2020

  # [[6]]
  # [1] 2000 2019

  # which vectors have the 2nd from last element > maxblockyr?
  (late_block_patterns <- inputs$ctl$Block_Design |>
    sapply(function(x) x[length(x) - 1] > maxblockyr) |>
    which())
  # [1] 2 3 4

  # which fleet's selectivity parameteters use those blocks?
  (fleets_with_late_blocks <- inputs$ctl$size_selex_parms |>
    dplyr::filter(Block %in% late_block_patterns) |>
    rownames() |>
    # extract fleet number, e.g. 29 from the parentheses at the end of the rowname
    # which is in the format "SizeSel_PFemOff_4_29_coastwide_Tri_early(29)"
    stringr::str_extract("\\(\\d+\\)") |>
    stringr::str_replace("\\)", "") |>
    stringr::str_replace("\\(", "") |>
    as.numeric() |>
    unique())
  # [1] 4 7 8

  # check fleet names to confirm they are fishing fleets and not indices
  inputs$dat$fleetnames[fleets_with_late_blocks]
  # [1] "4_CA_NTWL" "7_CA_REC"  "8_OR_REC"

  # combine the final block with the previous block from each block pattern
  new_Block_Design <- inputs$ctl$Block_Design
  for (iblock in 1:length(inputs$ctl$Block_Design)) {
    if (iblock %in% late_block_patterns) {
      if (length(inputs$ctl$Block_Design[[iblock]]) == 2) {
        cli::cli_alert_danger(
          "Block pattern {iblock} only has 1 blocks, can't modify"
        )
      }
      new_Block_Design[[iblock]] <- inputs$ctl$Block_Design[[iblock]][
        -(length(inputs$ctl$Block_Design[[iblock]]) - 1:2)
      ]
    }
  }
  # replace old with new block design
  inputs$ctl$Block_Design <- new_Block_Design
  inputs$ctl$blocks_per_pattern <- sapply(inputs$ctl$Block_Design, length) / 2

  # remove extra parameter lines for the blocks that were removed
  new_size_selex_parms_tv <- inputs$ctl$size_selex_parms_tv
  # filter out parameters for which the rowname ends in a year greater than maxblockyr
  new_size_selex_parms_tv <- new_size_selex_parms_tv |>
    dplyr::filter_out(
      rownames(new_size_selex_parms_tv) |> # extract year from rowname, e.g. "SizeSel_PFemOff_4_29_coastwide_Tri_early(29)_2017"
        stringr::str_extract("\\d{4}$") |> # extract the 4 digit year at the end of the rowname
        as.numeric() >
        maxblockyr
    )
  #  [1] "SizeSel_P_1_4_CA_NTWL(4)_BLK2repl_2020"       "SizeSel_P_3_4_CA_NTWL(4)_BLK2repl_2020"       "SizeSel_P_4_4_CA_NTWL(4)_BLK2repl_2020"
  #  [4] "SizeSel_PFemOff_3_4_CA_NTWL(4)_BLK2repl_2020" "SizeSel_P_1_7_CA_REC(7)_BLK3repl_2017"        "SizeSel_P_3_7_CA_REC(7)_BLK3repl_2017"
  #  [7] "SizeSel_P_4_7_CA_REC(7)_BLK3repl_2017"        "SizeSel_PFemOff_3_7_CA_REC(7)_BLK3repl_2017"  "SizeSel_P_1_8_OR_REC(8)_BLK4repl_2015"
  # [10] "SizeSel_P_3_8_OR_REC(8)_BLK4repl_2015"        "SizeSel_P_4_8_OR_REC(8)_BLK4repl_2015"        "SizeSel_PFemOff_3_8_OR_REC(8)_BLK4repl_2015"

  # which parameters were removed?
  setdiff(
    rownames(inputs$ctl$size_selex_parms_tv),
    rownames(new_size_selex_parms_tv)
  )
  # replace old with new time-varying parameters
  inputs$ctl$size_selex_parms_tv <- new_size_selex_parms_tv

  # write modified inputs
  r4ss::SS_write(
    inputs,
    dir = "data-raw/wcgbts_only_retro/Canary rockfish_revised_blocks",
    overwrite = TRUE
  )
  # run model
  r4ss::run(
    dir = "data-raw/wcgbts_only_retro/Canary rockfish_revised_blocks",
    show_in_console = TRUE
  )
}

if (make_tables) {
  cli::cli_alert_info("expanding info3_wcgbts table...")
  # change the directory in the info table
  info3_wcgbts$wcgbts_only_retro_dir[
    info3_wcgbts$Species == "Canary rockfish"
  ] <-
    "data-raw/wcgbts_only_retro/Canary rockfish_revised_blocks"
  # done fixing canary

  # additional analyses:
  # - removing WCGBTS ages but keeping lengths and indices (if we haven't prioritized ageing for some species)
  # - removing all WCGBTS data from recent years as well (like a catch-only projection),

  # first make directories for these additional analyses
  info3_wcgbts <- info3_wcgbts |>
    dplyr::mutate(
      wcgbts_no_ages_dir = stringr::str_replace(
        wcgbts_only_retro_dir,
        "wcgbts_only_retro",
        "wcgbts_no_ages"
      ),
      catch_only_dir = stringr::str_replace(
        wcgbts_only_retro_dir,
        "wcgbts_only_retro",
        "catch_only"
      )
    )

  # additional Canary models were already run when the code above was written, so need to restore their directories
  info3_wcgbts$wcgbts_no_ages_dir[
    info3_wcgbts$Species == "Canary rockfish"
  ] <-
    "data-raw/wcgbts_no_ages/Canary rockfish"
  # additional Canary models were already run when the code above was written, so need to restore their directories
  info3_wcgbts$catch_only_dir[
    info3_wcgbts$Species == "Canary rockfish"
  ] <-
    "data-raw/catch_only/Canary rockfish"
  # done fixing canary
}

if (modify_model_inputs) {
  cli::cli_alert_info("modifying model inputs for additional analyses...")
  # remove WCGBTS ages but keep lengths and indices (starting from the original retro directory)
  for (i in 1:nrow(info3_wcgbts)) {
    cli::cli_inform("Writing input files for {info3_wcgbts$Species[i]}...")
    remove_recent_data(
      species = info3_wcgbts$Species[i],
      dir = info3_wcgbts$wcgbts_only_retro_dir[i],
      newdir = info3_wcgbts$wcgbts_no_ages_dir[i],
      WCGBTS_fleet = info3_wcgbts$WCGBTS_number[i],
      remove_WCGBTS_ages = TRUE
    )
  }

  # remove all WCGBTS data from recent years as well (like a catch-only projection), starting from the original retro directory
  for (i in 1:nrow(info3_wcgbts)) {
    cli::cli_inform("Writing input files for {info3_wcgbts$Species[i]}...")
    dir <- remove_recent_data(
      species = info3_wcgbts$Species[i],
      dir = info3_wcgbts$wcgbts_only_retro_dir[i],
      newdir = info3_wcgbts$catch_only_dir[i],
      WCGBTS_fleet = info3_wcgbts$WCGBTS_number[i],
      remove_WCGBTS_ages = TRUE,
      remove_WCGBTS_lengths = TRUE,
      remove_WCGBTS_index = TRUE
    )
  }
}

if (run_models) {
  cli::cli_alert_info("running more models...")
  # run models for these additional analyses
  # run in parallel

  # set up parallel plan with 7 workers (one for each species)
  future::plan(future::multisession, workers = 7)
  furrr::future_map(
    .x = info3_wcgbts$wcgbts_no_ages_dir,
    .f = r4ss::run,
    skipfinished = FALSE,
    show_in_console = FALSE #,
    # extras = "-nohess -stopph 0"
  )
  # remove parallel plan
  future::plan(future::sequential)

  future::plan(future::multisession, workers = 7)
  furrr::future_map(
    .x = info3_wcgbts$catch_only_dir,
    .f = r4ss::run,
    skipfinished = FALSE,
    show_in_console = FALSE #,
    # extras = "-nohess -stopph 0"
  )
  # remove parallel plan
  future::plan(future::sequential)
}


### loading and comparing results
if (load_results_text) {
  cli::cli_alert_info("loading model output from text files ...")
  # load original models
  outputs_original <- r4ss::SSgetoutput(
    dirvec = info3_wcgbts$Model.archive.directory
  )

  # load wcgbts_only_retro models
  output_retro <- r4ss::SSgetoutput(
    dirvec = info3_wcgbts$wcgbts_only_retro_dir
  )
  # load wcgbts_no_ages models
  output_no_ages <- r4ss::SSgetoutput(
    dirvec = info3_wcgbts$wcgbts_no_ages_dir
  )
  # load catch-only models
  output_catch_only <- r4ss::SSgetoutput(
    dirvec = info3_wcgbts$catch_only_dir
  )

  save(
    info3_wcgbts,
    inputs_recent,
    outputs_original,
    output_retro,
    output_no_ages,
    output_catch_only,
    file = "data-raw/wcgbts_only_retro/saved_inputs_outputs.rda"
  )
}

if (load_results_rda) {
  cli::cli_alert_info("loading model output from rda files...")
  load(file = "data-raw/wcgbts_only_retro/saved_inputs_outputs.rda")
}

if (make_plots_tables) {
  cli::cli_alert_info("making plots and tables ...")

  dir.create(
    "figures/wcgbts_only_retro",
    recursive = TRUE,
    showWarnings = FALSE
  )

  make_plot <- function(
    filename,
    subplots,
    xlim = NULL,
    order = 1:nrow(info3_wcgbts)
  ) {
    filename = file.path("figures/wcgbts_only_retro", filename)
    png(
      filename,
      width = 10,
      height = 6,
      units = "in",
      res = 300
    )
    par(mfrow = c(2, 4), mar = c(2, 3, 2, .5), oma = c(0, 2, 0, 0))
    for (i in (1:nrow(info3_wcgbts))[order]) {
      r4ss::SSplotComparisons(
        r4ss::SSsummarize(
          list(
            outputs_original[[i]],
            output_retro[[i]],
            output_no_ages[[i]],
            output_catch_only[[i]]
          ),
          verbose = FALSE
        ),
        legendlabels = c(
          "Original full assessment",
          "Catch + WCGBTS (ind + len + age)",
          "Catch + WCGBTS (ind + len)",
          "Catch-only"
        ),
        new = FALSE,
        subplots = subplots,
        indexfleets = info3_wcgbts$WCGBTS_number[i],
        indexQlabel = TRUE,
        xlim = xlim,
        legend = (i == 1),
        legendloc = "bottomleft"
      )
      title(info3_wcgbts$Species[i])
    }
    dev.off()
  }

  # new table with quantities of interest
  frac_unfished_colname <- paste0("Bratio_", info3_wcgbts$Year + 1)
  OFL_colname <- paste0("OFLCatch_", info3_wcgbts$Year + 3)

  # make empty tables
  tab_frac <- info3_wcgbts |>
    dplyr::select(Species) |>
    dplyr::mutate(
      frac_unfished_original = NA,
      frac_unfished_wcgbts_only_retro = NA,
      frac_unfished_wcgbts_no_ages = NA,
      frac_unfished_catch_only = NA
    )
  tab_ofl <- info3_wcgbts |>
    dplyr::select(Species) |>
    dplyr::mutate(
      OFL_original = NA,
      OFL_wcgbts_only_retro = NA,
      OFL_wcgbts_only_retro = NA,
      OFL_wcgbts_no_ages = NA
    )

  # fill in tables
  for (i in 1:nrow(info3_wcgbts)) {
    # get fraction unfished
    tab_frac$frac_unfished_original[i] <-
      outputs_original[[i]]$derived_quants[frac_unfished_colname[i], "Value"] |>
      round(3)
    tab_frac$frac_unfished_wcgbts_only_retro[i] <- output_retro[[
      i
    ]]$derived_quants[
      frac_unfished_colname[i],
      "Value"
    ] |>
      round(3)
    tab_frac$frac_unfished_wcgbts_no_ages[i] <-
      output_no_ages[[i]]$derived_quants[frac_unfished_colname[i], "Value"] |>
      round(3)
    tab_frac$frac_unfished_catch_only[i] <-
      output_catch_only[[i]]$derived_quants[
        frac_unfished_colname[i],
        "Value"
      ] |>
      round(3)
    # get OFLs
    tab_ofl$OFL_original[i] <-
      outputs_original[[i]]$derived_quants[OFL_colname[i], "Value"] |>
      round(0)
    tab_ofl$OFL_wcgbts_only_retro[i] <- output_retro[[i]]$derived_quants[
      OFL_colname[i],
      "Value"
    ] |>
      round(0)
    tab_ofl$OFL_wcgbts_no_ages[i] <-
      output_no_ages[[i]]$derived_quants[OFL_colname[i], "Value"] |>
      round(3)
    tab_ofl$OFL_catch_only[i] <-
      output_catch_only[[i]]$derived_quants[OFL_colname[i], "Value"] |>
      round(3)
  }

  # calculate ratios
  tab_frac <- tab_frac |>
    dplyr::mutate(
      ratio_wcgbts_only_retro = round(
        frac_unfished_wcgbts_only_retro / frac_unfished_original,
        2
      ),
      ratio_wcgbts_no_ages = round(
        frac_unfished_wcgbts_no_ages / frac_unfished_original,
        2
      ),
      ratio_catch_only = round(
        frac_unfished_catch_only / frac_unfished_original,
        2
      )
    )
  tab_ofl <- tab_ofl |>
    dplyr::mutate(
      OFL_ratio_wcgbts_only_retro = round(
        OFL_wcgbts_only_retro / OFL_original,
        2
      ),
      OFL_ratio_wcgbts_no_ages = round(OFL_wcgbts_no_ages / OFL_original, 2),
      OFL_ratio_catch_only = round(OFL_catch_only / OFL_original, 2)
    )

  OFL_ratio_minus_one = abs(tab_ofl$OFL_ratio_wcgbts_only_retro - 1)
  order <- order(OFL_ratio_minus_one, decreasing = FALSE)
  # [1] 1 4 2 7 5 3 6

  tab_frac <- tab_frac |>
    dplyr::slice(order)
  tab_ofl <- tab_ofl |>
    dplyr::slice(order)

  # write tables
  write.csv(
    tab_frac,
    "reports/wcgbts_only_retro/retro_results_fraction_unfished.csv",
    row.names = FALSE
  )
  write.csv(
    tab_ofl,
    "reports/wcgbts_only_retro/retro_results_OFL.csv",
    row.names = FALSE
  )

  # make plots for all years
  make_plot("fit_to_WCGBTS.png", subplots = 13, order = order)
  make_plot("spawning_output.png", subplots = 2, order = order)
  make_plot("fraction_unfished.png", subplots = 4, order = order)
  make_plot("summary_biomass.png", subplots = 18, order = order)

  # plots with recent years only
  make_plot(
    "spawning_output_recent.png",
    subplots = 2,
    xlim = c(1990, 2025),
    order = order
  )
  make_plot(
    "fraction_unfished_recent.png",
    subplots = 4,
    xlim = c(1990, 2025),
    order = order
  )
  make_plot(
    "summary_biomass_recent.png",
    subplots = 18,
    xlim = c(1990, 2025),
    order = order
  )

  
  # plots with recent years only
  make_plot(
    "spawning_output_recent_no_uncertainty.png",
    subplots = 1,
    xlim = c(1990, 2025),
    order = order
  )
  make_plot(
    "fraction_unfished_recent_no_uncertainty.png",
    subplots = 3,
    xlim = c(1990, 2025),
    order = order
  )
}
