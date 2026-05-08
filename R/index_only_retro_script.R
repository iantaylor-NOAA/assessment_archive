# script to run retrospective analysis removing most recent 10 years of data except for the WCGBTS index

# run script which makes the info2_recent table with info on recent assessments
source("R/load_assessments.R")
# source function to modify model input files for the retrospective
source("R/index_only_retro.R")

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

info3_wcgbts$index_only_retro_dir <- file.path(
  "data-raw/index_only_retro",
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

dir.create("data-raw/index_only_retro", recursive = TRUE, showWarnings = FALSE)

if (FALSE) {
  # modify models to remove 10 years from all data sources other than the WCGBTS
  for (i in 1:nrow(info3_wcgbts)) {
    cli::cli_inform("Writing input files for {info3_wcgbts$Species[i]}...")
    dir <- index_only_retro(
      species = info3_wcgbts$Species[i],
      dir = info3_wcgbts$Model.archive.directory[i],
      newdir = info3_wcgbts$newdir[i],
      WCGBTS_fleet = info3_wcgbts$WCGBTS_number[i]
    )
  }

  for (i in 1:nrow(info3_wcgbts)) {
    cli::cli_inform(
      "Running index-only retrospective for {info3_wcgbts$Species[i]}..."
    )
    r4ss::run(
      dir = info3_wcgbts$index_only_retro_dir[i],
      show_in_console = TRUE
    )
  }
}

# manual changes to Sablefish and Dover models to get them to run with the new data

# problems:
# * Sablefish needs change to Q param section for Recruitment_Index
# * Dover has problem reading ss.par

# fix sablefish
inputs <- r4ss::SS_read("data-raw/index_only_retro/Sablefish")
inputs$start$init_values_src <- 0 # don't use par because the number of parameters has changed
inputs$ctl$Q_options <- inputs$ctl$Q_options |> dplyr::filter_out(fleet == 11)
inputs$ctl$Q_parms <- inputs$ctl$Q_parms |>
  dplyr::filter_out(grepl("Recruitment_Index", rownames(inputs$ctl$Q_parms)))
r4ss::SS_write(
  inputs,
  dir = "data-raw/index_only_retro/Sablefish",
  overwrite = TRUE
)
r4ss::run(dir = "data-raw/index_only_retro/Sablefish", show_in_console = TRUE)

# fix dover
inputs <- r4ss::SS_read("data-raw/index_only_retro/Dover sole")
inputs$start$init_values_src <- 0 # don't use par because the number of parameters has changed
r4ss::SS_write(
  inputs,
  dir = "data-raw/index_only_retro/Dover sole",
  overwrite = TRUE
)
r4ss::run(dir = "data-raw/index_only_retro/Dover sole", show_in_console = TRUE)


# load original models
outputs_original <- r4ss::SSgetoutput(
  dirvec = info3_wcgbts$Model.archive.directory
)

# load all models
output_retro <- r4ss::SSgetoutput(
  dirvec = info3_wcgbts$index_only_retro_dir
)

dir.create("figures/index_only_retro", recursive = TRUE, showWarnings = FALSE)

make_plot <- function(filename, subplots, xlim = NULL) {
  filename = file.path("figures/index_only_retro", filename)
  png(
    filename,
    width = 10,
    height = 6,
    units = "in",
    res = 300
  )
  par(mfrow = c(2, 4), mar = c(2, 3, 2, 1))
  for (i in 1:nrow(info3_wcgbts)) {
    r4ss::SSplotComparisons(
      r4ss::SSsummarize(list(
        output_retro[[i]],
        outputs_original[[i]]
      )),
      legendlabels = c("Index-only retro", "Original"),
      new = FALSE,
      subplots = subplots,
      indexfleets = info3_wcgbts$WCGBTS_number[i],
      indexQlabel = TRUE,
      xlim = xlim
    )
    title(info3_wcgbts$Species[i])
  }
  dev.off()
}

make_plot("index_only_retro_fit_to_WCGBTS.png", subplots = 13)
make_plot("index_only_retro_spawning_output.png", subplots = 2)
make_plot("index_only_retro_fraction_unfished.png", subplots = 4)
make_plot("index_only_retro_summary_biomass.png", subplots = 18)
# recent versions
make_plot(
  "index_only_retro_spawning_output_recent.png",
  subplots = 2,
  xlim = c(1990, 2025)
)
make_plot(
  "index_only_retro_fraction_unfished_recent.png",
  subplots = 4,
  xlim = c(1990, 2025)
)
make_plot(
  "index_only_retro_summary_biomass_recent.png",
  subplots = 18,
  xlim = c(1990, 2025)
)

# new table with quantities of interest
frac_unfished_colname <- paste0("Bratio_", info3_wcgbts$Year + 1)
OFL_colname <- paste0("OFLCatch_", info3_wcgbts$Year + 3)

tab <- info3_wcgbts |>
  dplyr::select(Species, WCGBTS_name) |>
  dplyr::mutate(
    frac_unfished_endyr_plus1 = NA,
    frac_unfished_endyr_plus1_original = NA,
    OFL_endyr_plus_3 = NA,
    OFL_endyr_plus_3_original = NA
  )

for (i in 1:nrow(info3_wcgbts)) {
  tab$frac_unfished_endyr_plus1[i] <- output_retro[[i]]$derived_quants[
    frac_unfished_colname[i],
    "Value"
  ] |>
    round(3)
  tab$frac_unfished_endyr_plus1_original[i] <- outputs_original[[
    i
  ]]$derived_quants[frac_unfished_colname[i], "Value"] |>
    round(3)
  tab$OFL_endyr_plus_3[i] <- output_retro[[i]]$derived_quants[
    OFL_colname[i],
    "Value"
  ] |>
    round(0)
  tab$OFL_endyr_plus_3_original[i] <- outputs_original[[i]]$derived_quants[
    OFL_colname[i],
    "Value"
  ] |>
    round(0)
}

tab <- tab |>
  dplyr::mutate(
    frac_unfished_ratio = round(
      frac_unfished_endyr_plus1 / frac_unfished_endyr_plus1_original,
      2
    ),
    OFL_ratio = round(
      OFL_endyr_plus_3 / OFL_endyr_plus_3_original,
      2
    )
  )
