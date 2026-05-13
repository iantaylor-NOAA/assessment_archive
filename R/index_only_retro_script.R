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

# changes to Sablefish and Dover models to get them to run with the new data

# problems:
# * Sablefish needs change to Q param section for Recruitment_Index
# * Dover has problem reading ss.par
# * Canary has some selectivity blocks in starting after the removed years of data, causing bad estimates of selectivity (which impacts the OFLs). Need to modify block design and remove parameters for the blocks that are being removed.

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

# fix dover (the ss.par file from the 2021 SS3 version included extra unused parameters for implementation error)
inputs <- r4ss::SS_read("data-raw/index_only_retro/Dover sole")
inputs$start$init_values_src <- 0 # don't use par 
r4ss::SS_write(
  inputs,
  dir = "data-raw/index_only_retro/Dover sole",
  overwrite = TRUE
)
r4ss::run(dir = "data-raw/index_only_retro/Dover sole", show_in_console = TRUE)

# fix canary
inputs <- r4ss::SS_read("data-raw/index_only_retro/Canary rockfish")
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
  dir = "data-raw/index_only_retro/Canary rockfish_revised_blocks",
  overwrite = TRUE
)
# run model
r4ss::run(
  dir = "data-raw/index_only_retro/Canary rockfish_revised_blocks",
  show_in_console = TRUE
)
# change the directory in the info table
info3_wcgbts$index_only_retro_dir[info3_wcgbts$Species == "Canary rockfish"] <-
  "data-raw/index_only_retro/Canary rockfish_revised_blocks"


# load original models
outputs_original <- r4ss::SSgetoutput(
  dirvec = info3_wcgbts$Model.archive.directory
)

# load retro models
output_retro <- r4ss::SSgetoutput(
  dirvec = info3_wcgbts$index_only_retro_dir
)

dir.create("figures/index_only_retro", recursive = TRUE, showWarnings = FALSE)

make_plot <- function(
  filename,
  subplots,
  xlim = NULL,
  order = 1:nrow(info3_wcgbts)
) {
  filename = file.path("figures/index_only_retro", filename)
  png(
    filename,
    width = 10,
    height = 6,
    units = "in",
    res = 300
  )
  par(mfrow = c(2, 4), mar = c(2, 3, 2, 1))
  for (i in (1:nrow(info3_wcgbts))[order]) {
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

# new table with quantities of interest
frac_unfished_colname <- paste0("Bratio_", info3_wcgbts$Year + 1)
OFL_colname <- paste0("OFLCatch_", info3_wcgbts$Year + 3)

tab <- info3_wcgbts |>
  dplyr::select(Species) |>
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

OFL_ratio_minus_one = abs(tab$OFL_ratio - 1)
order <- order(OFL_ratio_minus_one, decreasing = FALSE)

tab <- tab[order, ]
write.csv(
  tab,
  "data-raw/index_only_retro/index_only_retro_results.csv",
  row.names = FALSE
)


make_plot("index_only_retro_fit_to_WCGBTS.png", subplots = 13, order = order)
make_plot("index_only_retro_spawning_output.png", subplots = 2, order = order)
make_plot("index_only_retro_fraction_unfished.png", subplots = 4, order = order)
make_plot("index_only_retro_summary_biomass.png", subplots = 18, order = order)
# recent versions
make_plot(
  "index_only_retro_spawning_output_recent.png",
  subplots = 2,
  xlim = c(1990, 2025),
  order = order
)
make_plot(
  "index_only_retro_fraction_unfished_recent.png",
  subplots = 4,
  xlim = c(1990, 2025),
  order = order
)
make_plot(
  "index_only_retro_summary_biomass_recent.png",
  subplots = 18,
  xlim = c(1990, 2025),
  order = order
)
