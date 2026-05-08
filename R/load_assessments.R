# needs working directory to be location of this repository
# can be improved by here::here() or some other good practice

read_output <- FALSE
read_input <- TRUE

info <- read.csv("inst/extdata/Stock Assessment History - main 2026-05-01.csv")
(species <- sort(unique(info$Species)))
#  [1] "Arrowtooth flounder"              "Aurora rockfish"
#  [3] "Big skate"                        "Black rockfish"
#  [5] "Blackgill rockfish"               "Blue rockfish"
#  [7] "Blue/Deacon rockfish"             "Bocaccio"
#  [9] "Brown rockfish"                   "Cabezon"
# [11] "California scorpionfish"          "Canary rockfish"
# [13] "Chilipepper rockfish"             "China rockfish"
# [15] "Copper rockfish"                  "Cowcod"
# [17] "Darkblotched rockfish"            "Data-Limited"
# [19] "Dover sole"                       "English sole"
# [21] "Gopher rockfish"                  "Gopher/black-and-yellow rockfish"
# [23] "Greenspotted rockfish"            "Greenstriped rockfish"
# [25] "Kelp greenling"                   "Lingcod"
# [27] "Longnose skate"                   "Longspine thornyhead"
# [29] "Pacific hake"                     "Pacific ocean perch"
# [31] "Pacific sanddab"                  "Pacific spiny dogfish"
# [33] "Petrale sole"                     "Quillback rockfish"
# [35] "Rex sole"                         "Rougheye/blackspotted rockfish"
# [37] "Sablefish"                        "Sharpchin rockfish"
# [39] "Shortbelly rockfish"              "Shortspine thornyhead"
# [41] "Splitnose rockfish"               "Squarespot rockfish"
# [43] "Starry flounder"                  "Stripetail rockfish"
# [45] "Vermilion rockfish"               "Vermilion/sunset rockfish"
# [47] "Widow rockfish"                   "Yelloweye rockfish"
# [49] "Yellowtail rockfish"

# figure out the most recent assessment for each species
info$newest <- FALSE
for (s in species) {
  info$newest[
    info$Species == s &
      info$Year == max(info$Year[info$Species == s])
  ] <- TRUE
}

# filter for models from 2015 onward which also have an archive folder listed
info2 <- info |>
  dplyr::filter(
    newest &
      Year >= 2015 &
      nchar(info$Model.archive.directory) > 0
  ) |>
  dplyr::select(
    Species,
    Stock,
    Year,
    Type,
    Model.archive.directory
  )

## # copy files to local directories (potential future option)
##
## info2$local.dir <- paste(info2$Species,
##                          info2$Stock,
##                          info2$Year,
##                          info2$Type)
## info2$local.dir <- gsub(pattern = "[[:blank:]]+",
##                         replacement = "_",
##                         x = info2$local.dir)
##
## # copy models to local file
## for (irow in 1:nrow(info2)) {
##   r4ss::copy_SS_inputs(dir.old=info2$Model.archive.directory[irow],
##                           dir.new=file.path(mydir, "models", info2$local.dir[irow]),
##                           create.dir=TRUE,
##                           overwrite=TRUE,
##                           recursive=FALSE,
##                           use_ss_new=FALSE,
##                           copy_exe=TRUE,
##                           copy_par=TRUE,
##                           verbose=TRUE)
## }

# exclude models with no report files
info2 <- info2 |>
  dplyr::filter(Species != "Cabezon" | Stock != "WA") |> # Cabezon WA used SSS
  dplyr::filter(Species != "Pacific hake") # Hake used MCMC

for (irow in 1:nrow(info2)) {
  if (
    !file.exists(file.path(info2$Model.archive.directory[irow], "Report.sso"))
  ) {
    cat("missing report in ", info2$Model.archive.directory[irow], "\n")
  }
}

# load model output
if (read_output) {
  outputs <- r4ss::SSgetoutput(
    dirvec = info2$Model.archive.directory
  )
}

# filter for only recent full assessments
info2_recent <- info2 |>
  dplyr::filter(Year %in% 2021:2026, Type == "Full")

# load model output for recent full assessments
if (read_output) {
  outputs_recent <- r4ss::SSgetoutput(
    dirvec = info2_recent$Model.archive.directory
  )
}

# load model input for recent full assessments
if (read_input) {
  inputs_recent <- list()
  for (i in 1:nrow(info2_recent)) {
    # if (
    #   info2_recent$Year[i] %in%
    #     2021:2026 &
    #     info2$Species[i] != "Arrowtooth flounder" &
    #     !grepl("BlackRF_2023", info2$Model.archive.directory[i])
    # ) {
    dir <- info2_recent$Model.archive.directory[i]
    print(i)
    print(dir)
    inputs_recent[[i]] <- r4ss::SS_read(dir)
    # }
  }
}

### don't run the following by accident, included here only as an example of what I did
if (FALSE) {
  ## save resulting giant list as Rdata file (ignored by git thanks to .gitignore)
  save(
    info2_recent,
    outputs_recent,
    inputs_recent,
    file = "data/models_recent_full_2026-05-01.Rdata"
  )
  ## load giant list rather than each individual model
  load("data/models_recent_full_2026-05-01.Rdata")
}
