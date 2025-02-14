# needs working directory to be location of this repository
# can be improved by here::here() or some other good practice
info <- read.csv("inst/extdata/Stock Assessment History - main 14 Feb 2025.csv")
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
  info$newest[info$Species == s &
    info$Year == max(info$Year[info$Species == s])] <- TRUE
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

for(irow in 1:nrow(info2)) {
  if (!file.exists(file.path(info2$Model.archive.directory[irow], "Report.sso"))) {
    cat("missing report in ", info2$Model.archive.directory[irow], "\n")
  }
}

models <- r4ss::SSgetoutput(
  dirvec = info2$Model.archive.directory,
  getcovar = FALSE,
  getcomp = FALSE,
  forecast = FALSE
)

input_files <- list()
for (i in 1:nrow(info2)) {
  if (info2$Year[i] %in% 2021:2024 & 
    info2$Species[i] != "Arrowtooth flounder" &
    !grepl("BlackRF_2023", info2$Model.archive.directory[i])
    ) {
    dir <- info2$Model.archive.directory[i]
    print(i); print(dir)
    input_files[[i]] <- r4ss::SS_read(dir)
  }
}

# # Note sure what I was doing with this incomplete code
# info2$do_recdev <- 
# for(i in c(25:29, 34:40)) {
#    <- input_files
#   print(mean(models[[i]]$recruit$dev, na.rm = TRUE))
# }


### don't run the following by accident, included here only as an example of what I did
if (FALSE) {
  ## save resulting giant list as Rdata file (ignored by git thanks to .gitignore)
  save(info2, models, file = "data/models_14Feb2025.Rdata")
  ## load giant list rather than each individual model
  load("data/models_14Feb2025.Rdata")
}
