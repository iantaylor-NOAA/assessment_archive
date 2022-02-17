# needs working directory to be location of this repository
# can be improved by here::here() or whatever is some other good practice
info <- read.csv("inst/extdata/Stock Assessment History - master 5 Aug 2021.csv")
(species <- sort(unique(info$Species)))
##  [1] "Arrowtooth flounder"              "Aurora rockfish"                 
##  [3] "Big skate"                        "Black rockfish"                  
##  [5] "Blackgill rockfish"               "Blue rockfish"                   
##  [7] "Blue/Deacon rockfish"             "Bocaccio"                        
##  [9] "Brown rockfish"                   "Cabezon"                         
## [11] "California scorpionfish"          "Canary rockfish"                 
## [13] "Chilipepper rockfish"             "China rockfish"                  
## [15] "Copper rockfish"                  "Cowcod"                          
## [17] "Darkblotched rockfish"            "Data-Limited"                    
## [19] "Dover sole"                       "English sole"                    
## [21] "Gopher rockfish"                  "Gopher/black-and-yellow rockfish"
## [23] "Greenspotted rockfish"            "Greenstriped rockfish"           
## [25] "Kelp greenling"                   "Lingcod"                         
## [27] "Longnose skate"                   "Longspine thornyhead"            
## [29] "Pacific hake"                     "Pacific ocean perch"             
## [31] "Pacific sanddab"                  "Pacific spiny dogfish"           
## [33] "Petrale sole"                     "Quillback rockfish"              
## [35] "Rex sole"                         "Rougheye/blackspotted rockfish"  
## [37] "Sablefish"                        "Sharpchin rockfish"              
## [39] "Shortbelly rockfish"              "Shortspine thornyhead"           
## [41] "Splitnose rockfish"               "Squarespot rockfish"             
## [43] "Starry flounder"                  "Stripetail rockfish"             
## [45] "Vermilion rockfish"               "Vermilion/sunset rockfish"       
## [47] "Widow rockfish"                   "Yelloweye rockfish"              
## [49] "Yellowtail rockfish"             

# figure out the most recent assessment for each species
info$newest <- FALSE
for(s in species){
  info$newest[info$Species == s &
              info$Year == max(info$Year[info$Species == s])] <- TRUE
}

# filter for models from 2015 onward which also have an archive folder listed
info2 <- info[info$newest &
              info$Year >= 2015 &
              nchar(info$Model.archive.directory) > 0,]
info2 <- info2[,c("Species",
                  "Stock",
                  "Year",
                  "Type",
                  "Model.archive.directory")]

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

# skipping model with no report file (used SSS and has different format for output)
info2 <- info2[!(info2$Species == "Cabezon" & info2$Stock == "WA"),]

models <- r4ss::SSgetoutput(dirvec = info2$Model.archive.directory,
                            getcovar = FALSE,
                            getcomp = FALSE,
                            forecast = FALSE)

### don't run the following by accident, included here only as an example of what I did
if (FALSE) {
  ## save resulting giant list as Rdata file
  save(info2, models, file = 'assessment_histories/models_17Feb2022.Rdata')
  ## load giant list rather than each individual model
  load('assessment_histories/models_17Feb2022.Rdata')
}
