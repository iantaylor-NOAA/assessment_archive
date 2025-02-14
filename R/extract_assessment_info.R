# 20 June 2024 getting info on extraSD parameters
# incomplete work to get names of all fleets with indices
# figure out which ones have Q_extraSD parameters
# create table of surveys with and without Q_extraSD parameters

models[[1]]$index_variance_tuning_check |> 
  dplyr::filter(`Input+VarAdj` > 0) |>
  dplyr::mutate(extraSD = as.numeric(`Input+VarAdj+extra`) - 
    as.numeric(`Input+VarAdj`)) |> 
  dplyr::select(extraSD, Name)

# NEXT: rename fleetname to Name
for(i in 1:40) {
  print(i)
  models[[i]]$index_variance_tuning_check |> 
    dplyr::filter(`Input+VarAdj` > 0) |>
    dplyr::mutate(extraSD = as.numeric(`Input+VarAdj+extra`) - 
      as.numeric(`Input+VarAdj`)) |> 
    dplyr::rename(Name = tidyselect::any_of(fleetname))
    dplyr::select(extraSD, Name)
}
#### old stuff not updated on August 5, 2021
#### depends on the info2 list created by the load_assessments.R script

info_small <- info2 [,c("Species",
                       "Stock",
                       "Type"),]
phase.info <- NULL
for(irow in 1:nrow(info.small)){
  if (info.small$Type[irow] == "Full") {
    message(paste(irow, paste(info.small[irow,], collapse = " ")))
    pars <- models[[irow]]$parameters
    pars <- pars[!is.na(pars$Label) & pars$Phase > 0,]
    phase.info <- rbind(phase.info,
                        data.frame(info.small[irow,], pars[,c("Label", "Phase")]))
  }
}

# get table with phase info
info.small <- info2[,c("Species",
                       "Stock",
                       "Type"),]
phase.info <- NULL
for(irow in 1:nrow(info.small)){
  if (info.small$Type[irow] == "Full") {
    message(paste(irow, paste(info.small[irow,], collapse = " ")))
    pars <- models[[irow]]$parameters
    pars <- pars[!is.na(pars$Label) & pars$Phase > 0,]
    phase.info <- rbind(phase.info,
                        data.frame(info.small[irow,], pars[,c("Label", "Phase")]))
  }
}
                                 


  # get names of data files from starter files
  #source('c:/SS/R/SS_readstarter2.R')
  n <- length(newestfull)
  dat <- rep(NA,n)
  ctl <- rep(NA,n)
  for(i in 1:n){
    #start1 <- SS_readstarter2(paste(newestfull[i]))
    start1 <- SS_readstarter(paste(newestfull[i]))
    dat[i] <- start1$datfile
    ctl[i] <- start1$ctlfile
  }

  # get directories
  newest1 <- strsplit(newestfull,'starter')
  newest2 <- rep(NA,n)
  for(i in 1:n) newest2[i] <- newest1[[i]][1]
  newest3 <- strsplit(newest2,'Starter')
  newest2 <- rep(NA,n)
  for(i in 1:n) newest2[i] <- newest3[[i]][1]

  # this file had ss2.rep while the automatically selected one didn't
  black <- grep("Black_RF",newest2)
  newest2[black] <- "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Black_RF\\Blackrf_2007\\ModelFiles_North_and_South\\BlackRF_N\\"
  dat[black] <- "BRFBase.dat"
  # manually changing to directory with base model files
  newest2[grep("Dover",newest2)] <- "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\DoverSole\\DoverSole_2011\\DoverSole2011_baseModelFiles\\"
  # manually changing hake away from retrospective directory to base model
  newest2[grep("Hake",newest2)] <- "\\\\nwcfile\\fram\\Assessments\\Archive temp\\Pacific Hake_Whiting\\PWhiting2014\\2014hake_21_TVselex1991start_MLE\\"
  
  datfile <- paste0(newest2,dat)
  ##  [1] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\ArrowtoothFl/Arrowtooth_2007/ArrowtoothFl_2007_ModelFiles/arrowtooth.dat"                                 
  ##  [2] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Aurora/AuroraRF_2013/Aurora2013_BaseModelInputFiles/ARRA_dat3.ss"                                         
  ##  [3] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Black_RF\\Blackrf_2007\\ModelFiles_North_and_South\\BlackRF_N\\BRFBase.dat"                               
  ##  [4] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\BlackgillRF/BlackgillRockfish_2011/ModelFiles_BlackgillRF_2011/bgill.star36.dat"                          
  ##  [5] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Blue_Rockfish/ModelFiles_BlueRf/finalblue.DAT"                                                            
  ##  [6] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Bocaccio/Bocaccio_2013/ModelFiles_Bocaccio2013/boc8.dat"                                                  
  ##  [7] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Cabezon/Cabezon_2009/Cabezon_ModelFiles_2009/SCS/BC_tuned.dat"                                            
  ##  [8] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\CanaryRF/CanaryRf_2011/Canary base case model files/Canary_data.SS"                                       
  ##  [9] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Chilipepper/ChilipepperRf_2007/ChilipepperRF_ModelFiles_2007_CorrectedBase/chilibase.dat"                 
  ## [10] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Cowcod/Cowcod_2009/ModelFiles_Cowcod_2009/cowcod_2009/base_model/moo4_base.dat"                           
  ## [11] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Darkblotched/Darkblotched_2013/Model_files/darkblotched_data.SS"                                          
  ## [12] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\DoverSole\\DoverSole_2011\\DoverSole2011_baseModelFiles\\dover2011.dat"                                   
  ## [13] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\EnglishSole/EnglishSole_2007/EnglishSole_2007ModelFiles/english.dat"                                      
  ## [14] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\GreenspottedRF/GreenspottedRF_2011/GreenspottedRF_ModelFiles_2011/GSPT_South_files_for_Archive/GSPT_S.dat"
  ## [15] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\GreenstripedRF/GreenstripedRF_modelfiles_2009/gsrk09.dat"                                                 
  ## [16] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Lingcod/Lingcod_2009/Lingcod_ModelFiles_2009/North/LingN_data.SS"                                         
  ## [17] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\LongnoseSkate/LongnoseSkate_2007/LongnoseSkate_ModelFiles/Vlada.dat"                                      
  ## [18] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\LongspineTH/Longspine Thornyhead_2013/Base/LST_data.SS"                                                   
  ## [19] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Pacific Hake_Whiting/PWhiting2014/2014hake_21_retros/2014hake_21_retro02/2014hake_data.SS"                
  ## [20] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\PacificOceanPerch/POP_2011/ModelFiles_POP_2011/POP_data.SS"                                               
  ## [21] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\PetraleSole/PetraleSole_2013/PetraleBaseModelFiles_2013/petrale13.dat"                                    
  ## [22] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\RougheyeRF_BlackspottedRF/base/REYE.dat"                                                                  
  ## [23] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Sablefish/Sablefish2011/2011 sablefish model files/2011_sablefish_data.SS"                                
  ## [24] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\ShortspineThornyhead/SST_2013/base_model_files/SST_data.SS"                                               
  ## [25] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\SpinyDogfish/SpinyDogfish_2011/ModelFiles_Spinydogfish_2011/Spiny_Dogfish.DAT"                            
  ## [26] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\SplitnoseRf/Splitnose_2009/SS_files/Splitnose_Rf.dat"                                                     
  ## [27] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\Widow/WidowRF_2011/Widow2011base/wdw1.dat"                                                                
# note: fix yelloweye
  ## [28] "\\\\nwcfile\\FRAM\\Assessments\\Archive temp\\YelloweyeRF/YelloweyeRF_2011/rebuilding/states/cat100_hhigh/yelloweye_data.SS"                            

  ctlfile <- paste0(newest2,ctl)

  # stuff for binning group
  compression <- rep(NA, n)
  constant <- rep(NA, n)
  for(i in 1:n){
    lines <- readLines(datfile[i])
    print(dirsallunique[i])
    tmp <- print(lines[grep("compress", lines)])
    #compression[i] <- as.numeric(strsplit(tmp, "
    #print(lines[grep("compress", lines)+1])
    
  }    
    
    

  
  n <- length(newest2)

  fishnames <- c("Arrowtooth",
                 "Aurora",
                 "Black",
                 "Blackgill",
                 "Blue",
                 "Bocaccio",
                 "Cabezon",
                 "Canary",
                 "Chilipepper",
                 "Cowcod",
                 "Darkblotched",
                 "Dover",
                 "English",
                 "Greenspotted",
                 "Greenstriped",
                 "Lingcod",
                 "Longnose",
                 "Longspine",
                 "Hake",
                 "POP",
                 "Petrale",
                 "Rougheye",
                 "Sablefish",
                 "Shortspine",
                 "Dogfish",
                 "Splitnose",
                 "Widow",
                 "Yelloweye")


  # read data files
  myreadLines <- list()
  good <- rep(NA,n)
  for(i in 1:n) good[i] <- !is.na(file.info(datfile[i])$size)
  #for(i in 1:n) if(!good[i]) print(dir(newest2[i]))
  for(i in 1:n) myreadLines[[i]] <- readLines(datfile[i],warn=F)

  # make table
  info <- data.frame(name=fishnames, NfishFleets=NA, NsurveyFleets=NA)

  # summarize fleets
  for(i in 1:n){
    # figure out number of fishing fleets
    fishFleetsLine <- grep("fishing",tolower(myreadLines[[i]][3:15]))[1]
    if(length(fishFleetsLine)==0 || is.na(fishFleetsLine)){
      fishFleetsLine <- grep("fleet",tolower(myreadLines[[i]][3:15]))[1]
    }
    if(length(fishFleetsLine)==0 || is.na(fishFleetsLine)){
      fishFleetsLine <- grep("fisheries",tolower(myreadLines[[i]][3:15]))[1]
    }
    if(length(fishFleetsLine)>0 && !is.na(fishFleetsLine)){
      info$NfishFleets[i] <- as.numeric(strsplit(tolower(myreadLines[[i]][3:15])[fishFleetsLine],"[[:blank:]]+")[[1]][1])
    }
    # figure out number of surveys
    surveyFleetsLine <- grep("surv",tolower(myreadLines[[i]][3:15]))[1]
    ## if(length(surveyFleetsLine)==0 || is.na(surveyFleetsLine)){
    ##   surveyFleetsLine <- grep("surv",tolower(myreadLines[[i]])[3:15])[1]
    ## }
    if(length(surveyFleetsLine)>0 && !is.na(surveyFleetsLine)){
      info$NsurveyFleets[i] <- as.numeric(strsplit(tolower(myreadLines[[i]][3:15])[surveyFleetsLine],"[[:blank:]]+")[[1]][1])
    }
  }
  # manual entry for Blue RF
  info$NfishFleets[info$name=="Blue"] <- 3
  info$NsurveyFleets[info$name=="Blue"] <- 5
  # total fleets
  info$Nfleets <- info$NfishFleets + info$NsurveyFleets

  png('C:/Users/tayloria/Documents/talks/ProgrammaticReview2014/Nfleets.png',
      res=300,units='in',width=6,height=6)
  par(mar=c(4.1,4.1,.3,.3))
  plot(0,xlim=c(0,12.9),ylim=c(0,12.9),type='n',xaxs='i',yaxs='i',las=1,
       xlab="Number of modeled fishing fleets",ylab="Number of modeled indices")
  #abline(h=1:12,v=1:12,col=grey(.8),lty=3)
  grid()
  points(info$NfishFleets,info$NsurveyFleets,pch=16,cex=.8)
  text(jitter(info$NfishFleets),jitter(info$NsurveyFleets),labels=info$name,cex=.8,pos=3)
  dev.off()
  
  # read Report files
  repfile <- paste(newest2,"Report.sso",sep='')
  myreadRep <- list()
  good <- rep(NA,n)
  for(i in 1:n) good[i] <- !is.na(file.info(repfile[i])$size)
  for(i in 1:n) if(!good[i]){
    repfile[i] <- paste(newest2[i],"ss2.rep",sep='')
  }
  repfile[grep("Aurora",repfile)] <-
    "\\\\nwcfile\\fram\\Assessments\\Archive temp\\Aurora\\AuroraRF_2013\\Aurora2013_BaseModelOutput\\Report.sso"
  repfile[grep("Darkblotched",repfile)] <-
    "\\\\nwcfile\\fram\\Assessments\\Archive temp\\Darkblotched\\Darkblotched_2013\\Model_output\\Report.sso"
  repfile[grep("Hake",repfile)] <-
    "\\\\nwcfile\\fram\\Assessments\\Archive temp\\Pacific Hake_Whiting\\PWhiting2014\\2014hake_21_TVselex1991start_MLE\\Report.sso"
  
  good <- rep(NA,n)
  for(i in 1:n) good[i] <- !is.na(file.info(repfile[i])$size)
  for(i in 1:n) if(good[i]) myreadRep[[i]] <- readLines(repfile[i],warn=F)

  info$time <- NA
  info$year <- NA
  for(i in 1:n){
    time1 <- myreadRep[[i]][grep("Time: ",myreadRep[[i]][1:10])[1]]
    time2 <- strsplit(time1,": ")[[1]][2]
    info$time[i] <- time2
    info$year[i] <- as.numeric(substring(time2,21))
  }


  outputdirs <- newest2
  outputdirs[grep("Aurora",outputdirs)] <-
    "\\\\nwcfile\\fram\\Assessments\\Archive temp\\Aurora\\AuroraRF_2013\\Aurora2013_BaseModelOutput\\"
  outputdirs[grep("Darkblotched",outputdirs)] <-
    "\\\\nwcfile\\fram\\Assessments\\Archive temp\\Darkblotched\\Darkblotched_2013\\Model_output\\"

  # get output
  from2011plus <- which(info$year >= 2011 & info$name!="Chilipepper")
  info$from2011plus <- (1:nrow(info)) %in% from2011plus
  
  #from2011 <- c(3,5,7,10,11,16,17,18,19,20,22,23)
  fishnames[from2011plus]
  ##  [1] "Aurora"       "Blackgill"    "Bocaccio"     "Canary"       "Darkblotched" "Dover"       
  ##  [7] "Greenspotted" "Longspine"    "Hake"         "POP"          "Petrale"      "Rougheye"    
  ## [13] "Sablefish"    "Shortspine"   "Dogfish"      "Widow"        "Yelloweye"
  alloutput <- list()
  for(i in from2011plus){
    alloutput[[i]] <- SS_output(dir=outputdirs[i],NoCompOK=TRUE,covar=FALSE,forecast=FALSE)
  }

  ####
  ### stuff added 11/4/2014 to export output
  # 
  names(alloutput) <- fishnames
  output2011to2014 <- alloutput[from2011plus]
  #save(output2011to2014, file="c:/SS/R/assessment_histories2011to2014.Rdata")
  
  load("c:/SS/R/assessment_histories/assessment_histories2011to2014.Rdata")
  
  getselage <- function(replist){
    Z_at_age <- replist$Z_at_age[replist$Z_at_age$Gender==1 & replist$Z_at_age$Bio_Pattern==1,]
    M_at_age <- replist$M_at_age[replist$M_at_age$Gender==1 & replist$M_at_age$Bio_Pattern==1,]
    yrs <- Z_at_age$Year
    Z_at_age <- as.matrix(Z_at_age[yrs<2012,-(1:3)])
    M_at_age <- as.matrix(M_at_age[yrs<2012,-(1:3)])
    yrs <- yrs[yrs<2012]
    nyrs <- length(yrs)
    F_at_age <- Z_at_age - M_at_age
    Sel_at_age <- 0*Z_at_age
    for(iyr in 1:nyrs){
      Sel_at_age[iyr,] <- F_at_age[iyr,]/max(F_at_age[iyr,],na.rm=TRUE)
    }
    maxage <- ncol(F_at_age)-1
    #image(0:(maxage-1), yrs, t(Sel_at_age[,1:maxage]),
    #      xlab="Age",ylab="Year",col=grey(seq(1,0,-.1)))
    Sel_ave <- apply(Sel_at_age,2,mean,na.rm=TRUE)
    Sel_ave <- Sel_ave/max(Sel_ave,na.rm=TRUE)
    plot(0:maxage,Sel_ave,type='l',lwd=3,ylim=c(0,1.1),yaxs='i',xaxs='i',axes=FALSE)
    axis(1)
    matplot(0:maxage,t(Sel_at_age),col=rgb(0,0,0,.1),lty=1,type='l',add=TRUE)
  }

    png("Scaled_Fs_for_recent_west-coast_assessments.png",width=10,height=14,res=300,units='in')
    par(mfrow=c(5,4),mar=rep(2,4),oma=c(3,3,3,3))
    for(i in 1:length(output2011to2014)){
      getselage(output2011to2014[[i]])
      #text(x=0,y=1.05,fishnames[i],pos=4,cex=2)
      text(x=0,y=1.05,names(output2011to2014)[i],pos=4,cex=2)
      box()
      axis(2)
    }
    mtext(side=3,"Scaled F-at-age summed across fleets by year (with average in bold)",line=0,outer=TRUE,cex=1.5)
    mtext(side=1,"Age",line=1.5,outer=TRUE,cex=1)
    mtext(side=2,"Scaled F-at-age",line=1.5,outer=TRUE,cex=1)
    dev.off()

    png("Selectivity-at-length_for_recent_west-coast_assessments.png",width=10,height=14,res=300,units='in')
    par(mfrow=c(5,4),mar=rep(2,4),oma=c(3,3,3,3))
    for(i in 1:length(output2011to2014)){
      lbinspop <- as.numeric(names(output2011to2014[[i]]$sizesel[,-(1:5)]))
      plot(0,type='n',ylim=c(0,1.1),xlim=range(lbinspop))
      SSplotSelex(output2011to2014[[i]],subplot=1,add=TRUE,legend=FALSE)
      #text(x=min(lbinspop),y=1.05,fishnames[i],pos=4,cex=2)
      text(x=min(lbinspop),y=1.05,names(output2011to2014)[i],pos=4,cex=2)
      box()
    }
    mtext(side=3,"Length-based selectivity by fleet",line=0,outer=TRUE,cex=1.5)
    mtext(side=1,"Length",line=1.5,outer=TRUE,cex=1)
    mtext(side=2,"Selectivity",line=1.5,outer=TRUE,cex=1)
    dev.off()
    
    
  # get recruitment section
  recruits <- list()
  for(i in 1:n){
    if(good[i]){
      mylines <- myreadRep[[i]]
      reclines <- max(grep("SPAWN_RECRUIT",mylines))+0:200
      headline <- grep("pred",mylines)
      headline <- headline[headline %in% reclines]
      virgline <- grep("Virg",mylines)
      virgline <- min(virgline[virgline %in% reclines])
      endline <- grep("$^",mylines)
      endline <- min(endline[endline %in% reclines])-1
      print(c(headline,virgline,endline))
      ## sr <- mylines[max(grep("SPAWN_RECRUIT",mylines))+0:200]
      ## sr <- sr[min(grep("Virg",sr)):(min(grep("$^",sr))-1)]
      ## sr <- as.data.frame(
      ## print(head(sr))
      recruits[[i]] <- read.table(repfile[i],skip=virgline-1,nrow=(endline-virgline),
                                  colClasses="character",fill=TRUE,row.names=NULL,col.names=1:50)
      ## for(i in 1:ncol(recruits[[i]]))
      ## names(recruits[[i]]) <- read.table(repfile[i],skip=headline-1,nrow=1,
      ##                                   stringsAsFactors=FALSE,fill=TRUE,row.names=NULL)
    }
  }

  pdf('c:/SS/R/Jake_Rice_style_recruit_dev_analysis_2013.pdf',width=8.5,height=11)
  devtable <- data.frame(species=fishnames,positive=NA,negative=NA,fraction_positive=NA)
  par(mfrow=c(6,4),mar=rep(0.5,4),oma=c(5,5,5,1))
  for(i in 1:n){
    if(!is.null(recruits[[i]])){
      recruits[[i]][recruits[[i]]==""] <- NA
      recruits[[i]][recruits[[i]]=="_"] <- NA
      recruits[[i]][recruits[[i]]=="-"] <- NA
      rec <- as.numeric(recruits[[i]][,7])
      spb <- as.numeric(recruits[[i]][,2])
      low <- spb/spb[1]<.4
      plot(spb/spb[1],rec,
           xlim=c(0,1.3),
           xaxs='i',
           #ylim=c(-1,1.1)*(.1+max(abs(rec),na.rm=TRUE)),
           ylim=c(-2.3,2.3),
           type='p',pch=16,cex=1,col=rgb(ifelse(low,1,0),0,ifelse(low,0,1),.5),
           axes=FALSE)
      #print(table(rec[rec!=0 & spb/spb[1]<.4]>0))
      mtext(side=3,line=-1.5,fishnames[i])
      hi40 <- sum(rec>0 & spb/spb[1]<.4)
      lo40 <- sum(rec<0 & spb/spb[1]<.4)
      if(length(hi40)==0 | is.na(hi40)) hi40 <- 0
      if(length(lo40)==0 | is.na(lo40)) lo40 <- 0
      devtable$positive[i] <- hi40
      devtable$negative[i] <- lo40
      if(hi40+lo40>0){
        devtable$fraction_positive[i] <- round(100*hi40/(lo40+hi40))/100
        text(x=.65,y=2,paste("+",hi40,", -",lo40,", ",round(100*hi40/(lo40+hi40)),"%>0",sep=""),pos=1)
      }else{
        text(x=.65,y=2,"No non-zero samples",pos=1)
      }
      abline(h=0)
      #grid()
      abline(v=c(.4,1),col='grey',lty=2)
      box(lwd=1,col='grey')
      if(par()$mfg[1]==par()$mfg[3]) axis(1,col='grey')
      if(par()$mfg[1]==par()$mfg[3]-1 & par()$mfg[2]>3) axis(1,col='grey')
      if(par()$mfg[2]==1) axis(2,las=1,col='grey')
    }
  }
  mtext("Spawning depletion",side=1,line=3,outer=TRUE)
  mtext("Recruitment deviation",side=2,line=3,outer=TRUE)
  devtable <- rbind(devtable,data.frame(species="TOTAL",
                                        positive=sum(devtable$positive,na.rm=TRUE),
                                        negative=sum(devtable$negative,na.rm=TRUE),
                                        fraction_positive=sum(devtable$positive,na.rm=TRUE)/(sum(devtable$positive,na.rm=TRUE)+sum(devtable$negative,na.rm=TRUE))))
  plot(0,type='n',axes=F,
       xlim=c(0,1.3),
       xaxs='i',
       ylim=c(-2.4,2.4))
  mtext(side=3,line=-5.5,"TOTAL",font=2)
  hi40 <- devtable$positive[devtable$species=="TOTAL"]
  lo40 <- devtable$negative[devtable$species=="TOTAL"]
  text(x=.65,y=0,paste("+",hi40,", -",lo40,", ",round(100*hi40/(lo40+hi40)),"%>0",sep=""),pos=1)
  title(main="Analysis of recruitment deviations at spawning depletion < 40%", outer=TRUE)
  dev.off()
  # end Jake Rice style plot

  # list recent models
  recent <- sort(unique(c((1:n)[grep("2011",datfile)],
                          (1:n)[grep("2012",datfile)],
                          (1:n)[grep("2013",datfile)],
                          (1:n)[grep("2014",datfile)],
                          (1:n)[grep("Rougheye",datfile)]
                          )))
  
  N <- length(recent)

  # read data files
  datlists <- list()
  replists <- list()
  recentnames <- NULL

  for(imodel in 1:N){
    i <- recent[imodel]
    recentnames[imodel] <- fishnames[i]
    cat(i," ",recentnames,"\n")
    datlists[[imodel]] <- SS_readdat(datfile[i])
    replists[[imodel]] <- SS_output(outputdirs[i],covar=FALSE,forecast=FALSE,NoCompOK=TRUE,printstats=FALSE)
  }

} # end if(FALSE)

getsolo <- function(yrs,yval=1){
  # function to identify solo years separately from adjacent years
  x <- min(yrs):max(yrs)
  n <- length(x)
  y <- rep(yval,n)
  y[!x%in%yrs] <- NA
  # identify solo points (no data from adjacent years)
  solo <- rep(FALSE,n)
  if(n==1) solo <- 1
  if(n==2 & yrs[2]!=yrs[1]+1) solo <- rep(TRUE,2)
  if(n>=3){
    for(i in 2:(n-1)) if(is.na(y[i-1]) & is.na(y[i+1])) solo[i] <- TRUE
    if(is.na(y[2])) solo[1] <- TRUE
    if(is.na(y[n-1])) solo[n] <- TRUE
  }
  return(data.frame(x,y,solo))
}

#pdf('c:/SS/R/recruit_analysis_2014.pdf',width=8.5,height=11)
order <- 1:N

npars.all <- rep(NA,N)
for(i in 1:N){
  endyr <- replists[[i]]$endyr
  goodpars <- replists[[i]]$parameters[!is.na(replists[[i]]$parameters$Active_Cnt) & 
                                       replists[[i]]$parameters$Active_Cnt > 0 &
                                       !replists[[i]]$parameters$Label %in%
                                       paste("ForeRecr_",endyr+1:30,sep="") &
                                       !replists[[i]]$parameters$Label %in%
                                       paste("Impl_err_",endyr+1:30,sep=""),]
  npars.all[i] <- nrow(goodpars)
}
order <- rev(order(npars.all))

version.info <- data.frame(name=recentnames,version=NA)
for(i in 1:N){
  version.info$version[i] <- replists[[i]]$SS_version
}
version.info$short <- substring(version.info$version,1,9)
version.info$date <- NA
for(i in 1:N){
  tmp <- strsplit(version.info$version[i],";_",fixed=TRUE)[[1]][2]
  date <- as.Date(tmp,format="%m/%d/%Y")
  version.info$date[i] <- date
  version.info$date.char[i] <- as.character(date)
}


### make plot

# options related to 2 pages for powerpoint or just 1 page
option <- 1
if(option==1){
  order2 <- order
  png('c:/SS/R/recruit_analysis_2014b.png',width=8.5,height=11,units='in',res=400)
  par(mfrow=c(N,1),mar=rep(0,4),oma=c(3,2,7,1))
}
if(option==2){
  order2 <- order[1:8]
  png('c:/SS/R/recruit_analysis_2014c.png',width=8.5,height=6,units='in',res=400)
  par(mfrow=c(N/2,1),mar=rep(0,4),oma=c(3,2,7,1))
}
if(option==3){
  order2 <- order[9:16]
  png('c:/SS/R/recruit_analysis_2014d.png',width=8.5,height=6,units='in',res=400)
  par(mfrow=c(N/2,1),mar=rep(0,4),oma=c(3,2,7,1))
}

colvec <- rich.colors.short(N,alpha=.7)
even <- function(x){x%%2==0}

# empty table to store parameter counts
par.table <- data.frame(name=rep(NA,length(order2)))
  
for(j in 1:length(order2)){
  i <- order2[j]
  yr <- replists[[i]]$sprseries$Year
  cat <- replists[[i]]$sprseries$Dead_Catch
  cat <- cat/max(cat)
  plot(0,xlim=c(1850,2033),ylim=c(-3,1.3),xaxs='i',yaxs='i',type='n',xlab='',ylab='',axes=FALSE)
  #abline(v=seq(1850,2010,10),col='grey',lty='11')
  #for(y in seq(1860,2010,20)) rect(y,-5,y+10,5,col='grey98',border=NA)
  abline(v=seq(1860,2010,10), lty=3, col='grey90')
  if(even(j)) rect(1840,-5,2030,5,col=grey(.5,alpha=.1),border=NA)
  if(all(par()$mfg[1]==1)) axis(3,at=c(seq(1850,2014,10),2014))

  # add catch series
  yy <- 2014
  polygon(c(yr[1],yr,rev(yr)[1]),c(0,cat,0),
          col=colvec[j+ifelse(option==3,8,0)],border='grey20')
  cex.val <- 0.8
  text("Catch",     x=yy+0.2,y= 0.6, adj=c(0,.5),cex=cex.val)
  text("Lengths",   x=yy+0.2,y=-0.25,adj=c(0,.5),cex=cex.val)
  text("Ages",      x=yy+0.2,y=-0.8, adj=c(0,.5),cex=cex.val)
  text("Rec. devs.",x=yy+0.2,y=-1.6, adj=c(0,.5),cex=cex.val)
  # add points and lines for length and age data
  lenyrs <- datlists[[i]]$lencomp$Yr
  if(!is.null(lenyrs)) lenyrs <- sort(unique(lenyrs))
  ageyrs <- datlists[[i]]$agecomp$Yr
  if(!is.null(ageyrs)) ageyrs <- sort(unique(ageyrs))
  
  if(!is.null(lenyrs)){
    lenyrs <- getsolo(yrs=lenyrs,yval=-0.3)
    #print(recentnames[i])
    #print(lenyrs)
    points(lenyrs$x[lenyrs$solo], lenyrs$y[lenyrs$solo], pch=16, cex=1,col='grey70')
    lines(lenyrs$x, lenyrs$y, lwd=5, col='grey70')
  }
  if(!is.null(ageyrs)){
    ageyrs <- getsolo(yrs=sort(unique(datlists[[i]]$agecomp$Yr)),yval=-0.7)
    points(ageyrs$x[ageyrs$solo], ageyrs$y[ageyrs$solo], pch=16, cex=1,col='grey40')
    lines(ageyrs$x, ageyrs$y, lwd=5, col='grey40')
  }
  endyr <- datlists[[i]]$endyr
  recruit <- replists[[i]]$recruit
  recruit <- recruit[!is.na(recruit$dev) & recruit$year<=endyr,]
  devyr   <- recruit$year
  dev     <- recruit$dev
  devera  <- recruit$era
  #if(any(dev!=0)) dev <- dev/max(abs(dev))
  dev <- dev/2

  yval <- -1.5
  lines(devyr,rep(yval,length(devyr)),lwd=0.5,col=1)
  lines(devyr,yval+0.9*dev,lwd=2,col="red3")
  lines(x=devyr[devera=="Main"],y=yval+0.9*dev[devera=="Main"],lwd=2,col=2)
#  if(par()$mfg[1]!=par()$mfg[3]) abline(h=-2.8,col='grey50',lwd=2)
  #text(1850,.4,paste(recentnames[i]," (",endyr+1,")",sep=""),pos=4,cex=2)
  # add name of species and other info
  text1 <- recentnames[i]
  endyr <- replists[[i]]$endyr
  text1 <- paste(text1," (",endyr+1,")",sep="")
  text(1850, 0.4, text1, pos=4, cex=2)
  par.table$name[j] <- recentnames[i]
  par.table$year[j] <- endyr+1
  
  
  #info on estimated parameters
  maxpars <- max(replists[[i]]$parameters$Active_Cnt,na.rm=TRUE)
  goodpars <- replists[[i]]$parameters[!is.na(replists[[i]]$parameters$Active_Cnt) & 
                                       replists[[i]]$parameters$Active_Cnt > 0 &
                                       !replists[[i]]$parameters$Label %in%
                                       paste("ForeRecr_",endyr+1:30,sep="") &
                                       !replists[[i]]$parameters$Label %in%
                                       paste("Impl_err_",endyr+1:30,sep=""),]
  goodpars$selex <- FALSE
  goodpars$selex[grep("Sel",goodpars$Label)] <- TRUE
  goodpars$selex[grep("Spline",goodpars$Label)] <- TRUE
  goodpars$retain <- FALSE
  goodpars$retain[grep("Retain",goodpars$Label)] <- TRUE
  goodpars$recdev <- FALSE
  goodpars$recdev[grep("RecrDev",goodpars$Label)] <- TRUE
  goodpars$recdev[grep("InitAge",goodpars$Label)] <- TRUE
  goodpars$recdev[grep("ForeRecr",goodpars$Label)] <- TRUE
  goodpars$growth <- FALSE
  goodpars$growth[grep("L_at_Am",goodpars$Label)] <- TRUE
  goodpars$growth[grep("VonBert",goodpars$Label)] <- TRUE
  goodpars$growth[grep("CV_young",goodpars$Label)] <- TRUE
  goodpars$growth[grep("CV_old",goodpars$Label)] <- TRUE
  goodpars$SR <- FALSE
  goodpars$SR[grep("SR_",goodpars$Label)] <- TRUE
  goodpars$SR[grep("RecrDist_",goodpars$Label)] <- TRUE
  goodpars$M <- FALSE
  goodpars$M[grep("NatM_",goodpars$Label)] <- TRUE
  goodpars$Q <- FALSE
  goodpars$Q[grep("Q_",goodpars$Label)] <- TRUE
  npars <- nrow(goodpars)
  npars.selex <- sum(goodpars$selex)
  npars.retain <- sum(goodpars$retain)
  npars.recdev <- sum(goodpars$recdev)
  npars.growth <- sum(goodpars$growth)
  npars.SR <- sum(goodpars$SR)
  npars.M <- sum(goodpars$M)
  npars.Q <- sum(goodpars$Q)
  # fixing yelloweye since model chosen is state of nature with fixed h
  if(i==16){
    npars <- npars+1
    npars.SR <- npars.SR+1
  }
print(goodpars$Label[!goodpars$selex & !goodpars$retain & !goodpars$recdev &
                     !goodpars$growth & !goodpars$SR & !goodpars$M & !goodpars$Q])
  text2 <- paste("Estimated parameters: ",npars)
  text3 <- paste("S-R: ",npars.SR,
                 ", M: ",npars.M,
                 ", Growth: ", npars.growth,
                 ", Sel/Ret/Q: ", npars.selex+npars.retain+npars.Q,
                 ", Rec. devs: ",npars.recdev, sep="")

  par.table$npars[j] <- npars
  par.table$npars.SR[j] <- npars.SR
  par.table$npars.M[j] <- npars.M
  par.table$npars.growth[j] <- npars.growth
  par.table$npars.selex[j] <- npars.selex
  par.table$npars.retain[j] <- npars.retain
  par.table$npars.Q[j] <- npars.Q
  par.table$npars.recdev[j] <- npars.recdev
  
  
  text(1850, -1, text2, pos=4, cex=1,font=2)
  text(1850, -2, text3, pos=4, cex=1)
}
if(option==1) mtext(side=3,"West Coast Benchmark Assessments 2011-2014",outer=TRUE,cex=2,line=4)
if(option==2) mtext(side=3,"West Coast Benchmark Assessments 2011-2014 (pg. 1/2)",outer=TRUE,cex=1.7,line=4)
if(option==3) mtext(side=3,"West Coast Benchmark Assessments 2011-2014 (pg. 2/2)",outer=TRUE,cex=1.7,line=4)
axis(1,at=c(seq(1850,2014,10),2014))
dev.off()

for(i in 1:nrow(par.table)){
  par.table$NfishFleets[i] <- info$NfishFleets[info$name %in% par.table$name[i]]
}
write.csv(par.table, file='c:/SS/R/assessment_histories_2011plus_parameter_counts.csv')
