# loads the csv file that contains the catalaog of the gcm archive:
rm(list=ls())  #clear variables in R space

source('/Users/nicholasdepsky/Desktop/GCM_ClimTool_NUEVO/GCMClimTool_function_library_NUEVO.R')
setwd("/Users/nicholasdepsky/Desktop/GCM_ClimTool_NUEVO/Required_Files")
df_dumm = read.csv('Archive_list_batch_20161205_MagCauca.csv',sep=',',header=T) # path of the spreadsheet which controls which GCM models to be analyzed

df = subset(df_dumm,df_dumm$DoIt.==1)
df$BaseDay <- NA
df$BaseMonth <- NA
df$BaseYear <- NA
df$LeapFix <- NA

models = unique(df$Model)

fut_all <- unique(df$Experiment)
fut_all <- fut_all[fut_all != "Historical"]
var_all <- unique(df$Variable)

BootstrapRuns <- 1 # number of iterations to perform of the bootstrapping routine

# sets up the main variables of the tool:
mFolder = "/Users/nicholasdepsky/Desktop/GCM_ClimTool_NUEVO"   # main folder
enso_signal_file = "/Required_Files/OBSERVATIONS/ENSO_SIGNAL/ONI_SIGNAL.csv" # ENSO EL NINO-LA NINA CYCLE DATA
svd_analysis_file <- "/Required_Files/SVD/SVD_Table_Comp1_Pac_Atl.csv" # Pacific and Atlantic Ocean Sea-Surface Temp Oscillation Data

### ------------------------------
### Region (Magdalena-Cauca Basin, Colombia)
### ------------------------------
observation_files_lock = c(daily = "/Required_Files/OBSERVATIONS/SERIES_VARIABLE_DIARIO_v4.csv",
                           monthly = "/Required_Files/OBSERVATIONS/SERIES_VARIABLE_v4.csv",
                           stationsCatalog = "/Required_Files/OBSERVATIONS/CNE_IDEAM_V5_MAYO_6.csv")

shp_file_background = "/Required_Files/SHP/WEAP_Catchment_4.shp" # Basin Shapefile

region <- "Antioquia"  # name of your basin or region for use in output graphics

rFolder = "/Results/"  # Results Folder
save(df, file = paste0(mFolder,rFolder,"df.Rda"))

bBox = c(latMin = 5.5,latMax = 7.2,lonMin = -76.2,lonMax = -74.5)  # Bounding Box of Study Region

#-----------------------------------------------------------
# now loops through the catalog and reads the netcdfs:
#-----------------------------------------------------------

for (model in models) {

  df1 = subset(df, Model == model)
  exp_model = unique(df1$Experiment)
  
  for (experiment in exp_model){
    
    df2 = subset(df1, Experiment == experiment)
    var_exp = unique(df2$Variable)
    
    for (varName in var_exp) {
    
      df3 = subset(df2, Variable == varName)
      mv = df3$Moving.Year[1]
      varLabel = df3$Alias[1]
      gcmFolders = df3$File[1]
      observation_files = observation_files_lock 
      observation_files[1] = gsub("VARIABLE",varName,observation_files_lock[1])
      observation_files[2] = gsub("VARIABLE",varName,observation_files_lock[2])
      load(paste0(mFolder,rFolder,"df.Rda"))
      
      res <- try(read_GCM_NCDF(boundingBox = bBox,                       
                      movingDateStep = mv,                                                 
                      main_folder = mFolder,       
                      results_folder = rFolder,                                        
                      GCM_folders = gcmFolders,
                      variableName = varName,
                      variableLabel = varLabel,
                      shp_file = shp_file_background,                                   
                      obs_files = observation_files,
                      sta_bbox = "region", # can choose whether to have historical data chosen from the bbox of the selected GCM pixel ("gcm pixel"), or from the user-defined bbox ("region")
                      minObsYear = 1980,
                      minGCMYear = 2011)) 
    
      if(inherits(res, "try-error"))
      {
        try(dev.off())
      }
      load(paste0(mFolder,rFolder,"df.Rda"))
      print(varName)
    } # next variable
  } # next experiment 
} #next model  

#-----------------------------------------------------------------------------
#---------------------------------------------------------------------------
# loop for comparisons and bootstrapping:
for (model in models) {
  
  load(paste0(mFolder,rFolder,"df.Rda"))
  df1 = subset(df, Model == model)
  exp_model = unique(df1$Experiment)
  exp_model <- exp_model[exp_model != "Historical"] # excludes historical experiment from following functions to just have futures
  
  for(experiment in exp_model)
  {
    rcp_pr_ind <- which(df$Model == model & df$Variable == "pr" & df$Experiment == experiment) # just to find base date of a future scenario
    rcp_tas_ind <- which(df$Model == model & df$Variable == "tas" & df$Experiment == experiment)
    
    try(compare_GCM(refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind]),
                    main_folder = mFolder,       
                    results_folder = rFolder,
                    enso_signal_file = enso_signal_file,
                    svd_analysis_file = svd_analysis_file,
                    modelName = model,                     
                    futures = experiment,
                    varName = "pr",
                    varLabel = "Precipitation (mm)",
                    minObsYear = 1980,                     
                    maxGCMYear = 2040,
                    alignHistYears = TRUE))
    
    ## Optional: update of extremes Obs based (data_d) of GCM Results.
    ## The future reference year reprents the upper boundary year for which new extremes are derived
    try(extremePR_GPD(modelName = model, 
                      futures = experiment,
                      main_folder = mFolder,
                      results_folder = rFolder,
                      minObsYear = 1980,
                      minGCMYear = 2015,
                      maxGCMYear = 2040))
    
    for (n in 1:BootstrapRuns){
      try(knn_bootstrap(refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind]),
                        main_folder = mFolder,
                        results_folder = rFolder,
                        modelName = model,              #"CCSM4",
                        futures = experiment,
                        varName = "pr",
                        varLabel = "Precipitation [mm]",
                        minObsYear = 1980,
                        nearWindow = 15,   
                        minGCMYear = 2015,
                        maxGCMYear = 2040,
                        expNumber = n,
                        JPmode = "Window",
                        alignHistYears = TRUE,
                        HistRepro = FALSE))
      
      
      try(heatSignal(refDate = c(month = df$BaseMonth[rcp_tas_ind], day = df$BaseDay[rcp_tas_ind], year = df$BaseYear[rcp_tas_ind]), #rDate,
                     main_folder = mFolder,
                     results_folder = rFolder,
                     modelName = model,                     #"CESM1-CAM5",
                     futures = experiment,
                     varName = "tas",
                     varLabel = "Temperature [c]",
                     minObsYear = 1980,
                     maxObsYear = 2015,
                     minGCMYear = 2015,
                     maxGCMYear = 2040))
  
  
      try(sinteticSeries(refDate = c(month = df$BaseMonth[rcp_pr_ind], day = df$BaseDay[rcp_pr_ind], year = df$BaseYear[rcp_pr_ind]),
                         modelName = model,
                         futures = experiment,
                         main_folder = mFolder,
                         results_folder = rFolder,
                         expNumber = n,
                         fullObsCatalog = FALSE,
                         HistRepro = FALSE,
                         minObsYear = 1980,
                         maxObsYear = 2015))  
      
    } # next n ensemble member
  } # next experiment
} # next model