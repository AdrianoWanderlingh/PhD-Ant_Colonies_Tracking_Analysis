rm(list=ls())
gc()

########################################################################
#ASSIGN METADATA FROM METADATA-RICH FILES
# To ensure a correct orientation, AntsCreated files should be used as none of them have an assigned pose/orientation/etc.
# AntsCreated files can inherit Metadata info using CopyMetadata.R from unoriented files ("AntsCreated_DeathRecord_NoOrient.myrmidon") or from manually oriented files ("Oriented.myrmidon")


######load necessary libraries
library(FortMyrmidon) ####R bindings
library(Rcpp)
library(data.table)
library(circular)
library(R.utils)

## TRANSFORM MANUALLY ORIENTED FILES IN CLEAN FILES KEEPING Metadata ONLY
#COPY metadata info from man_oriented_myr_file to AntsCreated_myrmidon_file

#this script should be executed after the orientation of the files

##############################################################################################################################
##############################################################################################################################
### TO DO: 
# 1. COPY LOOP FROM AUTO-ORIENTATION LOOP TO ACCESS EACH FOLDER
# 2. AntsCreated_myrmidon_file SHOULD ACCESS ALL THE AUTO_ORIENT.MYRMIDON FILES
# 3. Metadata_myr_file SHOULD ACCESS THE CORRESPONDING FILESNAME (REP-PERIOD) FROM WHICH TO INHERIT THE METADATA

# temporary file/vars assignment
AntsCreated_myrmidon_file <- "/media/eg15396/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated.myrmidon"
Metadata_myr_file <- "/media/eg15396/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated_DeathRecord_NoOrient.myrmidon"
# AntsCreated_myrmidon_file <- "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated.myrmidon"
# Metadata_myr_file <- "/media/cf19810/DISK4/ADRIANO/EXPERIMENT_DATA/REP4/R4BP_20-03-21_AntsCreated_DeathRecord_NoOrient.myrmidon"


Metadata_exp <- fmExperimentOpen(Metadata_myr_file) 
Metadata_ants <- Metadata_exp$ants

### open antsCreated base file
AntsCreated  <- fmExperimentOpen(AntsCreated_myrmidon_file)
creat_exp_name <- unlist(strsplit(AntsCreated_myrmidon_file,split="/"))[length(unlist(strsplit(AntsCreated_myrmidon_file,split="/")))]
AntsCreated_ants <- AntsCreated$ants
#save
AntsCreated$save(paste0(sub("\\..*", "", AntsCreated_myrmidon_file),"_withMetaData.myrmidon")) # file now exists



### CHECK VISUALLY IF RETAG HAPPENED OR FROM METADATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## THEN EXPLORE THE STRUCURE, MAY BE NEEDED LOOP OVER LIST ELEMENTS

# check if retags happened:
for (ant in Metadata_ants){
 # individual  <- ant$ID
  print( ant$ID)
  print(ant$identifications)
}


  
for (ant in Metadata_ants){
  individual  <- ant$ID
#print(ant$identifications)

#AntsCreated_ants[[individual]]$identifications[[1]]$start <- Metadata_ants[[individual]]$identifications[[1]]$start
#if AntID matches
if (AntsCreated_ants[[individual]]$identifications[[1]]$targetAntID==Metadata_ants[[individual]]$identifications[[1]]$targetAntID){
  #remove identification
  AntsCreated$deleteIdentification(AntsCreated_ants[[individual]]$identifications[[1]])
  
  AntsCreated$addIdentification(antID= Metadata_ants[[individual]]$identifications[[1]]$targetAntID,
                              tagID= Metadata_ants[[individual]]$identifications[[1]]$tagValue,
                              start= Metadata_ants[[individual]]$identifications[[1]]$start,
                              end=   Metadata_ants[[individual]]$identifications[[1]]$end)
  
}
}



#CHECK copy
for (ant in AntsCreated_ants){
  print(ant$identifications)
  ant$identifications
}






list_keys <- list()
#assign metadata keys
for (KEY in   1:length(Metadata_exp$metaDataKeys)) {
  key <- names(Metadata_exp$metaDataKeys[KEY])
  defaultvalue <- unname(Metadata_exp$metaDataKeys[KEY][[1]])
  AntsCreated$setMetaDataKey(key,defaultvalue)
  #check
  #AntsCreated$metaDataKeys[KEY]
  list_keys <- c(list_keys,key)
}


for (ant in Metadata_ants){
  individual  <- ant$ID
  #extract metadata info per key
  for (METADATA_KEY in list_keys) {
    for (ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
      # AntsCreated$ants[[individual]]$deleteValue(key=METADATA_KEY,time=NA) #not possible to delete the first row of metadata as NA is not a valid time
      # if ("metadata value is equal in both files, skip assign") { #not possible as above
      if ( is.na(ant$getValues(METADATA_KEY)[ROW,"times"])){ # if there is an NA, set time as epoch start (fmTimeSinceEver doesn't work as it assigns currrent time)
        AntsCreated$ants[[individual]]$setValue(key=METADATA_KEY,value=ant$getValues(METADATA_KEY)[ROW,"values"],time=fmTimeCreate(offset = 0)) #fmTimeSinceEver()
      }else{
        AntsCreated$ants[[individual]]$setValue(METADATA_KEY,value=ant$getValues(METADATA_KEY)[ROW,"values"],time=fmTimeCreate(offset=ant$getValues(METADATA_KEY)[ROW,"times"]))
      }#assign value
     # }#skip if same
    }#ROW
  }#METADATA_KEY
}#ant

#save
AntsCreated$save(paste0(sub("\\..*", "", AntsCreated_myrmidon_file),"_withMetaData.myrmidon")) # file now exists

#check
# for (ant in AntsCreated_ants){
#   individual  <- ant$ID
#   #extract metadata info per key
#   for (METADATA_KEY in list_keys) {
#     cat(paste("ANT",individual,"\n","key:",METADATA_KEY,sep=" "))
#     
#     print(ant$getValues(METADATA_KEY))
#   }}

# cleaning: for value = base.value and time is 1970-01-01 01:00:00, delete event 
for (ant in AntsCreated_ants){
  individual  <- ant$ID
  #extract metadata info per key
  for (METADATA_KEY in list_keys) {
    for (ROW in 1:nrow(ant$getValues(METADATA_KEY))) {
      #if there is an assigned time
      if (!is.na(ant$getValues(METADATA_KEY)[ROW,"times"])){
      #if value is equal to base.value 
        if (ant$getValues(METADATA_KEY)[ROW,"values"]==AntsCreated$metaDataKeys[[METADATA_KEY]]) {
          AntsCreated$ants[[individual]]$deleteValue(METADATA_KEY,time=fmTimeCreate(offset = 0))
        }}}}}

cat("NOTE: \nAs it is not possible to assign a metadata value change without time info (NA), value for \"Is queen\" results as a timed change starting on epoch start (1st Jan 1970)")

#assign zones
list_zones <- NULL
for (ZONE in   1:length(Metadata_exp$spaces[[1]]$zones)) {
  zone <- Metadata_exp$spaces[[1]]$zones[[ZONE]]$name
  #defaultvalue <- unname(Metadata_exp$metaDataKeys[ZONE][[1]])
  AntsCreated$spaces[[1]]$createZone(zone)
  #check
  #AntsCreated$spaces[[1]]$zones
  list_zones <- c(list_zones,zone)
}


#Assign SHAPE to ZONE (geometry to the nest zone)
#extract metadata info per key
for (ZONE.1 in   1:length(AntsCreated$spaces[[1]]$zones)) {
#for (ZONE_KEY in list_zones) {
  zone_definition <- Metadata_exp$spaces[[1]]$zones[[ZONE.1]]$definitions
  #assign shapes
  AntsCreated$spaces[[1]]$zones[[ZONE.1]]$addDefinition(shapes=zone_definition[[1]][["shapes"]],start= fmTimeSinceEver(),end=fmTimeForever())
  }#ZONE.1


#save
AntsCreated$save(paste0(sub("\\..*", "", AntsCreated_myrmidon_file),"_withMetaData.myrmidon")) # file now exists



##################CHECK IDENTIFICATIONS AND OTHER METADATA INFO
  
  #make check on identifications to make sure there is a exact correspondance between ants in both Ants Created files and metadata provided files
  
  #pass all files through this system to make them all the same
  
