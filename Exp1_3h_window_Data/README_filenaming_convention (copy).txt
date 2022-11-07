Types of myrmidon files generated during the pre-processing:

RXtreatment: short for REPLICATE-PERIOD name
TSname: tracking system name (relevant only for "-base" files)

FILES:
RXtreatment_DATE_base : first file created in FortStudio, to load the experiment folders
RXtreatment_DATE_AntsCreated : file were Ants are Created, produced by Define_Ant_Identifications_standalone.R
RXtreatment_DATE_AntsCreated_DeathRecord_NoOrient : non oriented files with Metadata information
RXtreatment_DATE_ManOriented : manually oriented files with Metadata information
RXtreatment_DATE_ManOriented_TSname-base : reference files for ant length by Tracking System, base for DataPrep1_Clone-capsule-manual-to-manual_v082.R
RXtreatment_DATE_AntsCreated_AutoOriented : automatically oriented files produced by DataPrep2_auto-orientation-loop_v082.R
RXtreatment_DATE_AntsCreated_AutoOriented_withMetaData: files to which metadata are assigned, produced by DataPrep3_CopyMetadata_v082.R

more information here: https://docs.google.com/document/d/1NXH4hZ97UtmRAYg50thtf_ObNDTt6loE/edit?usp=sharing&ouid=104843790314847559053&rtpof=true&sd=true
