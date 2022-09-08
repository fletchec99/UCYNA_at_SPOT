#Parse through SPOT CTD data that have already had temp, O2, salinity, and depth extracted
##And that have problematic symbols removed from headers
#Non-for-loop version 04.01.2021

getwd()
setwd("SPOT_PROCESSED_CTD_Data/formatted_CTD_files_all_dates")
list.files(pattern="_temp_df.tsv")


data=read.table("d20000918C_temp_df.tsv", header=T, stringsAsFactors = F, sep="\t")
names(data)

#Format name --> Sampling date
if(data$Name[1]!=unique(data$Name)){
  stop("There is something wrong with file: it has multiple dates.n", call=FALSE)
}

a=strsplit(data$Name[1], split="", fixed=T)
#Year
data$Year=paste(a[[1]][2], a[[1]][3], a[[1]][4], a[[1]][5], sep="")
#Month
if(as.numeric(paste(a[[1]][6], a[[1]][7], sep=""))<10){
  data$Month=a[[1]][7]
}else{
  data$Month=paste(a[[1]][6], a[[1]][7], sep="")
}
#Day
data$Date=paste(a[[1]][8], a[[1]][9], sep="")

#Try to match SPOT DCM
DCM=liv_data$Depth_m[grep(paste(data$Year[1], data$Month[1], data$Date[1], "DCM", sep="."), liv_data$SampleID)]

#Now find values for each SPOT depth 
SPOT_depths=c(5.00, DCM, 150.00, 500.00, 885.00)
include_rows=c()
ERROR=c()
for(depth in SPOT_depths){
  #print(depth)
  ROW=which(data$Depth>(depth-0.5)&data$Depth<(depth+0.5))[1]
  #print(ROW)
  if(is.na(ROW)==TRUE){
    ERROR=c(ERROR, paste(data$Name[1], "is missing data from", depth, "m", sep=" "))
  }else{
    include_rows=c(include_rows, ROW)
  }
}

#Subset data
data_subs=data[include_rows,]
#dim(data_subs)

#Write out data
write.table(x=data_subs, file=paste(data$Name[1], "temp_subsetted.tsv", sep="_"), row.names=F, sep="\t", quote=F)

#Also write out ERROR (date + missing depths)
write.table(x=ERROR, file=paste(data$Name[1], "error_subsetted.txt", sep="_"), quote=F)
