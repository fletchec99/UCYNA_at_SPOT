#Format SPOT abiotic data
#Format the following data: SPOT data from Liv (1_...UnderProtection in /OriginalFiles), NOAA BEUTI upwelling data, NOAA CUTI upwelling data, NOAA Bakun upwelling data, NOAA MEI data
#Ver. 02.2021

getwd()
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

#PART I: Chemical and biological data (nutrients, bacterial prod.)-----

#1. Liv's SPOT Data-----
#Already removed some columns and rows that were hashed out using LibreOffice/ Excel
#Read it in 
Liv_data <- read.csv("abio_data_unformatted/SPOT_Nutrients_Fuhrman_Wrigley_2000_2019_Formatted.csv", header=T, stringsAsFactors = F)
dim(Liv_data) #890 9
names(Liv_data)

#Adjust column Depth_m to match DepthBin
#Some days 5m is zero, some days 5, some days "?"
length(grep("5m", Liv_data$DepthBin)) #199
Liv_data$Depth_m[grep("5m", Liv_data$DepthBin)]=5
length(grep("DCM", Liv_data$DepthBin)) #179
Liv_data$Depth_m[grep("DCM", Liv_data$DepthBin)] #This is gonna take some doing 
length(grep("150", Liv_data$DepthBin)) #177
Liv_data$Depth_m[grep("150m", Liv_data$DepthBin)]=150
length(grep("500", Liv_data$DepthBin)) #182
Liv_data$Depth_m[grep("500", Liv_data$DepthBin)]=500
length(grep("890", Liv_data$DepthBin)) #153 why?
Liv_data$Depth_m[grep("890", Liv_data$DepthBin)]=890
#Write this back out and manually fix DCM depths

(x=Liv_data, file="abio_data_unformatted/SPOT_nutrients_formatted_02.13.2021.csv", row.names=F, quote=F, sep=",")
#Write back in 
liv_data=read.csv("abio_data_unformatted/B_SPOT_nutrients_formatted_02.13.2021.csv", header=T, stringsAsFactors = F)
rm(Liv_data)

#Paste together year and month
dim(liv_data)
names(liv_data)
#liv_data$Year_Month=paste(liv_data$Year, liv_data$Month, sep="_")
head(liv_data$Year_Month)
#Also unique sample IDs with year, date, month, and depthBin
liv_data$SampleID=paste(liv_data$Year, liv_data$Month, liv_data$Day, liv_data$DepthBin, sep=".")
head(liv_data$SampleID)
length(unique(liv_data$SampleID)) #890

#Set up "Sample ID" column #Oh forgot I did that
#liv_data$SampleID=paste(liv_data$Year, liv_data$Month, liv_data$Day, liv_data$DepthBin, sep=".")

liv_data$Year.Mo.Dep=paste(liv_data$Year, liv_data$Month, liv_data$DepthBin, sep=".")
length(liv_data$Year.Mo.Dep) #890
length(unique(liv_data$Year.Mo.Dep)) #860 #Crap, there are 30 duplicates
liv_data$SampleID[which(duplicated(liv_data$Year.Mo.Dep))] #Shoot, there are 30 days where sampling was done twice in the same month, for some reason 
#Remove duplicate rows?? 

#Remove duplicated samples: double check their sample ID against the SPOT 16S data
proks_data=read.table("ModifiedFiles/2.SPOT_16S_no_chloro_counts.tsv", header=T, stringsAsFactors = F, sep=" ") #For some 
dim(proks_data)
colnames(proks_data)[1:20]
DATES=colnames(proks_data)[2:ncol(proks_data)]
DATES 
#Remove proks_data as a variable, because I think it's slowing R down
rm(list=ls(pattern="proks_data"))

#Gonna have a lot of missing samples if we go through all the rows in liv_data -- just do the duplicates
du_rows=c()
for(i in liv_data$Year.Mo.Dep[which(duplicated(liv_data$Year.Mo.Dep))]){
  print(i)
  du_rows=c(du_rows, grep(i, liv_data$Year.Mo.Dep))
  print(liv_data$SampleID[grep(i, liv_data$Year.Mo.Dep)])
}
length(du_rows) #72
length(unique(du_rows)) #54

#Grep liv_data rows in SPOT DATES variable
#If the row is missing from SPOT data, remove the ID of that row from liv_data
exclude_rows=c()
for(i in du_rows){
  #print(liv_data$SampleID[i])
  if(liv_data$Month[i]<10){
    ID=paste(liv_data$Year[i], paste("0", liv_data$Month[i], sep=""), liv_data$Day[i], liv_data$DepthBin[i], sep=".")
    no=grep(ID, DATES)
    #print(no)
    if(length(no)<1){
      exclude_rows=c(exclude_rows, i)
    }
  }else{
    ID=liv_data$SampleID[i]
    no=grep(ID, DATES)
    #print(no)
    if(length(no)<1){
      exclude_rows=c(exclude_rows, i)
    }
  }
}
exclude_rows
length(unique(exclude_rows)) #26

#Manually pick which rows we want to keep 
liv_data$SampleID[exclude_rows]
liv_data$SampleID[grep("2000.10", liv_data$SampleID)] #keep 2000.10.9, exclude 10.27
exclude_rows=exclude_rows[-grep("2000.10.9", liv_data$SampleID[exclude_rows])]
liv_data$SampleID[grep("2012.7", liv_data$SampleID)] #Keep 2012.7.9, get rid of 2012.7.21?
liv_data$SampleID[grep("2014.7", liv_data$SampleID)] #Get rid of both the 17th and the 23rd
liv_data$SampleID[grep("2014.8", liv_data$SampleID)] #Get rid of the 14th
liv_data$SampleID[grep("2015.10", liv_data$SampleID)] #Get rid of the 7th

#Subset liv_data accordingly 
dim(liv_data)
liv_subs=liv_data[-exclude_rows,]
dim(liv_subs)

#Check that there are no more duplicates 
liv_subs$SampleID[which(duplicated(liv_subs$Year.Mo.Dep))]
sort(liv_subs$SampleID[grep("2012.6", liv_subs$SampleID)]) #Pick 6.20.2012, not 6.13.2012
liv_subs=liv_subs[-grep("2012.6.12|2012.6.13", liv_subs$SampleID),]
liv_subs$SampleID[grep("2013.12", liv_subs$SampleID)] #Get rid of 12.19.2013
liv_subs$SampleID[grep("2013.12.19", liv_subs$SampleID)]
liv_subs=liv_subs[-grep("2013.12.19", liv_subs$SampleID),]

which(duplicated(liv_subs$Year.Mo.Dep)) #NONE, F YEAH
#Now skip down to #3, reading in the Leu data

#2. SPOT CTD data-----
#First step: Read in the 175 files, grep only the columns I want, write out files
directory <- getwd()
filenames=system(paste('ls ', directory, '/SPOT_PROCESSED_CTD_Data/CTD_files', sep=''), intern=TRUE)

for(i in c(grep("asc.txt", filenames))){
  print(filenames[i])
  temp=read.table(file=paste("SPOT_PROCESSED_CTD_Data/CTD_files/", filenames[i], sep=""), header=T, stringsAsFactors = F, sep="\t")
  #Make new df
  temp_df=as.data.frame(matrix(nrow=nrow(temp), ncol=5))
  colnames(temp_df)=c("Name", "Depth", "Temp", "Sal", "Oxy")
  #Fill it with variables by grepping out columns corresponding to the variables I want: Name, PrDM/DepSM/Depth, T090/CTDTMP, Sal00/Sal11, OxML, put them in a new df
  temp_df$Name=strsplit(filenames[i], split=".", fixed=T)[[1]][1]
  temp_df$Depth=temp[,grep("PrDM|Depth|DepSM", colnames(temp), ignore.case=T)[1]] #Just take whichever comes first
  temp_df$Temp=temp[,grep("T090|CTDTMP", colnames(temp), ignore.case=T)[1]]
  temp_df$Sal=temp[,grep("Sal00|CTDSAL", colnames(temp), ignore.case=T)[1]]
  temp_df$Oxy=temp[,grep("Sbeox0mL|Oxsol|Oxsat", colnames(temp), ignore.case=T)[1]]
  #Order columns 
  temp_df=temp_df[,order(colnames(temp_df))]
  #Write out 
  
  (x=temp_df, file=(paste("SPOT_PROCESSED_CTD_Data/temp_CTD_files/", temp_df$Name[1], "_temp_df.tsv", sep="")), quote=F, row.names=F, sep="\t")
}
#OMG it worked!! :D :D :D 

#Write a separate R script to loop through all abbreviated CTD files, labeled *_temp_df.tsv, and pull out SPOT depths within a meter of each depth. (See SPOT_PROCESSED_CTD_Data/formatted_CTD_files_all_dates"01_parsing_CTD_data_v_04.02.2021.R" & "02_combine_subsetted_data.sh")
#Concatenate all subsetted files into one file, concat_subsetted_CTD_files_04.02.21.tsv
#Copy/ paste this file into "abio_data_unformatted" and read it in!! (Also copy/ paste error file, don't read that in)

#Removing all other variables called "CTD_*"
rm(list=ls(pattern="CTD"))
ls(pattern="CTD")

#Read in concatenated file
CTD_data=read.table("abio_data_unformatted/concat_subsetted_CTD_files_04.02.21.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(CTD_data) #955 8
names(CTD_data)

#Get rid of 174 times that the rownames appear
length(grep("Depth", CTD_data$Depth)) #174
CTD_data=CTD_data[-grep("Depth", CTD_data$Depth),]
dim(CTD_data) #781 8
names(CTD_data)
grep("Depth", CTD_data$Depth) #integer(0)

#Convert Depth, Oxy, Sal, Temp, Year, Month, Date to be numeric, not characters
CTD_data$Depth=as.numeric(CTD_data$Depth)
head(CTD_data$Depth)
CTD_data$Oxy=as.numeric(CTD_data$Oxy)
CTD_data$Sal=as.numeric(CTD_data$Sal)
head(CTD_data$Sal)
CTD_data$Temp=as.numeric(CTD_data$Temp)
head(CTD_data$Temp)
CTD_data$Year=as.numeric(CTD_data$Year)
tail(CTD_data$Year)
CTD_data$Month=as.numeric(CTD_data$Month)
head(CTD_data$Month)
CTD_data$Date=as.numeric(CTD_data$Date)
head(CTD_data$Date)

#Write in Depth_bin variable 
CTD_data$DepthBin="NA"
length(which(CTD_data$Depth<6)) #171
#174, booya #Wait, actually there are a few samples at depth between 5-10m, just how shallow is the DCM?? 
sort(liv_data$Depth_m[which(liv_data$DepthBin=="DCM")]) #Ok so there are some dates where DCM was 7m. Hmm. 7m. Anything less than 6m shall be assigned to 5m depth bin 
CTD_data$DepthBin[which(CTD_data$Depth<6)]="5m"

length(which(CTD_data$Depth<145 & CTD_data$Depth>6)) #161 
CTD_data$DepthBin[which(CTD_data$Depth<145 & CTD_data$Depth>6)]="DCM"

length(which(CTD_data$Depth>145 & CTD_data$Depth<495)) #174
CTD_data$DepthBin[which(CTD_data$Depth>145 & CTD_data$Depth<495)]="150m"

length(which(CTD_data$Depth>495&CTD_data$Depth<880)) #159
CTD_data$DepthBin[which(CTD_data$Depth>495&CTD_data$Depth<880)]="500m"

length(which(CTD_data$Depth>880)) #116
CTD_data$DepthBin[which(CTD_data$Depth>880)]="890m"

#Double check number of rows 
length(which(CTD_data$Depth<6)) + length(which(CTD_data$Depth<145 & CTD_data$Depth>6)) + length(which(CTD_data$Depth>145 & CTD_data$Depth<495)) + length(which(CTD_data$Depth>495&CTD_data$Depth<880)) + length(which(CTD_data$Depth>880)) == nrow(CTD_data) #TRUE

#Create SampleID
head(liv_data$SampleID, n=15) #I guess August 2020 is missing from CTD data
CTD_data$SampleID=paste(CTD_data$Year, CTD_data$Month, CTD_data$Date, CTD_data$DepthBin, sep=".")
head(CTD_data$SampleID)

#Is SampleID unique? 
length(unique(CTD_data$SampleID)) #732
length(CTD_data$SampleID) #781 #Uh oh 
CTD_data$SampleID[which(duplicated(CTD_data$SampleID))]
#Just look at the first one
grep(CTD_data$SampleID[which(duplicated(CTD_data$SampleID))][1], CTD_data$SampleID)
CTD_data[grep(CTD_data$SampleID[which(duplicated(CTD_data$SampleID))][1], CTD_data$SampleID),] #Identical
CTD_data[grep(CTD_data$SampleID[which(duplicated(CTD_data$SampleID))][49], CTD_data$SampleID),] #Not identical, A and B casts from the same day 
for(i in which(duplicated(CTD_data$SampleID))){
  print(grep(CTD_data$SampleID[i], CTD_data$SampleID))
  #print(CTD_data$Name[grep(CTD_data$SampleID[i], CTD_data$SampleID)])
}
#Ok, take the last one out of each element #also set up a vector of rows to include to be sure we keep that one
exclude_rows=c()
include_rows=c()
for(i in which(duplicated(CTD_data$SampleID))){
  exclude_rows=c(exclude_rows, grep(CTD_data$SampleID[i], CTD_data$SampleID)[-length(grep(CTD_data$SampleID[i], CTD_data$SampleID))])
  include_rows=c(include_rows, grep(CTD_data$SampleID[i], CTD_data$SampleID)[length(grep(CTD_data$SampleID[i], CTD_data$SampleID))])
}
length(unique(exclude_rows))
sort(exclude_rows)
sort(include_rows)

#Subset CTD data to exclude the rows that are duplicated
CTD_subs=CTD_data[-exclude_rows,]
dim(CTD_subs) #732 10
dim(CTD_data)
nrow(CTD_data)-nrow(CTD_subs)==length(unique(exclude_rows)) #TRUE
length(unique(CTD_subs$SampleID)) #732

#Write dummy variables into liv_data for oxygen, temp, salinity 
names(liv_data)
liv_data$Oxyen="NA"
liv_data$Salinity="NA"
liv_data$Temperature="NA"

#Write variables from CTD_subs into Liv_data 
#Samples missing FROM liv_data
missing_samples=c()
for(i in c(1:nrow(CTD_subs))){
  a=grep(CTD_subs$SampleID[i], liv_data$SampleID)
  if(length(a)<1){
    missing_samples=c(missing_samples, CTD_subs$SampleID[i])
  }#else{
    #print(a)
    #liv_data[a,c(11:13)]=CTD_subs[i,c(3:5)]
  #}
}
missing_samples #Most of the missing samples are from 150m
length(missing_samples)
length(grep("150m", missing_samples)) #Yep, 11/15 missing samples

#Go through those missing samples and check them out
for(i in missing_samples){
  print(i)
  a=strsplit(i, split=".", fixed=T)
  b=paste(a[[1]][1], a[[1]][2], a[[1]][3], sep=".")
  print(liv_data$SampleID[grep(b, liv_data$SampleID)])
}
#All data from 9/12/2018 is missing
#Guess the others are too

#Convert Liv data to be numeric in the oxy, temp, sal departments
head(liv_data$Oxyen)
names(liv_data)[grep("Oxy", names(liv_data))]="Oxygen" #Forgot the G, oops! 
liv_data$Oxygen=as.numeric(liv_data$Oxygen)
head(liv_data$Oxygen)
liv_data$Salinity=as.numeric(liv_data$Salinity)
head(liv_data$Salinity)
liv_data$Temperature=as.numeric(liv_data$Temperature)
head(liv_data$Temperature)

length(which(is.na(liv_data$Temperature)))/nrow(liv_data) #Only 20% of temp data is missing
length(which(is.na(liv_data$Temperature[grep("5m", liv_data$DepthBin)])))/length(grep("5m", liv_data$DepthBin)) #22% of temp data at 5m is missing 
length(which(is.na(liv_data$Temperature)|is.na(liv_data$Oxygen)|is.na(liv_data$Salinity)))/nrow(liv_data) #173
#Overall, 20% of data missing - is it the same dates?? 
which(is.na(liv_data$Temperature) & liv_data$Oxygen>1) #yep, integer(0)

#COULD try to go to cruise logs and fill in manually 

#3. Bacterial production data------
#a. Leucine data-----
leu_data=read.csv("abio_data_unformatted/SPOT_bacterial_production_Leu_2000_2021_04.05.2021.csv", header=T, stringsAsFactors = F)
dim(leu_data) #974 4
head(leu_data)

#First, split apart the mo.yr column to generate month, year, date 
strsplit(head(leu_data$mo.yr), split="/", fixed=T)
leu_data$Month=01
leu_data$Day=01
leu_data$Year=2000

#Overwrite with real data
for(i in c(1:nrow(leu_data))){
  leu_data$Month[i]=strsplit(leu_data$mo.yr[i], split="/", fixed=T)[[1]][1]
  leu_data$Day[i]=strsplit(leu_data$mo.yr[i], split="/", fixed=T)[[1]][2]
  leu_data$Year[i]=strsplit(leu_data$mo.yr[i], split="/", fixed=T)[[1]][3]
}
sort(unique(leu_data$Year))
leu_data$mo.yr[which(is.na(leu_data$Year))] #WTF is 11/05? Eff that! 
leu_data$mo.yr[which(leu_data$Year==19)] #January through April of 2019
leu_data$Year[which(leu_data$Year==19)]=2019

#Fix problematic months
leu_data$Month[which(leu_data$Month=="01")]=1
leu_data$Month[which(leu_data$Month=="02")]=2
leu_data$Month[which(leu_data$Month=="03")]=3
leu_data$Month[which(leu_data$Month=="04")]=4
sort(unique(leu_data$Month))

#Then paste them together with DepthCategory to create unique sample IDs
leu_data$SampleID=paste(leu_data$Year, leu_data$Month, leu_data$Day, leu_data$DepthCategory, sep=".")
length(leu_data$SampleID) #969
length(unique(leu_data$SampleID)) #969
which(duplicated(leu_data$SampleID))
leu_data[grep(leu_data$SampleID[which(duplicated(leu_data$SampleID))][1], leu_data$SampleID),] #Can remove those data
leu_data[grep(leu_data$SampleID[which(duplicated(leu_data$SampleID))][2], leu_data$SampleID),] #And those too
#And year.month.depth



#Subsetting data to remove problematic samples 
leu_data=leu_data[-c(grep(leu_data$SampleID[which(duplicated(leu_data$SampleID))][1], leu_data$SampleID), grep(leu_data$SampleID[which(duplicated(leu_data$SampleID))][2], leu_data$SampleID), which(is.na(leu_data$Year))),]
dim(leu_data) #969 rows
length(unique(leu_data$SampleID)) #Why is this much longer than Liv's data? #We are including months that Liv doesn't have data for

#Sort the leu_data by year (I think it's currently ordered by least -> most productive)
dim(leu_data)
leu_data=leu_data[order(leu_data$SampleID),]
dim(leu_data)
tail(leu_data$SampleID)

#Remove months after end of Liv's dataset
tail(sort(liv_data$SampleID)) #Last sampling date is 2019.2.13
leu_data$SampleID[which(leu_data$Year==2019 & leu_data$Month>2)[1]] 

#Subset data
leu_subs=leu_data[c(1:(which(leu_data$Year==2019 & leu_data$Month>2)[1]-1)),]
dim(leu_subs) #886 9 
tail(leu_subs$SampleID)

#Figure out which rows are duplicates, and remove them
length(unique(leu_subs$SampleID))
nrow(leu_subs)
which(duplicated(leu_subs)) #None! :) 

#Create year.month.depth
leu_subs$Year.Mo.Depth=paste(leu_subs$Year, leu_subs$Month, leu_subs$DepthCategory, sep=".")
length(unique(leu_subs$Year.Mo.Depth)) #870 #16 duplicates
nrow(leu_subs)

#Set up dummy variables for leucine production in Liv's data  
names(liv_subs)
liv_subs$BactProd_Leu="NA"
#I THINK this is how you set it to be NA, not "NA"
liv_subs$BactProd_Leu=as.numeric(liv_subs$BactProd_Leu)
head(liv_subs$BactProd_Leu)
length(which(is.na(liv_subs$BactProd_Leu))) #860
nrow(liv_subs)

#Write leu data in based on year.month.depth -- days of sampling may vary 
missing_samples=c() #Samples missing FROM Liv's data
for(i in c(1:nrow(leu_subs))){
  a=grep(leu_subs$Year.Mo.Depth[i], liv_subs$Year.Mo.Dep)
  #print(a)
  if(length(a)<1){
    missing_samples=c(missing_samples, leu_subs$Year.Mo.Depth[i])
  }#else{
    #liv_subs$BactProd_Leu[a]=leu_data$Average_cell.ml.day[i]
  #}
}
missing_samples=sort(missing_samples)
length(missing_samples) #66
nrow(leu_subs)-nrow(liv_subs) #26
#So there are 40 samples that are truly missing
length(which(is.na(liv_subs$BactProd_Leu))) #47
liv_subs$SampleID[which(is.na(liv_subs$BactProd_Leu))]

#Check that they are really missing
for(i in missing_samples){
  a=strsplit(i, split=".", fixed=T)
  b=paste(a[[1]][1], a[[1]][2], sep=".")
  print(i)
  print(liv_subs$Year.Mo.Dep[grep(b, liv_subs$Year.Mo.Dep)])
}
#Some of these missing samples are from too late in 2019

#Fix: "2005.7.883", "2006.8.887"
unique(liv_subs$DepthBin)
unique(liv_subs$Depth_m)
grep("887|883", liv_subs$Year.Mo.Dep) #IDK what that's about.
#OKAY, it looks like they really are missing!! :) Moving on to Thy!! 

#b. Thymidine data-----
#Read in data
thy_data=read.csv("abio_data_unformatted/SPOT_bacterial_production_Thy_2000_2021_04.06.21.csv", header=T, stringsAsFactors = F)
dim(thy_data)
names(thy_data)
head(thy_data)

#Break apart that first column to create year, month, date
#Oh god, do all the months have leading zeros? #Nope, it's a mix 8)
strsplit(thy_data$mo.yr, split="/", fixed=T)
thy_data$Year=2000
thy_data$Month=1
thy_data$Date=1

#Fixing error message from for-loop real quick: 
thy_data$mo.yr[440]
thy_data=thy_data[-440,]
dim(thy_data)

#Overwrite with actual sampling dates
for(i in c(1:nrow(thy_data))){
  #print(thy_data$mo.yr[i])
  a=strsplit(thy_data$mo.yr[i], split="/", fixed=T)
  if(a[[1]][3]<2000){
    thy_data$Year[i]=paste(20,a[[1]][3], sep="")
  }else{
    thy_data$Year[i]=a[[1]][3]
  }
  thy_data$Date[i]=a[[1]][2]
  thy_data$Month[i]=a[[1]][1]
}

#Fix problematic dates
sort(unique(thy_data$Year))
thy_data$Year[which(thy_data$Year==21)]=2021
sort(unique(thy_data$Month))
thy_data$Month[which(thy_data$Month=="01")]=1
thy_data$Month[which(thy_data$Month=="02")]=2
thy_data$Month[which(thy_data$Month=="03")]=3
thy_data$Month[which(thy_data$Month=="04")]=4
thy_data$Month[which(thy_data$Month=="05")]=5
thy_data$Month[which(thy_data$Month=="06")]=6
thy_data$Month[which(thy_data$Month=="07")]=7
thy_data$Month[which(thy_data$Month=="08")]=8
thy_data$Month[which(thy_data$Month=="09")]=9
sort(unique(thy_data$Month))

#Paste into sample ID and year.mo.depth
unique(thy_data$DepthBIN)
#Remove samples from "0", "7", and "Bucket"
thy_data$DepthBIN[grep("5|DCM|150|500|890", thy_data$DepthBIN, invert=T)]
dim(thy_data) #590 7
thy_data=thy_data[-grep("5|DCM|150|500|890", thy_data$DepthBIN, invert=T),]
dim(thy_data) #575 7
unique(thy_data$DepthBIN)

thy_data$SampleID=paste(thy_data$Year, thy_data$Month, thy_data$Date, thy_data$DepthBIN, sep=".")
nrow(thy_data) #575
length(unique(thy_data$SampleID)) #575
thy_data$Year.Mo.Dep=paste(thy_data$Year, thy_data$Month, thy_data$DepthBIN, sep=".")
length(unique(thy_data$Year.Mo.Dep)) #563 - 12 duplicate months

#Sort by date
dim(thy_data)
thy_data=thy_data[order(thy_data$SampleID),]
dim(thy_data)

#Remove dates after end of liv_data
tail(liv_data$SampleID)
dim(thy_data) #575 9
thy_data=thy_data[c(1:(which(thy_data$Year==2019 & thy_data$Month>2)[1]-1)),]
dim(thy_data) #531 9
tail(thy_data$SampleID)
grep("2020", thy_data$SampleID)

#Remove duplicates
length(unique(thy_data$SampleID)) #531
length(unique(thy_data$Year.Mo.Dep)) #519 #12 duplicates

#Figure out which rows are duplicates
du_rows=c()
for(i in thy_data$SampleID[which(duplicated(thy_data$Year.Mo.Dep))]){
  #print(i)
  a=strsplit(i, fixed=T, split=".")
  b=paste(a[[1]][1], a[[1]][2], a[[1]][4], sep=".")
  du_rows=c(du_rows, grep(b, thy_data$Year.Mo.Dep))
}
du_rows
length(unique(du_rows)) #24, as expected
du_rows=unique(du_rows)

#Figure out which of the duplicates has 16S data
thy_data$SampleID[du_rows]

#Actually, I can just subset based on removing duplicates from Liv's data (see above)
grep("2000.10.27", thy_data$SampleID[du_rows]) #Get rid of 2000.10.27
DATES[grep("2012", DATES)] #Executive decision, keep 5.30.2012, remove 5.3.2012
grep("2012.5.3", thy_data$SampleID[du_rows])
liv_subs$SampleID[grep("2012.6", liv_subs$SampleID)] #Keep 6.20.2012, remove 6.13.2012
grep("2012.6.13", thy_data$SampleID[du_rows])
liv_subs$SampleID[grep("2012.9", liv_subs$SampleID)] #Keep 9.28.2012, remove 9.21.2012
grep("2012.9.21", thy_data$SampleID[du_rows])
liv_subs$SampleID[grep("2014.7", liv_subs$SampleID)] #Keep the 17th, get rid of the 24th
grep("2014.7.24", thy_data$SampleID[du_rows])
liv_subs$SampleID[grep("2015.10", liv_subs$SampleID)] #Keep 10.20.2015, Zachary's birthday :)
grep("2015.10.07", thy_data$SampleID[du_rows])
liv_subs$SampleID[grep("2016.5", liv_subs$SampleID)] #LAST ONE! Keep 5.18.2016, remove 5.16.2016
grep("2016.5.16", thy_data$SampleID[du_rows])

exclude_rows=du_rows[c(grep("2000.10.27", thy_data$SampleID[du_rows]), 
              grep("2012.5.3", thy_data$SampleID[du_rows]), 
              grep("2012.6.13", thy_data$SampleID[du_rows]),
              grep("2012.9.21", thy_data$SampleID[du_rows]),
              grep("2014.7.24", thy_data$SampleID[du_rows]),
              grep("2015.10.07", thy_data$SampleID[du_rows]),
              grep("2016.5.16", thy_data$SampleID[du_rows]))]
du_rows
exclude_rows
length(du_rows)
length(exclude_rows)
thy_data$SampleID[du_rows]
thy_data$SampleID[exclude_rows] #All the dates I want to get rid of 

#Subset thy_data
dim(thy_data) #531 9
thy_subs=thy_data[-exclude_rows,]
dim(thy_subs) #519 9

#Make sure rows really are all unique
nrow(thy_subs)
length(unique(thy_subs$Year.Mo.Dep)) #519 #good!

#Write in dummy variable in Liv's data
liv_subs$BactProd_Thy="NA"
liv_subs$BactProd_Thy=as.numeric(liv_subs$BactProd_Thy)

#Overwrite with actual data, based on Year.Mo.Dep
missing_thy=c() #thy data missing FROM liv's samples
for(i in c(1:nrow(thy_subs))){
  #print(thy_subs$Year.Mo.Dep[i])
  a=grep(thy_subs$Year.Mo.Dep[i], liv_subs$Year.Mo.Dep)
  #print(a)
  if(length(a)<1){
    missing_thy=c(missing_thy, thy_subs$Year.Mo.Dep[i])
  }else{
    liv_subs$BactProd_Thy[a]=thy_subs$Average_cells.ml.day[i]
  }
}
missing_thy=sort(missing_thy)
length(missing_thy) #36 #ok

length(which(is.na(liv_subs$BactProd_Thy))) #246
nrow(liv_subs) #860
nrow(liv_subs)-length(which(is.na(liv_subs$BactProd_Thy))) #614 rows have data?? 
length(which(liv_subs$BactProd_Thy>0)) #614 #But there are 100 fewer rows in thy data?
length(which(duplicated(liv_subs$BactProd_Thy))) #453
#Wait
length(which(liv_subs$BactProd_Thy>0 & liv_subs$DepthBin=="500m")) #168 #Boom
#All the 5m samples are hitting the 500m samples too
#Change "5" to "5m" on DepthBin, re-paste Year.Mo.Dep and re-run floop
unique(thy_subs$DepthBIN)
#thy_data$DepthBIN[which(thy_subs$DepthBIN==5)]="5m" #Oops!
thy_subs$DepthBIN[which(thy_subs$DepthBIN==5)]="5m"
thy_subs$Year.Mo.Dep=paste(thy_subs$Year, thy_subs$Month, thy_subs$DepthBIN, sep=".")
unique(thy_subs$DepthBIN)
length(unique(thy_subs$Year.Mo.Dep)) #519 #good

#Double-check again 
length(which(is.na(liv_subs$BactProd_Thy))) #377
nrow(liv_subs)-length(which(is.na(liv_subs$BactProd_Thy))) #483
nrow(thy_subs)-length(missing_thy) #483 
(nrow(liv_subs)-length(which(is.na(liv_subs$BactProd_Thy))))+length(missing_thy)==nrow(thy_subs)
#Length of liv_data that has bacterial production available + length missing samples

#Which of the samples missing from thy data are also missing from leu data? 
for(i in missing_thy){
  print(i)
  print(grep(i, missing_samples))
}
#Most of them! 

#4. Hammond data------
hammond_data=read.csv("abio_data_unformatted/Hammond_WIES_abio_data_formatted_01.22.21.csv", header=T, stringsAsFactors = F)
dim(hammond_data) #1053 14
head(hammond_data)
names(hammond_data)

#Figure out which depths correspond to SPOT depths
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==12)] #Take bottle 12 to be 5m
#DCM data not sampled -> need to linearly interpolate
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==5)] #Take bottle 5 to be 150m #Ranges from 100-150m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==3)] #Take bottle 3 to be 500m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==1)] #Take bottle 1 to be 890m

#depth.m and depth are the same except for these samples at 885m
hammond_data[which(hammond_data$Depth!=hammond_data$Depth.m),grep("Depth", names(hammond_data))]
#Ok to use interchangeably

#Samples collected in 2017 and 2016 are missing depth data
hammond_data$Year[which(hammond_data$Depth=="")]
length(which(hammond_data$Depth.m=="")) #137 samples are missing depth
nrow(hammond_data) #1053 
sort(hammond_data$CTD.Bottle.No[which(hammond_data$Depth.m=="")])

#Fix missing depths
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==1 & hammond_data$Depth.m==2)]=885 #Bottle1=885m #Accidentally set to 2m :( 
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==2 & hammond_data$Depth.m=="")]=750 #Bottle2=750m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==3 & hammond_data$Depth.m=="")]=500 #Bottle3=500m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==4 & hammond_data$Depth.m=="")]=250 #Bottle4=250m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==5 & hammond_data$Depth.m=="")]=150 #Bottle5=150m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==6 & hammond_data$Depth.m=="")]=60 #Bottle 6=60m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==7 & hammond_data$Depth.m=="")]=50 #Bottle 7=50m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==8 & hammond_data$Depth.m=="")]=40 #Bottle 8=40m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==9 & hammond_data$Depth.m=="")]=30 #Bottle 9=30m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==10 & hammond_data$Depth.m=="")]=20 #Bottle10=20m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==11 & hammond_data$Depth.m=="")]=10 #Bottle11=10m
hammond_data$Depth.m[which(hammond_data$CTD.Bottle.No==12 & hammond_data$Depth.m=="")]=2 #Bottle12=2m(
which(hammond_data$Depth.m=="") #integer(0)
which(is.na(hammond_data$Depth.m))

#Set unique IDs based on year, month, bottle number
head(paste(hammond_data$Year, hammond_data$Month, hammond_data$CTD.Bottle.No, sep="."))
hammond_data$SampleID=paste(hammond_data$Year, hammond_data$Month, hammond_data$CTD.Bottle.No, sep=".")
length(unique(hammond_data$SampleID)) #1034 #now 1053
nrow(hammond_data)
which(duplicated(hammond_data$SampleID)) #integer(0)

#There are ~20 duplicates <- remove em
hammond_data$SampleID[which(duplicated(hammond_data$SampleID))]
#There are X bottles from February 2016 that do not have a bottle number
View(hammond_data[which(hammond_data$Year==2016 & hammond_data$Month==2),])
#Ok fill in these bottle numbers  
hammond_data$CTD.Bottle.No[which(hammond_data$Year==2016 & hammond_data$Month==2)]=c(12:5)
#Next problematic date is 11/2012
View(hammond_data[which(hammond_data$Year==2012 & hammond_data$Month==11),])
#Ok so there are two dates in this month: 28th and the 12th #Wait, looking at the original data, this is a typo
#They sampled on 11/28/2012 and 12/11/2012 #s#!t #Change the month
hammond_data$Month[which(hammond_data$Year==2012 & hammond_data$Month==11 & hammond_data$Day==12)]=12
#Reset sample IDs

#Write in a column for Year.Mo, important for grepping DCM 
hammond_data$Year.Mo=paste(hammond_data$Year, hammond_data$Month, sep=".")

#Subset hammond_data to only take columns I want: depth.m, Si, fluor., Ph, etc
dim(hammond_data) #1053 17
grep("Depth.m|Si|Temp|Ph|Fluor|SampleID|Year.Mo", colnames(hammond_data))
hammond_subs=hammond_data[,grep("Depth.m|Si|Temp|Ph|Fluor|SampleID|Year.Mo", colnames(hammond_data))]
dim(hammond_subs) #1053 8
names(hammond_subs)

#Set up a new dataframe 
sort(unique(paste(hammond_data$Year, hammond_data$Month, sep="."))) #from 10.2008-5.2018
length(sort(unique(paste(hammond_data$Year, hammond_data$Month, sep=".")))) #91 unique sampling dates
names(hammond_data) #WANT: Si, pH, Fluorescence
#Really going to linearly interpolate across just 3 data points? Ok.

int_hammond=as.data.frame(matrix(ncol=1, nrow=(92*13))) #12 bottles + 1 for DCM
dim(int_hammond)
rep(unique(paste(hammond_data$Year, hammond_data$Month, sep=".")), each=13)
c(12:9, "DCM", 8:1)
length(paste(rep(unique(paste(hammond_data$Year, hammond_data$Month, sep=".")), each=13), c(12:9, "DCM", 8:1), sep="."))
#rownames(int_hammond)=paste(rep(unique(paste(hammond_data$Year, hammond_data$Month, sep=".")), each=13), c(12:9, "DCM", 8:1), sep=".")
#colnames(int_hammond)=c("depth", "si", "pH", "fluorescence")
#Set up SampleID and Year.Mo
dim(int_hammond) #1196 1
names(int_hammond)=c("SampleID", "Year.Mo")
int_hammond$SampleID=paste(rep(unique(paste(hammond_data$Year, hammond_data$Month, sep=".")), each=13), c(12:9, "DCM", 8:1), sep=".")
int_hammond$Year.Mo=paste(rep(unique(paste(hammond_data$Year, hammond_data$Month, sep=".")), each=13))
length(unique(int_hammond$SampleID)) #1196
length(unique(int_hammond$Year.Mo)) #92 #ok

#Write in depth data with tidyverse
library(dplyr)
hammond_fulljoin=full_join(as.data.frame(int_hammond), as.data.frame(hammond_subs), by="SampleID")
dim(hammond_fulljoin) #1196 9

#Write in DCM depth 
head(hammond_fulljoin$SampleID, n=20)
hammond_fulljoin$SampleID[which(is.na(hammond_fulljoin$Depth.m))] #Why are some of these not "DCM"??
length(grep("DCM", hammond_fulljoin$SampleID))
length(which(is.na(hammond_fulljoin$Depth.m))) #143
hammond_fulljoin$Depth.m[which(is.na(hammond_fulljoin$Depth.m))]
View(hammond_fulljoin[which(is.na(hammond_fulljoin$Depth.m)),])

nrow(int_hammond)-nrow(hammond_subs) #143, that looks ok

#Checking a specific sample
hammond_fulljoin[grep("2017.9.4",hammond_fulljoin$SampleID),] #Ok, no data
grep("2017.9.4", hammond_subs$SampleID) #integer(0)

#HMMM
for(i in hammond_fulljoin$SampleID[which(is.na(hammond_fulljoin$Depth.m))][grep("DCM", hammond_fulljoin$SampleID[which(is.na(hammond_fulljoin$Depth.m))], invert=T)]){
  print(grep(i, hammond_subs$SampleID))
} #a couple of samples are wrong - they are missing from hammond_fulljoin when they shouldn't be #There are 4 missing
hammond_subs$SampleID[962:964]
hammond_fulljoin[grep("2017.9.12", hammond_fulljoin$SampleID),]

#Check to see if data are all similar
sort(hammond_fulljoin$Temperature.C[which(as.numeric(hammond_fulljoin$Depth.m)<10)])
sort(hammond_fulljoin$Si.um[which(as.numeric(hammond_fulljoin$Depth.m)<10)])

#Ok write in DCM depth with a for-loop
for(i in grep("DCM", hammond_fulljoin$SampleID)){
  #print(hammond_fulljoin$SampleID[i])
  #print(hammond_fulljoin$Depth.m[i])
  a=grep(hammond_fulljoin$SampleID[i], liv_subs$Year.Mo.Dep)
  #print(a)
  if(length(a)<1){
    print(hammond_fulljoin$SampleID[i])
    hammond_fulljoin$Depth.m[i]=NA
  }else{
    hammond_fulljoin$Depth.m[i]=liv_subs$Depth_m[a]
  }
  #hammond_fulljoin$Depth.m[i]=liv_subs$Depth_m[grep(hammond_fulljoin$SampleID[i], liv_subs$Year.Mo.Dep)]
}
#There are 3 dates with no DCM data in liv_subs: 2010.1, 2010.3, 2010.4
sort(hammond_fulljoin$Depth.m[grep("DCM", hammond_fulljoin$SampleID)]) #Ok they are showing up as characters #But ok 

#Write out hammond_fulljoin data, use LibraOffice to sort by date, then by depth #also remove second year.mo column, not sure how that got in there
setwd("/Users/cfletcherh1/Desktop/DataAnalysis/SPOT_16S_18S_eLSA/")

(file="abio_data_unformatted/temp_hammond_data_for_interpolation_04.20.2021.tsv", x=hammond_fulljoin, quote=F, row.names = F, sep="\t")
#Read back in
hammond=read.table("abio_data_unformatted/temp_hammond_data_for_interpolation_04.20.2021.tsv", header=T, stringsAsFactors = F, sep="\t")
dim(hammond) #1196 9
head(hammond, n=20)
unique(hammond$Year.Mo.x) #wait #why TF have they put a 0 at the end of the months <12?? 
hammond$SampleID[grep(".10.", hammond$SampleID)]

#Ok, re-write depth (bottle #), year, month based on SampleID
hammond$Year=2008
hammond$Month=1
hammond$CTDNo=12

for(i in c(1:nrow(hammond))){
  a=strsplit(hammond$SampleID[i], split=".", fixed=T)
  hammond$Year[i]=a[[1]][1]
  hammond$Month[i]=a[[1]][2]
  hammond$CTDNo[i]=a[[1]][3]
}

hammond$Year.Mo=paste(hammond$Year, hammond$Month, sep=".")

#Double check
for(i in sample(x=c(1:nrow(hammond)), size=12)){
  print(hammond[i, grep("Sample|Year|Month|CTD", colnames(hammond))])
}

#Subset out fake Year.Mo data 
dim(hammond) #9 13
grep("Year.Mo.x|Year.Mo.y", names(hammond))
hammond=hammond[,-grep("Year.Mo.x|Year.Mo.y", names(hammond))]
dim(hammond) #1196 11

#For each DCM depth take the weighted average of one row above and below DCM depth
#For Si, temp (funsies), Ph, Fluor
for(i in grep("DCM", hammond$SampleID)){
  print(hammond$SampleID[i])
  #print(c(hammond$Depth.m[i-1], hammond$Depth.m[i+1]))
  #print(hammond$Depth.m[i])
  #print(c(hammond$Si.um[i-1], hammond$Si.um[i+1]))
  
  #Si
  hammond$Si.um[i]=weighted.mean(x=c(as.numeric(hammond$Si.um[i-1]), as.numeric(hammond$Si.um[i+1])), w=c(as.numeric(hammond$Depth.m[i-1]), as.numeric(hammond$Depth.m[i+1])))
  #Temp
  hammond$Temperature.C=weighted.mean(x=c(as.numeric(hammond$Temperature.C[i-1]), as.numeric(hammond$Temperature.C[i+1])), w=c(as.numeric(hammond$Depth.m[i-1]), as.numeric(hammond$Depth.m[i+1])))
  #Ph
  hammond$Ph=weighted.mean(x=c(as.numeric(hammond$Ph[i-1]), as.numeric(hammond$Ph[i+1])), w=c(as.numeric(hammond$Depth.m[i-1]), as.numeric(hammond$Depth.m[i+1])))
  #Fluorescence
  hammond$Fluorescence.mg.per.cubic.m=weighted.mean(x=c(as.numeric(hammond$Fluorescence.mg.per.cubic.m[i-1]), as.numeric(hammond$Fluorescence.mg.per.cubic.m[i+1])), w=c(as.numeric(hammond$Depth.m[i-1]), as.numeric(hammond$Depth.m[i+1])))
}
#FANTASTIC! 

#Write out this temporary file

(hammond, file="abio_data_unformatted/temp_hammond_data_interpolated_04.20.2021.tsv", quote=F, row.names=F, sep="\t")

#Write into liv_subs
names(liv_subs)
head(liv_subs$Year.Mo.Dep)

#Write in a column for SPOT depths 
hammond$SPOT.Dep="NA"
#Bottle 12=5m
length(which(hammond$CTDNo==12))
hammond$SPOT.Dep[which(hammond$CTDNo==12)]="5m"
hammond$Depth.m[which(hammond$CTDNo==12)]
#Bottle DCM= ...DCM?...
hammond$Depth.m[which(hammond$CTDNo=="DCM")]
length(which(hammond$CTDNo=="DCM")) #92
hammond$SPOT.Dep[which(hammond$CTDNo=="DCM")]="DCM"
#Bottle 5=150m
hammond$Depth.m[which(hammond$CTDNo==5)]
length(which(hammond$CTDNo==5)) #92
hammond$SPOT.Dep[which(hammond$CTDNo==5)]="150m" #Ranges by as much as 50m
#Bottle3=500m
length(which(hammond$CTDNo==3)) #92
hammond$Depth.m[which(hammond$CTDNo==3)]
hammond$SPOT.Dep[which(hammond$CTDNo==3)]="500m"
#Bottle1=890m
length(which(hammond$CTDNo==1))
hammond$Depth.m[which(hammond$CTDNo==1)]
hammond$SPOT.Dep[which(hammond$CTDNo==1)]="890m"

#Make new Sample IDs based on SPOT Dep
hammond$Year.Mo.SPOTDep=paste(hammond$Year, hammond$Month, hammond$SPOT.Dep, sep=".")

#Subset Hammond data
length(which(hammond$SPOT.Dep!="NA")) #460
dim(hammond) #1196 13
subs_hammond=hammond[which(hammond$SPOT.Dep!="NA"),]
dim(subs_hammond) #460 13
which(subs_hammond$SPOT.Dep=="NA") #of course this means I removed some real dates too
length(grep("150m", subs_hammond$SPOT.Dep))

#Write fake columns into liv_subs
liv_subs$Si=as.numeric("NA")
liv_subs$pH=as.numeric("NA")
liv_subs$Fluorescence.mg.per.cubic.m=as.numeric("NA")

#WRITE INTO LIV DATA
missing_samples=c()
for(i in 1:nrow(subs_hammond)){
  a=grep(subs_hammond$Year.Mo.SPOTDep[i], liv_subs$Year.Mo.Dep)
  if(length(a)<1){
    missing_samples=c(missing_samples,subs_hammond$Year.Mo.SPOTDep[i])
  }else{
    liv_subs$Si[a]=subs_hammond$Si.um[i]
    liv_subs$pH[a]=subs_hammond$Ph[i]
    liv_subs$Fluorescence.mg.per.cubic.m[a]=subs_hammond$Fluorescence.mg.per.cubic.m[i]
  }
}
missing_samples #AWESOME!! Yeah, that looks right

#Re-order columns in liv_subs, and also get rid of those we don't need
names(liv_subs)
names(liv_subs)[c(1:5, grep("SampleID", colnames(liv_subs)),6,7,9,8,11:13,17:19,14,16)]
dim(liv_subs) #860 19
liv_FINAL=liv_subs[,c(1:5, grep("SampleID", colnames(liv_subs)),6,7,9,8,11:13,17:19,14,16)]
dim(liv_FINAL) #860 18
names(liv_FINAL)

#WRITE IT OUT!! 

(liv_FINAL, "ModifiedFiles/SPOT_nutrient_CTD_prod_Hammond_data_04.20.2021.tsv", quote=F, row.names=F, sep="\t")
#FUCK YEAH!!!! 

#PART II: Physical oceanography data (upwelling indices, MEI)-----

#5. CUTI and BEUTI data-----
#a. CUTI data
CUTI_data=read.csv("abio_data_unformatted/NOAA_CUTI_monthly.csv", header=T, stringsAsFactors = F)
dim(CUTI_data)
names(CUTI_data) #Want 33N and want data from SPOT years only 
#Create a column for year_month to match that of Liv's data
CUTI_data$Year_Month=paste(CUTI_data$year, CUTI_data$month, sep="_")
head(CUTI_data$Year_Month)


#b. BEUTI data
BEUTI_data=read.csv("abio_data_unformatted/NOAA_BEUTI_monthly.csv", header=T, stringsAsFactors = F)
dim(BEUTI_data)
dim(CUTI_data)
BEUTI_data$Year_Month=paste(BEUTI_data$year, BEUTI_data$month, sep="_")
names(CUTI_data)==names(BEUTI_data) #All TRUE

#6. MEI data-----
MEI_data=read.csv("abio_data_unformatted/NOAA_MEI_Formatted.csv", header=T, stringsAsFactors = F)
dim(MEI_data) #42 13
names(MEI_data)
MEI_data$YEAR

#Gonna have to make this into a list with year_month, MEI ix
head(list(MEI_data))
MEI_list_1=as.list(as.numeric(t(as.matrix(MEI_data[,2:13])))) #Col 1= YEAR, remove that
head(MEI_list_1, n=24) #It goes down by row, not column 
MEI_data[1:2,]
length(MEI_list_1)
MEI_list_2=c()
for(i in c(1:length(MEI_list_1))){
  MEI_list_2=c(MEI_list_2, MEI_list_1[[i]][1])
}
length(MEI_list_2) #504
head(MEI_list_2, n=24)
MEI_data[1:2,]

#Now get a vector of dates up (month-year)
rep(MEI_data$YEAR[1]:MEI_data$YEAR[nrow(MEI_data)], each=12)
rep(1:12, times=length(MEI_data$YEAR))
MEI_yr_mo=paste(rep(MEI_data$YEAR[1]:MEI_data$YEAR[nrow(MEI_data)], each=12), rep(1:12, times=length(MEI_data$YEAR)), sep="_")
length(MEI_yr_mo) #504
length(unique(MEI_yr_mo)) #504
length(unique(MEI_yr_mo))==length(MEI_list_2) #TRUE

#Write into a new matrix/df
MEI_df=as.data.frame(matrix(ncol=2, nrow=504))
colnames(MEI_df)=c("Year_Month", "MEI_ix")
MEI_df$Year_Month=MEI_yr_mo
MEI_df$MEI_ix=MEI_list_2

#Check it over
head(MEI_data)
head(MEI_df, n=24)
MEI_df[sample(x=c(1:nrow(MEI_df)), size=10, replace=F),]
View(MEI_data)

#7. Bakun upwelling data-----
bakun_data=read.csv("abio_data_unformatted/NOAA_BakunUpwellingData_200_2018_33N_241E.csv", stringsAsFactors = F, header = T)
dim(bakun_data) #885 12
head(bakun_data)
#Get rid of the first row, which only has the units for each measurement, not actual data
#But store it
meas=bakun_data[1,]
bakun_data=bakun_data[-1,]
dim(bakun_data) #884 12

#Pull out only the lat and longs I am interested in
unique(bakun_data$latitude) #Want 33.5N
unique(bakun_data$longitude) #Want 241.5E

length(which(bakun_data$latitude==33.5 & bakun_data$longitude==241.5)) #221
bakun_subs=bakun_data[which(bakun_data$latitude==33.5 & bakun_data$longitude==241.5),]
dim(bakun_subs) #221 12

#Re-format "time" to be month-year in the same format as CUTI and BEUIT data
#So no leading zeros
mo_yr <- c()
for(i in c(1:nrow(bakun_subs))){
  a <- strsplit(bakun_subs$time[i], split="-", fixed=T)
  if(as.numeric(a[[1]][2])<10){
    b <- strsplit(a[[1]][2], split="0", fixed=T)
    mo_yr=c(mo_yr, paste(a[[1]][1], b[[1]][2], sep="_"))
  } else{
    #print("duck")  
    mo_yr <- c(mo_yr, paste(a[[1]][1], a[[1]][2], sep="_"))
  }
}
length(mo_yr) #221
length(unique(mo_yr)) #221, good
bakun_subs$Year_Month=mo_yr

#8. MODIS Satellite chlorophyll data-----
chl_data=read.table("OriginalFiles/TS_TMWchlaSmdayAverages_x-118.4_y33.55_t20020716120000_T20180716120000.asc.txt", header=T, stringsAsFactors = F, sep="\t")
dim(chl_data) #193 6
head(chl_data) #Time is reported as seconds since 1970-01-01T00Z

#Format time from seconds to date.time
as.POSIXct(chl_data$TIME, origin='1970-01-01 00:00:00') #THAT WORKS! 

#Write a fake column of data 
chl_data$Date="1970-01-01"

#Write this into a for-loop 
for(i in c(1:nrow(chl_data))){
  a=as.POSIXct(chl_data$TIME[i], origin='1970-01-01 00:00:00')
  b=strsplit(as.character(a), fixed=T, split=" ")
  chl_data$Date[i]=b[[1]][1]
}
head(chl_data$Date)
tail(chl_data$Date)

#Split again to create year and month 
chl_data$Year_Month="1970_01"

for(i in c(1:nrow(chl_data))){
  #print(chl_data$Date[i])
  a=strsplit(chl_data$Date[i], fixed=T, split="-")
  YEAR=a[[1]][1]
  b=a[[1]][2]
  #print(b)
  if(as.numeric(b)<10){
    d=strsplit(b, split="", fixed=T)
    MONTH=d[[1]][2]
  }else{
    MONTH=b
  }
  chl_data$Year_Month[i]=paste(YEAR, MONTH, sep="_")
}
unique(chl_data$Year_Month)
 
#9. MODIS SST data-----
sst_data=read.table("OriginalFiles/TS_TMWsstdSmdayAverages_x-118.4_y33.55_t20020716120000_T20180716120000.asc.txt", header=T, stringsAsFactors = F, sep="\t")
head(sst_data)
dim(sst_data) #191 6
names(sst_data)

#Write a fake column for date
sst_data$Date="1970-01-01"
dim(sst_data) #191 7

#Overwrite with actual data
for(i in c(1:nrow(sst_data))){
  a=as.POSIXct(sst_data$TIME[i], origin='1970-01-01 00:00:00')
  b=strsplit(as.character(a), split=" ")
  sst_data$Date[i]=b[[1]][1]
}
head(sst_data$Date)
tail(sst_data$Date)

#Split again to get month_date
sst_data$Year_Month="1970_01"

#With no leading zeros
for(i in c(1:nrow(sst_data))){
  a=sst_data$Date[i]
  b=strsplit(a, split="-")
  year=b[[1]][1]
  if(b[[1]][2]<10){
    d=strsplit(b[[1]][2], split="")
    month=d[[1]][2]
    #print(month)
  } else{
    month=b[[1]][2]
  }
  e=paste(year, month, sep="_")
  #print(e)
  sst_data$Year_Month[i]=e
}

#Check if chl data and sst data have same dates (different dimensions)
nrow(chl_data) #193
nrow(sst_data) #191
chl_data$Date==sst_data$Date #Ok until row 110

#Ok try this: 
missing_dates=c()
for(i in chl_data$Year_Month){
  a=(grep(i, sst_data$Year_Month)) #from longer to shorter
  print(a)
  if(length(a)<1)
    missing_dates=c(missing_dates, i)
}
length(missing_dates) #2
missing_dates #2011_9 2013_12 #ok #idk why those two are missing, but ok

#Set up a new data frame to write sst data into
sst_data_2=as.data.frame(matrix(nrow=nrow(chl_data), ncol=2))
dim(sst_data_2) #193 2
names(sst_data_2)=c("Year_Month", "SST")

sst_data_2$Year_Month=paste(c(rep(2002, times=length(c(7:12))), rep(2003:2017, each=12), rep(2018, times=length(c(1:7)))), c(7:12, rep(c(1:12), times=15), c(1:7)), sep="_")

library("dplyr")
SST_data=full_join(sst_data_2, sst_data, by="Year_Month")
dim(SST_data) #193 10 #Ok joined all the columns 
names(SST_data)
which(is.na(SST_data$TMWsstd)) #4?
SST_data[which(is.na(SST_data$TMWsstd)),] #ok, 4 #3 of 4 are from Fall 2011 (sept, oct, dec.)

#10. Upwelling-derived NO2 data----
no2_data=read.table("OriginalFiles/SPOT_5m_Dura_SST_Nitrate.csv", header=T, strings=F, sep=",")
dim(no2_data) #117 3
head(no2_data) #starts 2008 05
tail(no2_data) #ends 2018 07

no2_data=as.data.frame(as.matrix(no2_data))
typeof(no2_data)

no2_data$Year_Month="2008_5"

for(i in c(1:nrow(no2_data))){
  a=no2_data$Date[i]
  b=strsplit(a, split="-")
  year=b[[1]][1]
  if(b[[1]][2]<10){
    d=strsplit(b[[1]][2], split="")
    month=d[[1]][2]
  }else{
    month=b[[1]][2]
  }
  no2_data$Year_Month[i]=(paste(year, month, sep="_"))
}
head(no2_data$Year_Month)
tail(no2_data$Year_Month)

#Subset out to remove duplicates 
no2_data$Year_Month[which(duplicated(no2_data$Year_Month))] #Get rid of these
no2_data$Year_Month[which(duplicated(no2_data$Year_Month))-1] #prior rows

no2_data=no2_data[-which(duplicated(no2_data$Year_Month)),]
dim(no2_data) #113 4
names(no2_data)

#Write all the data from #5-on into a single matrix-----
PO_data=as.data.frame(matrix(ncol=length(c("Year", "Month", "Year_Month", "CUTI_ix", "BEUTI_ix", "MEI_ix", colnames(bakun_subs)[4:12], "Chl", "SST", "model_no2", "model_sst")), nrow=length(c(rep("2000", times=5), rep(c(2001:2020), each=12)))))
colnames(PO_data)=c("Year", "Month", "Year_Month", "CUTI_ix", "BEUTI_ix", "MEI_ix", colnames(bakun_subs)[4:12], "MODIS_Chl", "MODIS_SST", "Model_NO2", "Model_SST")
dim(PO_data) #245 18

#Write in dates data
PO_data$Year=c(rep("2000", times=5), rep(c(2001:2020), each=12))
PO_data$Month=c(8:12, rep(c(1:12), times=20))
PO_data$Year_Month=paste(PO_data$Year, PO_data$Month, sep="_")

#Add in CUTI data 
grep(PO_data$Year_Month[1], CUTI_data$Year_Month)
grep(PO_data$Year_Month[nrow(PO_data)], CUTI_data$Year_Month)
length(CUTI_data$X33N[grep(PO_data$Year_Month[1], CUTI_data$Year_Month):grep(PO_data$Year_Month[nrow(PO_data)], CUTI_data$Year_Month)])
PO_data$CUTI_ix=CUTI_data$X33N[grep(PO_data$Year_Month[1], CUTI_data$Year_Month):grep(PO_data$Year_Month[nrow(PO_data)], CUTI_data$Year_Month)]

#Add in BEUTI data
length(BEUTI_data$X33N[grep(PO_data$Year_Month[1], BEUTI_data$Year_Month):grep(PO_data$Year_Month[nrow(PO_data)], BEUTI_data$Year_Month)])
PO_data$BEUTI_ix=BEUTI_data$X33N[grep(PO_data$Year_Month[1], BEUTI_data$Year_Month):grep(PO_data$Year_Month[nrow(PO_data)], BEUTI_data$Year_Month)]

#Add in MEI data
length(MEI_df$MEI_ix[grep(PO_data$Year_Month[1], MEI_df$Year_Month):grep(PO_data$Year_Month[nrow(PO_data)], MEI_df$Year_Month)])
PO_data$MEI_ix=MEI_df$MEI_ix[grep(PO_data$Year_Month[1], MEI_df$Year_Month):grep(PO_data$Year_Month[nrow(PO_data)], MEI_df$Year_Month)]

#Add in Bakun data (need to add in X's as well)
length(grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs))
length(c(grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs), rep("na", times=31)))
names(bakun_subs)
names(PO_data)

PO_data$pmsl=c(bakun_subs$pmsl[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$u_mean=c(bakun_subs$u_mean[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$v_mean=c(bakun_subs$v_mean[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$uv_mag_mean=c(bakun_subs$uv_mag_mean[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$taux_mean=c(bakun_subs$taux_mean[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$tauy_mean=c(bakun_subs$tauy_mean[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$curl=c(bakun_subs$curl[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$ektrx=c(bakun_subs$ektrx[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))
PO_data$ektry=c(bakun_subs$ektry[grep(PO_data$Year_Month[1], bakun_subs$Year_Month):nrow(bakun_subs)], rep("na", times=31))

#Add in MODIS Chl data
length(chl_data$TMWchla[grep(chl_data$Year_Month[1], PO_data$Year_Month):grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month)])
length(c(rep("na", times=grep(chl_data$Year_Month[1], PO_data$Year_Month)-1), chl_data$TMWchla[grep(chl_data$Year_Month[1], PO_data$Year_Month):grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month)], rep("na", times=length(grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month):nrow(PO_data))-1))) #245

length(c(rep("na", times=grep(chl_data$Year_Month[1], PO_data$Year_Month)-1), chl_data$TMWchla, rep("na", times=length(grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month):nrow(PO_data))-1))) #245

PO_data$MODIS_Chl=c(rep("na", times=grep(chl_data$Year_Month[1], PO_data$Year_Month)-1), chl_data$TMWchla, rep("na", times=length(grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month):nrow(PO_data))-1))

#Double check that
PO_data$MODIS_Chl[grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month)] #Should be 0.7097582 #ok
tail(chl_data)
PO_data$MODIS_Chl[grep(chl_data$Year_Month[nrow(chl_data)], PO_data$Year_Month)-1] #ok
head(chl_data)
PO_data$MODIS_Chl[grep(chl_data$Year_Month[1], PO_data$Year_Month)] #ok
grep(PO_data$MODIS_Chl[grep(chl_data$Year_Month[1], PO_data$Year_Month)], chl_data$TMWchla) #ok

#Add in MODIS SST data
grep(SST_data$Year_Month[1], PO_data$Year_Month) #24
grep(SST_data$Year_Month[nrow(SST_data)], PO_data$Year_Month) #216
length(c(rep("na", times=grep(SST_data$Year_Month[1], PO_data$Year_Month)-1), SST_data$TMWsstd, rep("na", times=nrow(PO_data)-grep(SST_data$Year_Month[nrow(SST_data)], PO_data$Year_Month)))) #245

PO_data$MODIS_SST=c(rep("na", times=grep(SST_data$Year_Month[1], PO_data$Year_Month)-1), SST_data$TMWsstd, rep("na", times=nrow(PO_data)-grep(SST_data$Year_Month[nrow(SST_data)], PO_data$Year_Month)))
names(PO_data)

PO_data$MODIS_SST[which(is.na(PO_data$MODIS_SST))]="na"
which(is.na(PO_data$MODIS_SST))

#Add in model data (SST and NO2)
names(PO_data)
temp=full_join(PO_data, no2_data, by="Year_Month")
names(temp) #WHaT #it worked?? 
dim(temp) #245 22

names(PO_data)
dim(PO_data) #245 19
PO_data$Model_SST=temp$SST_Daily_Average
PO_data$Model_NO2=temp$Nitrate_Snyder_Daily_Average
names(PO_data)

length(which(is.na(PO_data$Model_SST))) #132
PO_data$Model_SST[which(is.na(PO_data$Model_SST))]="na"

PO_data$Model_NO2[which(is.na(PO_data$Model_NO2))]="na"
length(which(is.na(PO_data$Model_NO2)))

#Write it out!
dim(PO_data) #245 19
tail(PO_data)


write.table(x=PO_data, file="ModifiedFiles/Upwelling_Chl_NO2_Data_CA_Coast_2000_2018.tsv", row.names=F, sep="\t", quote=F)
