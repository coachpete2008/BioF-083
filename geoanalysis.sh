#!/bin/bash

#Author: Vijay Nagarajan PhD

#This script gets the disease name from user, gets all the variants associated with the users disease, identifies the unique genes associated with the disease and also save the counts of variants for each of the genes associated with the disease

#get the user input for their GEO Series of interest
echo "Enter your GEO Series ID (Example - GSE12494) : "

#read what the user types, and assign it to the variable 'geoseriesid'
read geoseriesid

#use $geoseriesid to get the corresponding series data from geo, stored as $geoseriesid.txt
wget -O "$geoseriesid".txt 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='"$geoseriesid"'&targ=self&view=full&form=text'

#extract sampleids, using grep and cut combination
#grep "!Series_sample_id" "$geoseriesid".txt | cut -d$' ' -f3

#create a directory on the series name to store all sample dataa
mkdir $geoseriesid

#iterate through the list of sampleids and download them one by one using wget, in a for loop
for sample in $( grep "!Series_sample_id" "$geoseriesid".txt | cut -d$' ' -f3 )
	do
		#remove any trailing spaces
		sampleid=$( echo $sample | tr -d [:space:] )
		#store output filepath as a variable
		outputname=$geoseriesid"/"$sampleid".txt"
		#download and store sample date within series folder, make sure to remove white space at the end of the sampleid, using tr
		wget -O $outputname 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='$( echo $sample | tr -d [:space:] )'&targ=self&view=full&form=text'
	done






