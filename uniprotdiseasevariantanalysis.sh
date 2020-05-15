#!/bin/bash

#Author: Vijay Nagarajan PhD

#This script gets the disease name from user, gets all the variants associated with the users disease, identifies the unique genes associated with the disease and also save the counts of variants for each of the genes associated with the disease

#get the user input for their disease of interest
echo "Enter your disease of interest: "

#read what the user types and assign it to the variable 'mydisease'
read mydisease

#use $mydisease to get the corresponding variant data from ebi, store it as usersdiseasevariants.txt
curl -X GET --header 'Accept:text/x-gff' 'https://www.ebi.ac.uk/proteins/api/variation?offset=0&size=100&disease='"$mydisease" > "$mydisease"variants.txt

#grep all lines other than comments
#cut the first column
#identify uniq genes and save it to usersdiseasegenes.txt file
grep -v "##" "$mydisease"variants.txt | cut -f1 | uniq > "$mydisease"genes.txt

#count the number of genes in the usersdiseasegenes.txt file, display only the line number
echo "There are" $( wc -l "$mydisease"genes.txt | cut -f1 -d' ' )  "genes associated with $mydisease"

echo "All the unique genes associated with $mydisease are saved in the $mydisease""genes.txt file"

#get all non comment lines, cut first field, sort, output counts of unique genes in to usersdiseasegenesvariantcounts.txt
grep -v "##" "$mydisease"variants.txt | cut -f1 | sort | uniq -c > "$mydisease"genesvariantcounts.txt







