#! /bin/bash
# This script downloads Chandra observation data based on a list of observational ID.

# Define the input file containing the list of observations to download.
input="/Users/juliaahlvind/Documents/projekt_1/Chandra/main_list_download.txt"

# Change directory to the location where the Chandra data should be stored.
cd ../../Chandra/data

while read Name obsid Type epoch dist expt dateSN dateXray ra dec
do
    # Call the function to download the Chandra data for the given obsid.
    download_chandra_obsid $obsid

# Redirect the input for the while loop from the specified input file.
done < "$input"

cd /Users/juliaahlvind/Documents/projekt_1/pipelines/chandra
