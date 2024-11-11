#!/bin/bash

# Directory to store downloaded files
DOWNLOAD_DIR="data"

# List of years to download
years=(1975 1980 1985 1990 1995 2000 2005 2010 2015 2020 2025)

# Load specific list of row-col tuples from a file
TUPLES_FILE="tuple.txt"

# Create the download directory if it doesn't exist
mkdir -p "$DOWNLOAD_DIR"

# Base URL template
BASE_URL="https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E{year}_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E{year}_GLOBE_R2023A_54009_100_V1_0_R{row}_C{col}.zip"

# Loop through each year, row, and col to download the files
for year in "${years[@]}"; do
  while IFS=" " read -r row col; do
    # Construct the URL
    url=${BASE_URL//\{year\}/$year}
    url=${url//\{row\}/$row}
    url=${url//\{col\}/$col}
    
    # Download the file
    echo "Downloading $url"
    wget -P "$DOWNLOAD_DIR" "$url" || echo "Failed to download $url"
  done < "$TUPLES_FILE"
done
