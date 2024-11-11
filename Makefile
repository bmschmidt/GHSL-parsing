data/*.tif:
  wget https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E${year}_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E${year}_GLOBE_R2023A_54009_100_V1_0_R${row}_C${col}.zip`;

prep-quadfeather:
	uv init
	uv add git+https://github.com/bmschmidt/quadfeather

quadtree:
	uv run quadtree.py ~/scrolly_tiles/everyone

points/%.feather: points/%.parquet
	uv run shuffle.py $<

DOWNLOAD_DIR := data
YEARS := 1975 1980 1985 1990 1995 2000 2005 2010 2015 2020 2025
TUPLES_FILE := tuples.txt

TUPLES_FILE := tuples.txt

# Create the download directory if it doesn't exist
$(DOWNLOAD_DIR):
	mkdir -p $(DOWNLOAD_DIR)

# Main target to download all files
download: $(DOWNLOAD_DIR) $(YEARS:%=$(DOWNLOAD_DIR)/%_downloaded)

# Target to download files for each year
$(DOWNLOAD_DIR)/%_downloaded: | $(DOWNLOAD_DIR)
	@year=$* ; \
	while IFS=" " read -r row col; do \
	  tif_file=$(DOWNLOAD_DIR)/GHS_POP_E$${year}_GLOBE_R2023A_54009_100_V1_0_R$${row}_C$${col}.tif ; \
	  if [ ! -f $$tif_file ]; then \
	    url=https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E$${year}_GLOBE_R2023A_54009_100/V1-0/tiles/GHS_POP_E$${year}_GLOBE_R2023A_54009_100_V1_0_R$${row}_C$${col}.zip ; \
	    zip_file=$(DOWNLOAD_DIR)/GHS_POP_E$${year}_GLOBE_R2023A_54009_100_V1_0_R$${row}_C$${col}.zip ; \
	    echo "Downloading $$url" ; \
	    wget -P $(DOWNLOAD_DIR) $$url || echo "Failed to download $$url" ; \
	    if [ -f $$zip_file ]; then \
	      echo "Extracting $$zip_file" ; \
	      unzip -o $$zip_file -d $(DOWNLOAD_DIR) && rm -f $$zip_file ; \
	    fi ; \
	  else \
	    echo "Skipping download, $$tif_file already exists" ; \
	  fi ; \
	done < $(TUPLES_FILE)
	@touch $@

.PHONY: download clean

# Clean up downloaded files
clean:
	rm -rf $(DOWNLOAD_DIR)