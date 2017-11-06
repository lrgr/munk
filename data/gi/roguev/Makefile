# Settings

# Directories
RAW_DIR=raw
SP_DIR=sp
ROGUEV_POMBE_EMAP_DIR=$(RAW_DIR)/raw_data/

# Input files
ROGUEV_EMAP_ZIP_URL=http://science.sciencemag.org/highwire/filestream/589554/field_highwire_adjunct_files/0/1162609sDataset_S1.zip
ROGUEV_EMAP_ZIP=$(RAW_DIR)/1162609sDataset_S1.zip
RAW_ROGUEV_GIS=$(RAW_DIR)/Dataset_S1.txt

# Processed files
ROGUEV_GIS=$(SP_DIR)/roguev-sp-emap-gis
ROGUEV_GIS_TSV=$(ROGUEV_GIS).tsv
ROGUEV_GIS_NPY=$(ROGUEV_GIS).npy
STD_ROGUEV_TSV=$(SP_DIR)/roguev-sp-emap-gis-std.tsv
STD_ROGUEV_GIS_PREFIX=$(SP_DIR)/roguev-sp-emap-gis-std

# Name mapper
NAME_MAPPER_DIR=../../name_mapping
GI_MAPPER=$(NAME_MAPPER_DIR)/standardize_gi_names.py
IDMAPPING=$(NAME_MAPPER_DIR)/mappings/sc-sp-idmapping.json

# General commands
all: $(STD_ROGUEV_TSV)

setup:
	mkdir -p $(RAW_DIR) $(SP_DIR)
clean:
	rm -r $(RAW_DIR) $(SP_DIR)

# Download and unzip data
download: $(ROGUEV_EMAP_ZIP)

$(ROGUEV_EMAP_ZIP): | setup
	wget $(ROGUEV_EMAP_ZIP_URL) -O $(ROGUEV_EMAP_ZIP)

$(RAW_ROGUEV_GIS): $(ROGUEV_EMAP_ZIP)
	unzip -p $(ROGUEV_EMAP_ZIP) raw_data/Dataset_S1.txt > $(RAW_ROGUEV_GIS)

# Process
process: $(ROGUEV_GIS_TSV)

$(ROGUEV_GIS_TSV): $(RAW_ROGUEV_GIS)
	python parse_roguev.py \
		-i $(RAW_ROGUEV_GIS) \
		-o $(ROGUEV_GIS)

# Standardize namees
standardize: $(STD_ROGUEV_TSV)

$(STD_ROGUEV_TSV): $(ROGUEV_GIS_TSV)
	python $(GI_MAPPER) \
		-sp sp \
		-if $(ROGUEV_GIS_TSV) \
		-of $(STD_ROGUEV_GIS_PREFIX) \
		-uf $(IDMAPPING)
