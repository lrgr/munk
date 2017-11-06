# Directories
SC_DIR=sc
RAW_DIR=raw

# Input files
RAW_COLLINS_DATA_URL=https://www.nature.com/nature/journal/v446/n7137/extref/nature05649-s2.zip
RAW_COLLINS_DATA=$(RAW_DIR)/nature05649-s2.zip
RAW_COLLINS_EMAP=$(RAW_DIR)/Collins2007-EMAP-data.tsv

# Output files
COLLINS_GIS=$(SC_DIR)/collins-sc-emap-gis
COLLINS_GIS_TSV=$(COLLINS_GIS).tsv

STD_COLLINS_GIS_PREFIX = $(SC_DIR)/collins-sc-emap-gis-std
STD_COLLINS_GIS_TSV = $(STD_COLLINS_GIS_PREFIX).tsv

# Name mapper
NAME_MAPPER_DIR=../../name_mapping
GI_MAPPER=$(NAME_MAPPER_DIR)/standardize_gi_names.py
IDMAPPING=$(NAME_MAPPER_DIR)/mappings/sc-sp-idmapping.json

# Rules
all: $(STD_COLLINS_GIS_TSV)

setup:
	mkdir -p $(SC_DIR) $(RAW_DIR)

clean:
	rm -r $(SC_DIR) $(RAW_DIR)

# Parse the raw Collins SLs
$(RAW_COLLINS_DATA): | setup
	wget -O $(RAW_COLLINS_DATA) $(RAW_COLLINS_DATA_URL)

$(RAW_COLLINS_EMAP): $(RAW_COLLINS_DATA)
	unzip -p $(RAW_COLLINS_DATA) "Supp\ data\ -\ strains\ and\ interaction\ data/Chromosome\ biology\ EMAP\ data.txt" > $(RAW_COLLINS_EMAP)

$(COLLINS_GIS_TSV): $(RAW_COLLINS_EMAP)
	python parse_collins.py -i $(RAW_COLLINS_EMAP) -o $(COLLINS_GIS)

#standardize: $(STD_COLLINS_GIS_TSV)

$(STD_COLLINS_GIS_TSV): $(COLLINS_GIS_TSV)
	python $(GI_MAPPER) \
		-sp sc \
		-if $(COLLINS_GIS_TSV) \
		-of $(STD_COLLINS_GIS_PREFIX) \
		-uf $(IDMAPPING)
