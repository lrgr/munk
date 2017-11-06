# Settings
ifndef SPECIES
	SPECIES=sp sc
endif

# Directories
RAW_DIR=raw

# Look up tables
files.sp := $(RAW_DIR)/BIOGRID-ORGANISM-Schizosaccharomyces_pombe_972h-3.4.149.mitab.txt
files.sc := $(RAW_DIR)/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.149.mitab.txt

# Files
BIOGRID_RAW=$(RAW_DIR)/BIOGRID-ORGANISM-3.4.149.mitab.zip
NETWORK_FILES=$(foreach S,$(SPECIES),$(S)/biogrid-$(S)-network.txt)
STD_NETWORK_FILES=$(foreach S,$(SPECIES),$(S)/biogrid-$(S)-std-network.txt)

# URLs
BIOGRID_URL=https://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.149/BIOGRID-ORGANISM-3.4.149.mitab.zip

# Name mapper
NAME_MAPPER_DIR=../../name_mapping
NETWORK_MAPPER=$(NAME_MAPPER_DIR)/standardize_network_names.py
IDMAPPING=$(NAME_MAPPER_DIR)/mappings/sc-sp-idmapping.json

# General commands
all: download process standardize

setup:
	mkdir -p $(RAW_DIR) $(SPECIES)

clean:
	rm -r $(RAW_DIR) $(SPECIES)

# Download data
download: $(BIOGRID_RAW)

$(BIOGRID_RAW): | setup
	wget -q $(BIOGRID_URL) -O $(BIOGRID_RAW);
	cd $(RAW_DIR) && unzip $(lastword $(subst /, ,$@))

# Process the data files
process: $(NETWORK_FILES) 

$(NETWORK_FILES): $(BIOGRID_RAW)
	INPUT_FILE=${files.${word 2,$(subst -, ,$(lastword $(subst /, ,$@)))}} && \
	python process_biogrid.py \
		-if $$INPUT_FILE \
		-of $@

standardize: $(STD_NETWORK_FILES)

$(STD_NETWORK_FILES): $(NETWORK_FILES)
	python $(NETWORK_MAPPER) \
		-sp ${word 2,$(subst -, ,$(lastword $(subst /, ,$@)))} \
		-if $(subst -std,,$@) \
		-of $@ \
		-uf $(IDMAPPING)

