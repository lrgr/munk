# Directories
RAW_DIR=raw
SC_SP_DIR=sc-sp

# Input files
HOMOLOGENE_RAW_URL=http://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.xml.gz
HOMOLOGENE_RAW=$(RAW_DIR)/homologene.xml.gz

# Output files
SC_SP_HOMOLOGS=$(SC_SP_DIR)/sc-sp-homologs.txt
SP_SC_HOMOLOGS=$(SC_SP_DIR)/sp-sc-homologs.txt

# Name mapper
NAME_MAPPER_DIR=../name_mapping
HOMOLOG_MAPPER=$(NAME_MAPPER_DIR)/standardize_homolog_names.py
IDMAPPING=$(NAME_MAPPER_DIR)/mappings/sc-sp-idmapping.json

# tempfiles
TMPFILE = $(RAW_DIR)/tmpfile_homologs

# Standard rules
all: download process

setup:
	mkdir -p $(RAW_DIR) $(SC_SP_DIR)

clean:
	rm -r $(RAW_DIR) $(SC_SP_DIR)

# Download the raw homolog data
download: $(HOMOLOGENE_RAW) $(IDMAPPING)

$(IDMAPPING): | setup
	wget $(IDMAPPING_URL) -O $(IDMAPPING)
$(HOMOLOGENE_RAW): | setup
	wget $(HOMOLOGENE_RAW_URL) -O $(HOMOLOGENE_RAW)

# Process the homolog data
process: $(SC_SP_HOMOLOGS) $(SP_SC_HOMOLOGS)
$(SC_SP_HOMOLOGS): $(HOMOLOGENE_RAW)
	/bin/rm -f $(TMPFILE) 
	python extract_homologs.py \
		-if $(HOMOLOGENE_RAW) \
		-of $(TMPFILE)
	python $(HOMOLOG_MAPPER) \
		-if $(TMPFILE) \
		-of $@ \
		-uf $(IDMAPPING)

$(SP_SC_HOMOLOGS): $(SC_SP_HOMOLOGS)
	sed 's/ /	/g' $(SC_SP_HOMOLOGS) > sc-sp/tmp.txt
	python swap_cols.py sc-sp/tmp.txt  $@
	rm sc-sp/tmp.txt
