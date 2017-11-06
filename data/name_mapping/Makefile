# Directories
DB_DIR=db
MAPPINGS=mappings

# Input files
SC_IDMAPPING_URL=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/YEAST_559292_idmapping.dat.gz
SP_IDMAPPING_URL=ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/SCHPO_284812_idmapping.dat.gz
SP_IDMAPPING=$(DB_DIR)/SCHPO_284812_idmapping.dat.gz
SC_IDMAPPING=$(DB_DIR)/YEAST_559292_idmapping.dat.gz

# Output files
IDMAPPING=$(MAPPINGS)/sc-sp-idmapping.json

# Standard rules
all: download process

setup:
	mkdir -p $(DB_DIR) $(MAPPINGS)

clean:
	rm -r $(DB_DIR) $(MAPPINGS)

# Download the raw homolog data
download: $(SC_IDMAPPING) $(SP_IDMAPPING)

$(SC_IDMAPPING): | setup
	wget $(SC_IDMAPPING_URL) -O $(SC_IDMAPPING)

$(SP_IDMAPPING): | setup
	wget $(SP_IDMAPPING_URL) -O $(SP_IDMAPPING)

# Process the map data
process: $(IDMAPPING)
$(IDMAPPING): $(SC_IDMAPPING) $(SP_IDMAPPING)
	python unify_id_mappings.py \
		-i $(SP_IDMAPPING) $(SC_IDMAPPING) \
		-o $(IDMAPPING)
