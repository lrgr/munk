# Setup
OUTPUT_DIR = output/feats
SCORES_DIR = output
LOG = log.txt

$(OUTPUT_DIR):
	mkdir -p $@

# Scripts
SRC_DIR = $(realpath ../../src)
FACTOR = $(SRC_DIR)/factorized_laplacian.py
EMBED = $(SRC_DIR)/handl_embed.py

# Input Data
DATA_DIR = $(realpath ../../data)

# Networks
SP_PPI_EDGELIST = $(DATA_DIR)/ppi/biogrid/sp/biogrid-sp-std-network.txt
SC_PPI_EDGELIST = $(DATA_DIR)/ppi/biogrid/sc/biogrid-sc-std-network.txt

################################################################################
#	HANDL
################################################################################

SP_C = $(OUTPUT_DIR)/sp-c.pkl
SP_D = $(OUTPUT_DIR)/sp-d.pkl
$(SP_C) $(SP_D): $(SP_PPI_EDGELIST) | $(OUTPUT_DIR)
	python $(FACTOR) \
		-e $(SP_PPI_EDGELIST)\
		-o $(SP_C) \
		-df $(SP_D) \
		-p 1.0 \
		-l 0.05
sp-c-d: $(SP_C) $(SP_D)

SC_D = $(OUTPUT_DIR)/sc-d.pkl
SC_C = $(OUTPUT_DIR)/sc-c.pkl
$(SC_D) $(SC_C): $(SC_PPI_EDGELIST) | $(OUTPUT_DIR)
	python $(FACTOR) \
		-e $(SC_PPI_EDGELIST) \
		-o $(SC_C)\
		-df $(SC_D)\
		-l 0.05 \
		-p 1.0

sc-c-d: $(SC_D) $(SC_C)

SP_C_HAT = $(SCORES_DIR)/sc-to-sp-scores-and-matrices.pkl
CP_HOMOLOGS = $(DATA_DIR)/homologs/sc-sp/sc-sp-homologs.txt

$(SP_C_HAT): $(SP_D) $(SC_C) $(CP_HOMOLOGS) | $(OUTPUT_DIR)
	python $(EMBED)\
		-sv $(SC_C)\
		-df $(SP_D)\
		-hf $(CP_HOMOLOGS)\
		-o $(SP_C_HAT)

sp-chat: $(SP_C_HAT)

SC_C_HAT = $(SCORES_DIR)/sp-to-sc-scores-and-matrices.pkl
PC_HOMOLOGS = $(DATA_DIR)/homologs/sc-sp/sp-sc-homologs.txt

$(SC_C_HAT): $(SC_D) $(SP_C) $(PC_HOMOLOGS) | $(OUTPUT_DIR)
	python $(EMBED)\
		-sv $(SP_C)\
		-df $(SC_D)\
		-hf $(PC_HOMOLOGS)\
		-o $(SC_C_HAT)

sc-chat: $(SC_C_HAT)

all: sc-chat sp-chat
clean:
	rm -r output
