#!/usr/bin/env python

# Genetic interactions
SL = 'SL'
NON_SL = 'Non-SL'
INCONCLUSIVE = 'Inconclusive'

# Feature functions
ADD_FEATURES = 'add'
MEAN_FEATURES = 'mean'
FEATURE_FUNCTIONS = { ADD_FEATURES, MEAN_FEATURES }

# Cross-validation choices
GENES = 'genes'
GENE_PAIRS = 'pairs'

# Normalization choices
CORAL_SRC_NORM = 'coral-src' # align source "into" target
CORAL_TGT_NORM = 'coral-tgt' # align target non-landmarks into target landmarks
NONE_NORM  = 'none'
NORM_MODES = [ CORAL_SRC_NORM, CORAL_TGT_NORM, NONE_NORM ]
