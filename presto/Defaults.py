"""
Default parameters
"""
# Info
__author__ = 'Jason Anthony Vander Heiden, Namita Gupta'
from presto import __version__, __date__

# Annotation parameters
default_delimiter = ('|', '=', ',')
default_separator = default_delimiter[2]

# Commandline argument defaults
choices_coord = ['illumina', 'solexa', 'sra', '454', 'presto']
default_coord = 'presto'
default_out_args = {'log_file':None,
                    'delimiter':default_delimiter,
                    'separator':default_separator,
                    'out_dir':None,
                    'out_name':None,
                    'out_type':None,
                    'failed':True}

# Fields
default_barcode_field = 'BARCODE'
default_primer_field = 'PRIMER'
default_cluster_field = 'CLUSTER'

# External applications
default_muscle_exec = 'muscle'
default_usearch_exec = 'usearch'
default_blastn_exec = 'blastn'
default_blastdb_exec = 'makeblastdb'

# Sequence sets
default_gap_chars = set(['-', '.'])
default_mask_chars = set(['n', 'N'])
default_missing_chars = set(['-', '.', 'n', 'N'])
default_missing_residues = set(['.', '-', 'x', 'N'])

# Consensus defaults
default_min_freq = 0.6
default_min_qual = 20

# Primer defaults
default_gap_penalty = (1, 1)
default_max_error = 0.2
default_max_len = 50
default_start = 0