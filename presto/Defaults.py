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
default_coord_choices = ['illumina', 'solexa', 'sra', '454', 'presto']
default_coord_type = 'presto'
default_barcode_field = 'BARCODE'
default_primer_field = 'PRIMER'
default_min_freq = 0.6
default_min_qual = 20
default_out_args = {'log_file':None,
                    'delimiter':default_delimiter,
                    'separator':default_separator,
                    'out_dir':None,
                    'out_name':None,
                    'out_type':None,
                    'failed':True}

# External applications
default_muscle_exec = r'/usr/local/bin/muscle'
default_usearch_exec = r'/usr/local/bin/usearch'
default_blastn_exec = r'/usr/bin/blastn'

# Sequence sets
default_missing_chars = set(['-', '.', 'N', 'n'])
default_missing_residues = set(['.', '-', 'X', 'x'])
