#!/bin/bash

#
# pangenome_jbrowse2.sh - Creates a Pangenome JBrowse2 for the 1001G+ Project  
# 
# Authors:
#   Joffrey Fitz (joffrey.fitz@tuebingen.mpg.de)
# 
# Copyright (c) 2024 Max Planck Institute for Biology, 
#   Tübingen, Germany, https://weigelworld.org
#  
# This file is part of 1001GpJBrowse2.
#   
# 1001GpJBrowse2 is free software: you can redistribute it 
# and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of 
# the License, or (at your option) any later version.
# 
# 1001GpJBrowse2 is distributed in the hope that it will be 
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with 1001GpJBrowse2.  
# If not, see <https://www.gnu.org/licenses/>.
# 


# Include modules
. config.sh
. modules/utils.sh
. modules/refseq.sh
. modules/annotations.sh
. modules/rnaseq.sh
. modules/methylation.sh


# Load Asemblies configuration
load_assemblies_table $assemblies_table


# 1 Create JBrowse2 directory
init_jbrowse


# 2 Pangenome maximum consensus sequence

# 2.1 Compress and index pangenome reference sequence
process_refseq

# 2.2 Load refseq to JBrowse2
load_refseq

# 2.3 Copy fasta sequence and indices to web server
deploy_refseq


# 3 Annotations

# 3.1 Copy to work dir
#copy_annotations
# Copy without 22001f (MPI internal only)
_copy_annotations

# Fix 22001f (MPI internal only)
_fix_annotation_22001f

# 3.2 Sort, bgzip, and tabix
process_annotations

# 3.3 Generate annotation tracks
load_gff_annotations
# remove_annotation_tracks

# 3.4 Copy annotations and indices to web server
deploy_annotations


# 4 RNASeq

# 4.1 Download FR files from 1001genomes.org
get_rnaseq_data

# 4.2 Adjust chromosome identifiers, translate coordinates
process_rnaseq_data

# 4.3 Generate RNASeq tracks
load_rnaseq

# 4.4 Upload data to web server
deploy_rnaseq_data


# 5 Methylation data

# 5.1 Download methylation files from 1001genomes.org
get_methylation_data

# 5.2 Adjust chromosome identifiers, translate coordinates
process_methylation_data

# 4.3 Generate methylation tracks
load_methylation

# 5.3 Upload data to web server
upload_methylation_data


# 6 Upload JBrowse2 configuration to webserver
upload_jbrowse2_config

exit 0
