#
# rnaseq.sh - Helper functions for RNASeq data.  
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


# get_rnaseq_data()
#
# Download RNASeq data from 10001genomes.org

function get_rnaseq_data() {
    rm -rf $rnaseq_work_dir
    mkdir -p $rnaseq_work_dir

    echo -n "Downloading RNASeq data from $rnaseq_url_prefix ... "

    (cd $rnaseq_work_dir && 
        wget -nd -q -r -l1 --no-parent -A "RNASeq_merged.*.bw" $rnaseq_url_prefix)

    echo "done."
}


# process_rnaseq_data()
#
# Adjust chromosome identifier and translate coordinates.

function process_rnaseq_data() {
    for i in ${!ass_names[@]}
    do
        acc_name=${acc_names[$i]}
        ass_name=${ass_names[$i]}
        jbrowse_enabled=${jbrowse_enabled_arr[$i]}
        file_name=${file_names[$i]}
        color="#00497F"

        if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
        then
        continue 
        fi
        echo -n "Processing $ass_name ... "


        bedgraph_file="$rnaseq_work_dir/RNASeq_merged.${ass_name}.bedGraph"

        # Convert to bedGraph
        bigWigToBedGraph $rnaseq_work_dir/RNASeq_merged.${ass_name}.bw $bedgraph_file
          
        # Adjust chromosome identifiers
        # 1741_Chr1 -> PanGen_Chr_1
        perl -p -i -e "s/${ass_name}_Chr/PanGen_Chr_/g" $bedgraph_file

        # Translate coordinates
        ### TBD

        # Determine chromosome sizes
        ### USE pangenome sequence ###
        cut -f1,2 $assemblies_dir/${ass_name}.scaffolds_corrected.v2.1.fasta.gz.fai \
        > $rnaseq_work_dir/${ass_name}_chrom.sizes
        ### TO BE REMOVED
        perl -p -i -e "s/${ass_name}_Chr/PanGen_Chr_/g" $rnaseq_work_dir/${ass_name}_chrom.sizes

        # Convert to bigwig
        wigToBigWig $bedgraph_file \
        $rnaseq_work_dir/${ass_name}_chrom.sizes \
        $rnaseq_work_dir/RNASeq_merged.${ass_name}.pangenome.bw

        # Clean up 
        rm $rnaseq_work_dir/RNASeq_merged.${ass_name}.bw \
        $bedgraph_file

        echo "done."
    done
}


# load_rnaseq()
#
# Load annotations to jbrowse with the "jbrowse add-track" command.
# Track data is a remote file defined by variables in config.sh
#
# ${server_url}/${server_rnaseq_dir}/RNASeq_merged.${ass_name}.pangenome.bw

function load_rnaseq() {
    for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}
      file_name=${file_names[$i]}
      color="#00497F"

      if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi
      echo "$ass_name"
      
      # Add all merged track
      merged_file="RNASeq_merged.${ass_name}.bw"
      name="$acc_name RNASeq"
      url="${server_url}/${server_rnaseq_dir}/RNASeq_merged.${ass_name}.pangenome.bw"

      json=$( 
        jq -c --null-input \
        --arg merged_file $merged_file \
        --arg name "$name" \
        --arg url $url \
        --arg ass_name $ass_name \
        --arg color $color \
        --arg acc_name $acc_name \
          '{
              "type": "MultiQuantitativeTrack",
              "trackId": $merged_file,
              "name": $name,
              "category": [$acc_name],
              "assemblyNames": ["pangenome"],
              "adapter": {
                "type": "MultiWiggleAdapter",
                "subadapters": [
                  {
                    "type": "BigWigAdapter",
                    "name": $name,
                    "bigWigLocation": {
                      "uri": $url
                    },
                    "color": $color
                  }                  
                ]
              }
            }'
      )

      # echo $json | jq -c

      jbrowse remove-track $merged_file --out $jbrowse_dir
      jbrowse add-track-json "$json" --out $jbrowse_dir
      
    done
}


# deploy_rnaseq_data()
#
# Rsync pangenome RNASeq to server

function deploy_rnaseq_data() {
    # Rsync to public server
    dest="${user}@${hostname}:${server_fs_path}/${server_rnaseq_dir}/"
    echo "Deploying to $dest ..."
    rsync -ravz $rnaseq_work_dir/*pangenome.bw "$dest" 
}




