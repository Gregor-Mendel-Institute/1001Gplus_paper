#
# methylation.sh - Helper functions for methylation data.  
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


# get_methylation_data()
#
# Download methylation data from 10001genomes.org

function get_methylation_data() {
    rm -rf $methylation_work_dir
    mkdir -p $methylation_work_dir

    echo -n "Downloading RNASeq data from $methylation_url_prefix ... "

    (cd $methylation_work_dir && 
        wget -nd -q -r -l1 --no-parent -A "*.bedgraph" $methylation_url_prefix)

    echo "done."
}


# process_methylation_data()
#
# Remove chromosomes C and M, adjust chromosome identifiers, 
# translate coordinates.

function process_methylation_data() {
    ass_name=${meth_ass[$i]}
    echo "Processing $ass_name ... "

    # Determine chromosome sizes
    ### USE pangenome sequence ###
    cut -f1,2 $assemblies_dir/${ass_name}.scaffolds_corrected.v2.1.fasta.gz.fai \
      > $methylation_work_dir/${ass_name}_chrom.sizes
    ### TO BE REMOVED
    perl -p -i -e "s/${ass_name}_Chr/PanGen_Chr_/g" $methylation_work_dir/${ass_name}_chrom.sizes

    for meth_type in "CG" "CHG" "CHH"
    do
        echo -en "\t$meth_type ... "

        # Remove ChrM/ChrC
        grep -v "Chr[MC]" \
        $methylation_work_dir/${ass_name}.${meth_type}meth.bedgraph \
          > $methylation_work_dir/${ass_name}.${meth_type}meth.noChrCM.bedgraph

        # Adjust chromosome identifiers
        # 1741_Chr1 -> PanGen_Chr_1
        perl -p -i -e "s/${ass_name}_Chr/PanGen_Chr_/g" \
          $methylation_work_dir/${ass_name}.${meth_type}meth.noChrCM.bedgraph

        # Translate coordinates
        ### TBD        

        bedGraphToBigWig $methylation_work_dir/${ass_name}.${meth_type}meth.noChrCM.bedgraph \
        $methylation_work_dir/${ass_name}_chrom.sizes \
        $methylation_work_dir/${ass_name}.${meth_type}meth.bw
        echo "done."
    done
}


# load_methylation()
# 
# Load annotations to jbrowse with the "jbrowse add-track" command.
# Track data is a remote file defined by variables in config.sh
#
# ${server_url}/${server_methylation_dir}/${ass_name}.CGmeth.bw"
# ${server_url}/${server_methylation_dir}/${ass_name}.CHGmeth.bw"
# ${server_url}/${server_methylation_dir}/${ass_name}.CHHmeth.bw"

function laod_methylation() {
    for i in ${!meth_ass[@]}
    do
        ass_name=${meth_ass[$i]}
        acc_name=${meth_acc[$i]}
        echo "$ass_name"

        color_aqua="#0096FF"
        color_clover="#008F00"
        color_tangerine="#FF9300"

        CG_url="${server_url}/${server_methylation_dir}/${ass_name}.CGmeth.bw"
        CHG_url="${server_url}/${server_methylation_dir}/${ass_name}.CHGmeth.bw"
        CHH_url="${server_url}/${server_methylation_dir}/${ass_name}.CHHmeth.bw"

        track_id="${ass_name}_methylation"

        json=$( 
            jq -c --null-input \
            --arg track_id $track_id \
            --arg CG_url $CG_url \
            --arg CHG_url $CHG_url \
            --arg CHH_url $CHH_url \
            --arg ass_name $ass_name \
            --arg acc_name $acc_name \
            --arg color_aqua $color_aqua \
            --arg color_clover $color_clover \
            --arg color_tangerine $color_tangerine \
              '{
                  "type": "MultiQuantitativeTrack",
                  "trackId": $track_id,
                  "name": "Methylation",
                  "category": [$acc_name],
                  "assemblyNames": ["pangenome"],
                  "adapter": {
                    "type": "MultiWiggleAdapter",
                    "subadapters": [
                      {
                        "type": "BigWigAdapter",
                        "name": "CG",
                        "bigWigLocation": {
                          "uri": $CG_url
                        },
                        "color": $color_aqua
                      },
                      {
                        "type": "BigWigAdapter",
                        "name": "CHG",
                        "bigWigLocation": {
                          "uri": $CHG_url
                        },
                        "color": $color_clover
                      },
                      {
                        "type": "BigWigAdapter",
                        "name": "CHH",
                        "bigWigLocation": {
                          "uri": $CHH_url
                        },
                        "color": $color_tangerine
                      }
                      
                    ]
                  }
                }'
        )

        # echo $json | jq -c

        jbrowse remove-track $track_id --out $jbrowse_dir
        jbrowse add-track-json $json --out $jbrowse_dir
    done
}


# deploy_methylation_data()
#
# Rsync pangenome methylation data to server

function deploy_rnaseq_data() {
    # Rsync to public server
    dest="${user}@${hostname}:${server_fs_path}/${server_methylation_dir}/"
    echo "Deploying to $dest ..."
    rsync -ravz $rnaseq_work_dir/*.bw "$dest" 
}



