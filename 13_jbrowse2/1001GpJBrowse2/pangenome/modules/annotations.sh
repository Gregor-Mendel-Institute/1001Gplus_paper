#
# annotations.sh - Helper functions for GFF annotations.  
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


# copy_annotations()
#
# Copy annotations to work dir

function copy_annotations() {
    # Copy to work dir

    rm -rf $ann_work_dir
    mkdir -p $ann_work_dir

    for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}

      if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi

      echo -n "Copying $ass_name to $ann_work_dir ... "

      cp $ann_dir/${annotations_prefix}${ass_name}.gff $ann_work_dir

      echo "done."
    done
}


# _copy_annotations()
# 
# Copy annotations to work dir (MPI internal only). 22001f is skipped

function _copy_annotations() {
    # Copy to work dir

    rm -rf $ann_work_dir
    mkdir -p $ann_work_dir

    for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}

      if [ "$ass_name" = "22001f" ] || [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi

      echo -n "Copying $ass_name to $ann_work_dir ... "

      cp $ann_dir/${annotations_prefix}${ass_name}.gff $ann_work_dir

      echo "done."
    done
}


# _fix_annotation_22001f()
# 
# Rename file and chromosome indentifiers

function _fix_annotation_22001f() {
    # Fix 22001 weirdo
    cp $ann_dir/${annotations_prefix}220011.gff $ann_work_dir/${annotations_prefix}22001f.gff
    perl -p -i -e "s/Consensus_220011/Consensus_22001f/g" $ann_work_dir/${annotations_prefix}22001f.gff
    perl -p -i -e "s/\.220011/\.22001f/g" $ann_work_dir/${annotations_prefix}22001f.gff
}


# process_annotations()
#
# Sort, bgzip, and tabix

function process_annotations() {
    # Sort GFF, bgzip, and tabix
    for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}

      if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi
      echo -n "Processing ${annotations_prefix}${ass_name}.gff ... "

      gt gff3 -sortlines -tidy -retainids $ann_work_dir/${annotations_prefix}${ass_name}.gff \
        > $ann_work_dir/${annotations_prefix}${ass_name}.sorted.gff 2>> $ann_work_dir/errors.txt
      bgzip $ann_work_dir/${annotations_prefix}${ass_name}.sorted.gff
      mv $ann_work_dir/${annotations_prefix}${ass_name}.sorted.gff.gz $ann_work_dir/${annotations_prefix}${ass_name}.gff.gz
      tabix $ann_work_dir/${annotations_prefix}${ass_name}.gff.gz

      echo "done."
    done

    # Add TAIR10
    cp $ann_dir/${annotations_prefix}tair.gff $ann_work_dir/
    ~/bin/gt gff3 -sortlines -tidy -retainids $ann_work_dir/${annotations_prefix}tair.gff \
      > $ann_work_dir/${annotations_prefix}tair.sorted.gff 2>> $ann_work_dir/errors.txt
    bgzip $ann_work_dir/${annotations_prefix}tair.sorted.gff
    mv $ann_work_dir/${annotations_prefix}tair.sorted.gff.gz $ann_work_dir/${annotations_prefix}tair.gff.gz
    tabix $ann_work_dir/${annotations_prefix}tair.gff.gz
}


# load_gff_annotations()
#
# Load annotations to jbrowse with the "jbrowse add-track" command.
# Track data is a remote file defined by variables in config.sh
#
# ${server_url}${server_url_annotations_path}${server_url_annotations_prefix}${ass_name}.gff.gz

function load_gff_annotations() {
    for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}

      if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi
      echo "$ass_name ($acc_name)"

      # Assemble URL
      url="${server_url}${server_url_annotations_path}${server_url_annotations_prefix}${ass_name}.gff.gz"
      echo "URL: $url"

      jbrowse add-track \
        $url \
        --assemblyNames pangenome \
        --category "$acc_name" \
        --name "$acc_name Genes" \
        --out $jbrowse_dir \
        --force
    done

    # Add TAIR10
    url="${server_url}${server_url_annotations_path}${annotations_prefix}tair.gff.gz"
      

    jbrowse add-track \
      $url \
      --assemblyNames pangenome \
      --category "TAIR10" \
      --name "TAIR10 Genes" \
      --out $jbrowse_dir \
      --force

    jbrowse text-index \
      --assemblies=pangenome \
      --attributes=Name \
      --tracks "genes_v05_pangen_tair.gff" \
      --force \
      --target=$jbrowse_dir

}


# remove_annotation_tracks()
#
# Remove annotation tracks using the "jbrowse remove-track" command.
# Note: The track ID is the filename without .gz suffix

function remove_annotation_tracks() {
    local trackID

    for i in ${!ass_names[@]}
    do
      ass_name=${ass_names[$i]}
      trackID="${annotations_prefix}${ass_name}.gff"

      echo "Removing $trackID"
      # jbrowse remove-track $trackID --out $jbrowse_dir
    done

    # And TAIR
    trackID="${annotations_prefix}tair.gff"
    # jbrowse remove-track $trackID --out $jbrowse_dir
}


# deploy_annotations_data()
#
# Rsync annotations indices to server

function deploy_annotations() {
    # Rsync to public server
    dest="${user}@${hostname}:${server_fs_path}/${server_annotations_dir}/"
    echo "Deploying to $dest ..."
    rsync -ravz $ann_work_dir/*.gz* "$dest" 
}
