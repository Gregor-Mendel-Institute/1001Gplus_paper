#
# utils.sh - Helper functions for methylation data.  
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


# Error handling
set -xe
trap 'previous_command=$this_command; this_command=$BASH_COMMAND' DEBUG
trap 'handle_exit $? "$previous_command"' EXIT


# handle_exit()
#
# Function to handle the exit signal.
#
# Parameters:
#   $1: Exit code
#   $2: Last command
#
# Example:
#      handle_exit $? $last_command

function handle_exit () {
    local exit_code=$1
    local command=$2

    # Check if the exit code is non-zero
    if [ $exit_code -ne 0 ]; then
        echo "[ERROR] ${command}"
        echo "failed with exit code $1."
    fi
}


# load_assemblies_table()
#
# Read columns from assemblies table into arrays
#
# Parameters:
#   $1: assemblies table
#
# Example:
#   load_assemblies_table path/to/assemblies_1001Gp.tsv

declare -a acc_names
declare -a ass_names
declare -a jbrowse_enabled_arr

function load_assemblies_table () {
    local assemblies_table=$1

    acc_names=($(awk '{ print $3 }' $assemblies_table))
    ass_names=($(awk '{ print $5 }' $assemblies_table))
    jbrowse_enabled_arr=($(awk '{ print $9 }' $assemblies_table))
}


# list_assemblies_names()
#
# Print all assemblies names

function list_assemblies_names() {
    for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}

      if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi

      echo $ass_name
    done
}


# list_accession_names()
#
# Print all accession names

function list_accession_names() {
        for i in ${!ass_names[@]}
    do
      acc_name=${acc_names[$i]}
      ass_name=${ass_names[$i]}
      jbrowse_enabled=${jbrowse_enabled_arr[$i]}

      if [ "$ass_name" = "Assembly_Name" ] || [ "$ass_name" = "TAIR10" ] || [ "$ass_name" = "Araport11" ] || [ "$jbrowse_enabled" = "n" ]
      then
        continue 
      fi

      echo $acc_name
    done
}


# init_jbrowse()
#
# Create and initialize a new JBrowse2 directory.
# THIS DELETES $jbrowse_dir

function init_jbrowse() {
    echo "Setting up $jbrowse_dir ... "
    rm -rf $jbrowse_dir
    mkdir -p $jbrowse_dir

    # Always do install for upgrade
    sudo npm install -g @jbrowse/cli
    jbrowse create $jbrowse_dir

    echo "done."
}


# upload_jbrowse2_config()
#
# Copy JBrowse2 directory to web server.

function upload_jbrowse2_config() {
    # Rsync to public server
    dest="${user}@${hostname}:server_jbrowse_dir/"
    echo "Deploying to $dest ..."

    rsync -ravz $jbrowse_dir/* "$dest"
}


