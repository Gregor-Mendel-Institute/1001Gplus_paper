#!/bin/bash

# Check if the db_file argument is present
if [ -z "$1" ] || [ -z "$2" ]; then
  echo "Usage: $0 <db_file> <query_file>"
  exit 1
fi

db_file="$1"
query_file="$2"



# Run makeblastdb with the provided db_file parameter
# Check if the database already exists
if [ ! -e "${db_file}.nhr" ]; then
  # If the database does not exist, create it using makeblastdb
  makeblastdb -in "$db_file" -dbtype nucl
else
  echo "Database already exists. Skipping makeblastdb."
fi

# Extract the directory path and base file name of the query_file
query_dir=$(dirname "$query_file")
query_base=$(basename "${query_file%.fasta}")

# Construct the output file path (out_file) in the same directory as query_file
out_file="${query_dir}/out_${query_base}.txt"

blastn -db ${db_file} -query ${query_file} -out ${out_file} -outfmt "7 qseqid qstart qend sstart send pident length sseqid"


# Remove the database files
if [ -e "${db_file}.nhr" ]; then
  #echo "Removing database files..."
  rm "${db_file}".*
fi

echo "Done!"
