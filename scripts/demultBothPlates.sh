#!/bin/bash

# created quiime1 as so, need earlier version of numpy for one of these things to work
# conda create -n qiime1 numpy=1.10 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda

# move out of scripts directory and into base
cd ..

echo "clear directories, if they exist"
rm -rf raw096seq
rm -rf proc096
rm -rf for_dada2

echo "make directories"
# for intermediatly processed sequence data
mkdir proc096 
mkdir proc096/plate65
mkdir proc096/plate66
# for demultiplexed and trimmed sequence data
mkdir for_dada2 

echo "extract sequence data from tar files"
# makes raw096seq file
tar -xjvf data/raw096seq.tar.bz2

echo "--turn on qiime1--"
source activate qiime1

echo "--save information about qiime1--"
print_qiime_config.py > proc096/qiime_config_output.txt

# Plate 65 Work

echo "--process plate 65--"

echo "--validate mapping files--"
validate_mapping_file.py -m data/mapping_file_096_plate65.tsv -o proc096/plate65/mapping_val/

echo "--convert fasta and qual file to fastaq file--"
time convert_fastaqual_fastq.py -f raw096seq/Plate\ 65/raw/data/1.TCA.454Reads.fna -q raw096seq/Plate\ 65/raw/data/1.TCA.454Reads.qual -o proc096/plate65/fastq

echo "--seperate barcodes from sequences--"
time extract_barcodes.py -f proc096/plate65/fastq/1.TCA.454Reads.fastq -o proc096/plate65/sepbc -l 6

echo "--demultiplex--"
time split_libraries_fastq.py -i proc096/plate65/sepbc/reads.fastq -b proc096/plate65/sepbc/barcodes.fastq -m proc096/plate65/mapping_val/mapping_file_096_plate65.tsv_corrected.txt --barcode_type 6 -o proc096/plate65/demult --phred_offset 33 --start_seq_id 100000 --store_demultiplexed_fastq

echo "--make one fasta file for each sample--"
time split_sequence_file_on_sample_ids.py -i proc096/plate65/demult/seqs.fastq -o proc096/plate65/split --file_type fastq

echo "--move those fasta files into a convenient directory--"
time cp proc096/plate65/split/* for_dada2

# Plate 66 Work

echo "--process plate 66--"

echo "--validate mapping files--"
validate_mapping_file.py -m data/mapping_file_096_plate66.tsv -o proc096/plate66/mapping_val/

echo "--convert fasta and qual file to fastaq file--"
time convert_fastaqual_fastq.py -f raw096seq/Plate\ 66/raw/data/1.TCA.454Reads.fna -q raw096seq/Plate\ 66/raw/data/1.TCA.454Reads.qual -o proc096/plate66/fastq

echo "--seperate barcodes from sequences--"
time extract_barcodes.py -f proc096/plate66/fastq/1.TCA.454Reads.fastq -o proc096/plate66/sepbc -l 6

echo "--demultiplex--"
time split_libraries_fastq.py -i proc096/plate66/sepbc/reads.fastq -b proc096/plate66/sepbc/barcodes.fastq -m proc096/plate66/mapping_val/mapping_file_096_plate66.tsv_corrected.txt --barcode_type 6 -o proc096/plate66/demult --phred_offset 33 --start_seq_id 100000 --store_demultiplexed_fastq

echo "--make one fasta file for each sample--"
time split_sequence_file_on_sample_ids.py -i proc096/plate66/demult/seqs.fastq -o proc096/plate66/split --file_type fastq

echo "--move those fasta files into a convenient directory--"
time cp proc096/plate66/split/* for_dada2

echo "turn off qiime"
source deactivate
