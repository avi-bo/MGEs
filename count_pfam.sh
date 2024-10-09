#!/bin/bash

transposon_pfam_list="path_to/transposon_PF_list"
plasmid_pfam_list="path_to/plasmid_PF_list"
defense_pfam_list="path_to/antiphage_PF_list"
phage_pfam_list="path_to/phage_PF_list"

pfam_dir="path_to_pfam_folder"

output_file="genomeId_categoryName_numPfam_totalPfam.tsv"
# >$output_file

while read -r id;
do
        if [ -f "$pfam_dir/$id.pfam" ]; then
                total_pfam=$(wc -l < "$pfam_dir/$id.pfam")
                for category_list in "$transposon_pfam_list" "$plasmid_pfam_list" "$defense_pfam_list" "$phage_pfam_list";
                do
                        category_name=$(basename "$category_list" | cut -d '_' -f 1)
                        category_pfam_count=$(grep -cf "$category_list" "$pfam_dir/$id.pfam")
                        echo -e "$id\t$category_name\t$category_pfam_count\t$total_pfam" >> "$output_file"
                done
        else
                echo "$id" >> no_pfam_file
        fi
done < "$1"
