#! /bin/bash

cd /Volumes/hjdrive/db/esm1b

date > 20230115_idmapping_refseq_Jan152023_runningtime_hj.txt
touch 20230115_idmapping_table.txt
while read col1 col2 col3 col4
do
	cat /Volumes/hjdrive/db/hg19.ncbiRefSeqLink.txt |\
	awk '{if ($3 == "'${col2}'") {print "'${col1}'""\t""'${col2}'""\t""'${col3}'""\t""'${col4}'""\t"$2"\t"$1}}' >> 20230115_idmapping_table.txt
done < /Volumes/hjdrive/db/20230115_idmapping_selected_edit_hj.tab

date >> 20230115_idmapping_refseq_Jan152023_runningtime_hj.txt


awk -F '\t' '{print $4}' 20230115_idmapping_table.txt | sort | uniq > test1
awk -F ',' '{print $2}' ALL_hum_isoforms_ESM1b_LLR/contents_u_df.csv | sed '1d' | sort | uniq > test2
comm -1 -2 test1 test2 | sort | uniq > test3
while read line; do cat 20230115_idmapping_table.txt | awk -F '\t' '{if ($4 == "'${line}'") print $0}' >> test4; done < test3

rm -rf test1; rm -rf test2; rm -rf test3

mv test4 20230115_idmapping_table_esm1b_match.txt
