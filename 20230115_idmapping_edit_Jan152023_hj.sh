#! /bin/bash

cd /Volumes/hjdrive/db
date > 20230115_idmapping_edit_Jan152023_runningtime_hj.txt

cat /Volumes/hjdrive/db/idmapping_selected.tab | \
grep "NP_[0-9|.]*" | \
awk 'BEGIN{FS="\t";ORS="\n";OFS="\t"}
{
split($2, human, "_")
if (($10 != "") && (human[2] == "HUMAN")) # col10 means UniRef50 id
	{
	gsub(" ", "", $4)
	split($4, a, ";") # split fourth column using ";", and then assign them to array a
	orig=$10
	gsub("UniRef50_", "", $10)
	for(i in a) {
	  split(a[i], b, "_")
	  if (b[1] == "NP") {print $2,a[i],orig,$10} # if a[1] is "NP", it means they have RefSeq NP id
	  }
	}
}' > /Volumes/hjdrive/db/20230115_idmapping_selected_edit_hj.tab

date >> 20230115_idmapping_edit_Jan152023_runningtime_hj.txt
