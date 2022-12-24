#! /bin/bash

# usage
# test.sh A0A024RBG1_LLR.csv

echo 
filename=`echo $1`
DIRECT=/ycga-gpfs/scratch60/lifton/jc2545/esm1b/content/ALL_hum_isoforms_ESM1b_LLR
output=/home/jc2545/scratch60/hj/esm1b_hj/esm1b_results
#output=/home/jc2545/scratch60/hj/esm1b_hj
uniprot=`echo "${filename}" | sed 's/_LLR.csv//g'`
hgnc_uniprotid=/home/jc2545/scratch60/hj/esm1b_hj/hgnc_uniprotid_new.txt

### col1 : aa position
### col2 : aa ref 
cat ${DIRECT}/${filename} | head -n 1 | tr "," "\n" | sed 's/ /\t/g' | awk '{print $2"\t"$1}' | sed '1d'  > ${output}/${uniprot}_temp.txt

### iter 20 
for i in `seq 1 20`
do
    cat ${output}/${uniprot}_temp.txt >> ${output}/${uniprot}_temp2.txt 
done


### alternative aa and Score
cat ${DIRECT}/${filename} | sed '1d' > ${output}/${uniprot}_temp3.txt

### func1 used in each score for alternative aa 
func1() {
    # set array
        ar=()
        sentence=`echo "$1" | sed 's/\,/ /g'`
    #    echo "${sentence}"

        for word in $sentence
        do
                ar+=($word)
        done
    #    echo "${ar[0]}"

    # array length
        arlength=${#ar[@]}
        aalength=`expr \( ${arlength} - 1 \)`

        altaa=${ar[0]} # assign alt aa to altaa
        unset ar[0] # remove alt aa (index 0)  

    # assign each item to aa_score and iter
    # col1 : altaa
    # col2 : Score for aa change ref to altaa  
        for aa_score in "${ar[@]}"
        do
                echo -e "${altaa}\t${aa_score}" >> ${output}/${uniprot}_temp4.txt
        done
}


touch ${output}/${uniprot}_temp4.txt

while read line
do
    func1 $line
done < ${output}/${uniprot}_temp3.txt

paste -d "\t" ${output}/${uniprot}_temp2.txt ${output}/${uniprot}_temp4.txt > ${output}/${uniprot}_temp5.txt


### expand Gene symbol

## handling an exception
# if the number of matched lines is equal or greater than 2,
# the line not including "-" is assigned to HGNC variable

grepline=`grep "${uniprot}" $hgnc_uniprotid | wc -l | cut -d ' ' -f 1`
echo "${grepline}"
if [ "$grepline" -ge 2 ]
then 
	HGNC=`grep "${uniprot}" $hgnc_uniprotid | awk '{print $1}' | grep -v '-'`
	#echo "${HGNC}"
else # 0 or 1
	HGNC=`grep "${uniprot}" $hgnc_uniprotid | awk '{print $1}'`
	#echo "${HGNC}"
fi

# used isoform_LLR.csv in a cammand (if statement is true), remove *temp.txt files. 
if [ "$HGNC" == "" ]
then
	# remove temp files
	rm -rf ${output}/${uniprot}_temp.txt ${output}/${uniprot}_temp2.txt ${output}/${uniprot}_temp3.txt ${output}/${uniprot}_temp4.txt ${output}/${uniprot}_temp5.txt
else 
	linenumb=`wc -l ${output}/${uniprot}_temp5.txt | cut -d ' ' -f 1`
	#echo "${linenumb}"

	touch ${output}/${uniprot}_temp6.txt

	for i in `seq 1 ${linenumb}`
	do 
		echo "${HGNC}" >> ${output}/${uniprot}_temp6.txt
	done

	paste -d "\t" ${output}/${uniprot}_temp6.txt ${output}/${uniprot}_temp5.txt > ${output}/${HGNC}.txt

	# remove temp files
	rm -rf ${output}/${uniprot}_temp.txt ${output}/${uniprot}_temp2.txt ${output}/${uniprot}_temp3.txt ${output}/${uniprot}_temp4.txt ${output}/${uniprot}_temp5.txt ${output}/${uniprot}_temp6.txt
fi
