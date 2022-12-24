#! /bin/bash

echo
vt=`echo $1`

DIRECT=/home/jc2545/scratch60/hj/esm1b_hj/test_esm1b_annotation

# make header file
cat ${DIRECT}/${vt} | head -n 1 | tr '\t' '\n' > ${DIRECT}/header.txt

# assign AAChange filed number to AAchangeline
AAchangeline=`cat ${DIRECT}/header.txt | grep -n "AAChange.refGeneWithVer" | cut -d':' -f1`

# if AAchange field are included more then one time, 
# use the last line of the AAchangeline variable. 

linenumb=`echo "${AAchangeline}" | wc -l | cut -d ' ' -f1`

if [ "${linenumb}" -gt 1 ]
then
	AAchangelastline=`echo "${AAchangeline}" | tail -n 1`
else
	:
fi

# select the AAChange field and remove header
cat ${DIRECT}/${vt} | awk -F'\t' '{print $'"${AAchangelastline}"'}' | sed '1d' > ${DIRECT}/aachange.col.txt

touch ${DIRECT}/aachange.col.table_temp.txt

# write header
echo -e "gene\trefaa\taltaa\tposition" >> ${DIRECT}/aachange.col.table_temp.txt

pchange_split() {
        # set array
        ar=()
        # grep AA from the argument
        sequence=`echo $1 | grep -o '[A-Z]*'`
        #echo "${sequence}"
        for AA in $sequence
        do
                ar+=($AA)
                #echo "${ar[@]}"
                #echo "next"
        done
        Refaa=${ar[0]}
        Altaa=${ar[1]}
        Position=`echo $1 | grep -o '[0-9]*'`
        Gene=`echo $2`
        # write down
        echo -e "${Gene}\t${Refaa}\t${Altaa}\t${Position}" >> ${DIRECT}/aachange.col.table_temp.txt
}

while read line
do
        if [ "${line}" == "." ]
        then
                echo -e ".\t.\t.\t." >> ${DIRECT}/aachange.col.table_temp.txt
        else
                genename=`echo "${line}" | cut -d ':' -f1`
                #echo "${genename}"
                # Because annovar annotation does not present canonical transcripts,
                # the last transcript information is used to make protein change table.
                pchange=`echo "${line}" | tr ':' '\n' | tail -n 1`
                #echo "$pchange"
                # grep protein level information of missense mutatation
                missense=`echo "${pchange}" | grep 'p\.[A-Z][0-9]*[A-Z]$'`
                #echo "${missense}"
                if [ "${missense}" == "" ]
                then
                        echo -e "${genename}\t.\t.\t." >> ${DIRECT}/aachange.col.table_temp.txt
                else
                        pchange_split ${missense} ${genename}
                fi
        fi
done < ${DIRECT}/aachange.col.txt

cp ${DIRECT}/aachange.col.table_temp.txt ${DIRECT}/aachange.table.txt
rm -rf ${DIRECT}/header.txt ${DIRECT}/aachange.col.txt ${DIRECT}/aachange.col.table_temp.txt