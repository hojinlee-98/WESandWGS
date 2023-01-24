#! /bin/bash


print_usage() {
        echo "$0 -i [input] -o [outdir] -h"
        echo "<options>"
	echo "  -i [input] : /path/to/input"
        echo "  -o [output] : /path/to/outdir"
        echo "  -h : show this message"
}

while getopts i:o:h opts; do
        case ${opts} in
		i) input=$OPTARG
			;;
                o) outdir=$OPTARG
                        ;;
                h) print_usage
                        exit;;
        esac
done

uniprot=`basename ${input} | sed 's/.csv//g'`

### col1 : aa position
### col2 : aa ref 
cat ${input} | head -n 1 | tr "," "\n" | sed 's/ /\t/g' | awk '{print $2"\t"$1}' | sed '1d'  > ${outdir}/${uniprot}_temp.txt

### iter 20 times 
for i in `seq 1 20`
do
    cat ${outdir}/${uniprot}_temp.txt >> ${outdir}/${uniprot}_temp2.txt 
done


### alternative aa and Score
cat ${input} | sed '1d' > ${outdir}/${uniprot}_temp3.txt

### func1 used in each score for alternative aa 
esm1b_converter() {
    # set array
        ar=()
        sentence=`echo "$1" | sed 's/\,/ /g'`
    #    echo "${sentence}"

        for word in ${sentence}
        do
                ar+=(${word})
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
		# if this script is run on linux, use option -e for echo command
                echo "${altaa}\t${aa_score}" >> ${outdir}/${uniprot}_temp4.txt
        done
}


touch ${outdir}/${uniprot}_temp4.txt

while read line
do
    esm1b_converter ${line}
done < ${outdir}/${uniprot}_temp3.txt

paste -d "\t" ${outdir}/${uniprot}_temp2.txt ${outdir}/${uniprot}_temp4.txt > ${outdir}/${uniprot}.table.txt

# remove temp files 
rm -rf ${outdir}/${uniprot}_temp*
