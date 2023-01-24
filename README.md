# esm1b_annotation
this repository includes scripts for annotation WES or WGS data with esm1b database.

### make idmapping file 
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/  
```shell
grep 'NM_' UP000005640_9606.idmapping | awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$3}' | awk 'BEGIN{FS="\t";OFS="\t"} {sub(".[0-9]*$","",$4); print $0}' | sed 's/-1//g' > 20230122_UP000005640_9606.idmapping_refseq_NM_hj.txt
```
### convert all isoforms 
```shell
awk '{print "source ~/.bashrc; sh 20230121_ruddle_esm1b_conversion_Jan202311_hj.sh -i ./content/ALL_hum_isoforms_ESM1b_LLR/"$1" -o esm1b_converted"}' filelist.txt > command.sh
dsq --job-file command.sh -J hj_esm1b -c 1 --mem=10G -t 120:00:00
sbatch dsq-command-2023-01-21.sh
```


```shell
sh 20230115_idmapping_refseq_Jan152023_hj.sh
```
20230115_idmapping_table.txt  
20230115_idmapping_refseq_Jan152023_runningtime_hj.txt  
20230115_idmapping_table_esm1b_match.txt  

### test
20230122_esm1b_annotation_Jan222023_hj.py is for variant table.  
20230124_esm1b_annotation_Jan222024_hj.py is for *vcf.gz file.  

```python
(py3.7)[jc2545@c15n08 esm1b_hj]$ python 20230122_esm1b_annotation_Jan222023_hj.py thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno_0.001_REVEL_plDiff8_dom_variant_removeDup.txt esm1b_converted uniprot_id/20230122_UP000005640_9606.idmapping_refseq_NM_hj.txt

(py3.7)[jc2545@c25n09 esm1b_hj]$ python 20230124_esm1b_annotation_Jan222024_hj.py ../thyroiditis/thy_genome_calls_pass_n17_decomposed_normalized_anno.hg19_multianno.vcf.gz ./esm1b_converted uniprot_id/20230122_UP000005640_9606.idmapping_refseq_NM_hj.txt
```

