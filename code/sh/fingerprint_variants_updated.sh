#!/usr/bin/bash

EXEC_NUCMER="/usr/bin/nucmer"
EXEC_SHOW_SNPS="/usr/bin/show-snps"
dir_data_fasta="/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants"
dir_data_computed="/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed"
fasta_reference="${dir_data_fasta}/MN908947.3.fasta"
file_all_snp_data="${dir_data_computed}/$(basename ${fasta_reference} .fasta)_all_snp.data"
file_all_snp_data_header="/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/column.header.csv"
file_all_metadata="/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/metadata.csv" # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/column.header.csv"

echo ${EXEC_NUCMER}
echo ${EXEC_SHOW_SNPS}
echo ${dir_data_fasta}
echo ${dir_data_computed}
echo ${fasta_reference}
echo ${file_all_snp_data}
echo ${file_all_snp_data_header}
ls -hals ${fasta_reference} ${file_all_snp_data_header} ${file_all_snp_data}
ls ${dir_data_fasta}/*.fasta | sort | uniq
echo

touch ${file_all_snp_data}
head -n 1 ${file_all_snp_data_header} > ${file_all_snp_data}
cat ${file_all_snp_data}

for file_query in $(ls ${dir_data_fasta}/*.fasta | sort | uniq); 
do 
fasta_query=${file_query}
prefix="${dir_data_computed}/$(basename ${fasta_reference} .fasta)_$(basename ${fasta_query} .fasta)"
tag_reference="$(basename ${fasta_reference} .fasta)"
tag_query="$(basename ${fasta_query} .fasta)"
file_delta="${prefix}.delta"
file_snp="${prefix}.snp"
file_snp_raw="${prefix}.snp_raw.txt"
file_tsv="${prefix}.tsv"
echo ${prefix} 
echo ${fasta_reference} 
echo ${fasta_query} 
echo ${file_delta} 
echo ${file_snp}
echo ${file_snp_raw}
echo ${file_tsv}
echo ${prefix} 
echo ${fasta_reference} 
echo ${fasta_query}
echo "${EXEC_NUCMER} --prefix=${prefix} ${fasta_reference} ${fasta_query}"
${EXEC_NUCMER} --prefix=${prefix} ${fasta_reference} ${fasta_query}
# done
${EXEC_SHOW_SNPS} -Clr ${file_delta} > ${file_snp_raw}
${EXEC_SHOW_SNPS} -HClr ${file_delta} > ${file_snp}
cat ${file_snp} | tr '|' ' ' | tr -s ' ' | tr -s ' ' | tr ' ' '\t' > ${file_tsv}
sed -i '1s/^/'${tag_query}'\t/; 2,$s/^/'${tag_query}'\t/' ${file_tsv}
sed -i '1s/^/'${tag_reference}'\t/; 2,$s/^/'${tag_reference}'\t/' ${file_tsv}
done



cat ${file_snp} | tr '|' ' ' | tr -s ' ' | tr -s ' ' | tr ' ' '\t' > ${file_tsv}

sed -i '1s/^/'${tag_query}'\t/; 2,$s/^/'${tag_query}'\t/' ${file_tsv}
sed -i '1s/^/'${tag_reference}'\t/; 2,$s/^/'${tag_reference}'\t/' ${file_tsv}

# done
cat $(ls ${dir_data_computed}/*.tsv | sort | uniq) >> ${file_all_snp_data}
touch ${file_all_snp_data}
cat ${file_all_snp_data}
# sed -i.bak '1,$s/^/${tag_reference}\t/' ${file_tsv}

# sed -i "s/$/\t${tag_reference}/" ${file_tsv}
# sed -i "s/$/\t${tag_query}/" ${file_tsv}

# ${EXEC_SHOW_SNPS} -HClr ${file_delta} > ${file_snp}
# cat ${file_snp} | grep -v "========================================================================================" | tr '|' ' ' | tr '[' ' ' | tr ']' ' ' | tail -n +2 | tr -s ' ' | tr -s ' ' | tr -s ' ' | tr ' ' '\t' > ${file_tsv}
# done
# awk 'NR==1 {print}  NR>1 {printf("%s\t%s\n", $0, "*") }' input.csv > newoutput.csv


# cat $(ls /home/bgilot/Downloads/SarsCovControls/*.snp | sort | uniq) | sort | uniq

# tail -n +4 $(ls /home/bgilot/Downloads/SarsCovControls/*.snp | sort | uniq) | uniq | less


# nucmer --prefix=ref_qry ref.fasta qry.fasta

# show-snps -Clr ref_qry.delta > ref_qry.snps

# https://mummer.sourceforge.net/manual/#snpdetection

#  2001  nucmer --prefix=MN908947.3_MT007544.1 /home/bgilot/Downloads/SarsCovControls/MN908947.3.fasta '/home/bgilot/Downloads/SarsCovControls/MT007544.1.fasta' 
#  2002  show-snps -Clr /home/bgilot/Downloads/SarsCovControls/MN908947.3_MT007544.1.delta > /home/bgilot/Downloads/SarsCovControls/MN908947.3_MT007544.1.snps
#  2003  history


# fasta_query="/home/bgilot/Downloads/SarsCovControls/EPI_ISL_6841980.fasta"
# ls /home/bgilot/Downloads/SarsCovControls/*.fasta | sort | uniq

# fasta_reference="/home/bgilot/Downloads/SarsCovControls/MN908947.3.fasta"
# fasta_query="/home/bgilot/Downloads/SarsCovControls/EPI_ISL_6841980.fasta"
# prefix="$(basename ${fasta_reference} .fasta)_$(basename ${fasta_query} .fasta)"
# file_delta="$(dirname ${fasta_reference})/${prefix}.delta"
# file_snp="$(dirname ${fasta_reference})/${prefix}.snp"

# echo ${prefix} 
# echo ${fasta_reference} 
# echo ${fasta_query} 
# echo ${file_delta} 

# nucmer --prefix=${prefix} ${fasta_reference} ${fasta_query} 