#!/bin/bash

kos_list=${1}
kos_w_desc=${2}

kegg_genes=$(curl -# https://rest.kegg.jp/list/ko)

while read line
do
	hit=$(echo -e "${kegg_genes}" | grep "${line}")
	echo -e "${hit}" | tee -a "${kos_w_desc}"
done < ${kos_list}