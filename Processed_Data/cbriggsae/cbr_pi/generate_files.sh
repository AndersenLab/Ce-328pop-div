# get strain list 
bcftools query -l cbriggsae_cutter.vcf.gz > cbr_strains.txt

# make windows
bedtools makewindows -g cbr_chrom_len -w 1000 > windows.bed

#
while read p; do
	bcftools view -s $p cbriggsae_cutter.vcf.gz |\
	bcftools filter -i 'GT="alt"' -Ov |\
	bedtools coverage -a windows.bed -b stdin -counts > variant_counts/$p.variant_counts.txt
done <cbr_strains.txt