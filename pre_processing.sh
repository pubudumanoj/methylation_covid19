#this is the first script you need to run

#first select vcf files from WGS processing samples set
#put them into a vcf_selected folder
cd vcf_selected


module load mugqic/bedops/v2.4.35
module load mugqic/bedtools/2.27.0
mkdir ../filtered
for file in *
  do
zless $file | vcf2bed --do-not-split --deletions > ../filtered/$file.deletions.bed
done 

for file in *
  do
zless $file | vcf2bed --do-not-split > ../filtered/$file.snvs.bed
done 
cd ..

awk '{print $0}' filtered/*.bed | sort -k1,1 -k2,2n | bedtools merge -i - | awk -v OFS="\t" '{print $1,$2,$3}'  > wgs.snps.all.samples.union.txt

mkdir DSS_input

#prepare data for DSS

Rscript prepare_DSS.R

cd DSS_input

for filename in *.txt; do
sed '/_/d' $filename | egrep "^chr[0-9XYM][0-9T]?" | awk -v OFS="\t" '{if(NR!=1) {print $1,$2-1,$2,$3,$4}}' | bedtools intersect -a - -b ../wgs.snps.all.samples.union.withmed8.txt -v | awk -v OFS="\t" 'BEGIN{print "chr","pos","N","X"} { print $1,$3,$4,$5}' > ../DSS_input_wgs_all_corrected/$filename
done

cd ..
mdkir wgs_snpwith_med8_bed

cd DSS_input_wgs_all_corrected
for filename in *.txt; do
awk -v OFS="\t" '{if(NR!=1) {print $1,$2,$2}}' $filename > ../wgs_snpwith_med8_bed/$filename
done

sbatch day0.2021.oct.slurm
sbatch day5_15.2021.oct.slurm