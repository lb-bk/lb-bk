
########################################################################################################################
# Documents and scripts were written by: Bing Li & Zhaofu Yang
# For manuscript: Bing, et al. (2025). "Adaptive introgression provides insights into the evolution and novel sympatric habitats of Ostrinia furnacalis and Ostrinia nubilalis" (draft). 
# email: lb2021@nwafu.edu.cn
# Zhaofu Yang Lab. 
########################################################################################################################


# Part 1 ABBA-BABA test 
# This analysis was performed using the python script 'ABBABABAwindows.py' ,'calculate_abba_baba.r' and Fst values for the same windows were calculated using the script 'popgenWindows.py'
# The script is sourced at the following URL
# https://github.com/simonhmartin/genomics_general
# https://github.com/palc/tutorials-1
######################################################################################################################## 
## The D-statistic was quantified at whole-genome level
python /home/pyraloidea/workshop_material/genomics_general-master/VCF_processing/parseVCF.py -i Evo.recodeidtihuan.vcf.gz | gzip > Evo.recodeidtihuan.geno.gz
spc1=AJFX56; spc2=AHC16; spc3=EHC86; spc4=Outgroup
python /home/pyraloidea/workshop_material/genomics_general-master/freq.py -g Evo.recodeidtihuan.geno.gz -p ${spc1} -p ${spc2} -p ${spc3} -p ${spc4} --popsFile id.txt --target derived | grep -v nan | gzip > Evo.recodeidtihuan.tsv.gz
Rscript calculate_abba_baba.r Evo.recodeidtihuan.tsv.gz /home/pyraloidea/LBSNPanalysis/genhuocheng/AJF56_AHC16_EHC86.abba_baba.txt ${spc1} ${spc2} ${spc3} ${spc4} chr_lengths.txt
### To identify specific genomic regions subjected to introgression by nonoverlapping windows
## step1 Zscore D fd fdm caculate
#100kb
python3 ABBABABAwindows.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P1 LFS -P2 TBS -P3 TB -O Outgroup --popsFile all10pop_8qa_3qu.txt -w 100000 -m 200 --T 70 -f phased --writeFailedWindows -o Quercus_all267_allchr.ABBABABA_LFS_TBS_TB_out.w100km200.csv.gz 
#50kb
python3 ABBABABAwindows.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P1 LFS -P2 TBS -P3 TB -O Outgroup --popsFile all10pop_8qa_3qu.txt -w 50000 -m 100 --T 70 -f phased --writeFailedWindows -o Quercus_all267_allchr.ABBABABA_LFS_TBS_TB_out.w50km100.csv.gz 
#10kb
python3 ABBABABAwindows.py -g Quercus_all267_allchr_flt2SNP.geno2.gz -P1 LFS -P2 TBS -P3 TB -O Outgroup --popsFile all10pop_8qa_3qu.txt -w 10000 -m 20 --T 70 -f phased --writeFailedWindows -o Quercus_all267_allchr.ABBABABA_LFS_TBS_TB_out.w10km20.csv.gz 
## step2 Fst caculate
Fst caculate
#100kb
python /home/pyraloidea/ig_analysis/analysis/analysis1/genomics_general/popgenWindows.py -w 100000 -p ASY -p ADK -p EDK -o Outgroup --popsFile id.txt -g NEvo.recode.geno.gz -f phased -T 30 --analysis {popFreq,popDist,popPairDist} -o /home/pyraloidea/LBSNPanalysis/Fst/sympatric/1_ADK_EDKw100k.csv.gz
50kb
python /home/pyraloidea/ig_analysis/analysis/analysis1/genomics_general/popgenWindows.py -w 50000 -p ASY -p ADK -p EDK -o Outgroup --popsFile id.txt -g NEvo.recode.geno.gz -f phased -T 30 --analysis {popFreq,popDist,popPairDist} -o /home/pyraloidea/LBSNPanalysis/Fst/sympatric/1_ADK_EDKw100k.csv.gz
10kb
python /home/pyraloidea/ig_analysis/analysis/analysis1/genomics_general/popgenWindows.py -w 10000 -p ASY -p ADK -p EDK -o Outgroup --popsFile id.txt -g NEvo.recode.geno.gz -f phased -T 30 --analysis {popFreq,popDist,popPairDist} -o /home/pyraloidea/LBSNPanalysis/Fst/sympatric/1_ADK_EDKw100k.csv.gz
genome
python vcftools --vcf /home/pyraloidea/LBSNPanalysis/Fst/NEvo.recode.vcf --weir-fst-pop IMDK.txt --weir-fst-pop XJGL.txt --out p_1_2—single

###F-branch test
# step1 obtain tree.txt file
Dsuite Dtrios Quercus_all284_allchr_flt2SNPnosigmiss2.vcf.gz LFS_DWX_DW.txt -t LFS_DWX_DW.nwk -o ./LFS_DWX_DW 
# step2 f-branch analysis
Dsuite Fbranch LFS_DWX_DW.nwk LFS_DWX_DW__tree.txt > LFS_DWX_DW_branch.txt
# step3 f-branch statistic be visualized using dtools.py script
python3 ./dtools.py ./LFS_DWX_DW_branch.txt ./LFS_DWX_DW.nwk 

###Fst heatmap comparison
vcftools --vcf /home/pyraloidea/VCFwenjian/1_31.recode.vcf --weir-fst-pop XJYN.txt --weir-fst-pop IMZH.txt --out XJYN_IMZH
vcftools --vcf /home/pyraloidea/VCFwenjian/filter_1_31.recode.vcf --weir-fst-pop IMDK.txt --weir-fst-pop XJGL.txt --out /home/pyraloidea/Fst_heatmap/neturaltongyuheatmap/filter_IMDK_XJGL

###The constructed haplotype network 
plink --vcf xxxx.vcf --maf 0.05 --geno 0.2 --recode vcf-iid -out xxxx-maf0.05 --allow-extra-chr
plink --vcf xxxx-maf0.05.vcf --indep-pairwise 100 50 0.2 -out xxxx-maf0.05-LD --allow-extra-chr --make-bed
plink --noweb --bfile xxxx-maf0.05-LD --blocks no-pheno-req --allow-extra-chr
java -Xss5m -Xmx4g -jar YourWayToBeagle/beagle.jar gt=xxxx-maf0.05.vcf out=phased chrom=Chr:start-end
vcf2phylip.py -i phased.vcf.gz

###PSMC
time bwa mem -t 12 genome.fa sample_R1.clean.u.fq sample_R2.clean.u.fq | samtools sort -@ 12 -m 12G > sample_sort.bam
samtools mpileup -C50 -uf genome.fa sample_sort.bam | bcftools call --threads 16 -c - | vcfutils.pl vcf2fq -d 10 -D 100 |pigz -p 12 >diploid.fq.gz
fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa
nohup psmc -p "4+25*2+4+6" -t15 -N30 -r4 -o diploid.psmc diploid.psmcfa &
splitfa diploid.psmcfa > split.psmcfa
for i in {1..100};do nohup psmc -b -p "6+25*2+6+8" -t15 -N30 -r4 -o round${i}.psmc split.psmcfa &
cat diploid.psmc round*.psmc >combine.psmc
psmc_plot.pl -p -g 1 -u 4.79e-9 combine combine.psmc
cat /home/pyraloidea/LBSNPanalysis/zhongqundongtai/ADK11/ADK11.psmc /home/pyraloidea/LBSNPanalysis/zhongqundongt
ai/AHS21/AHS21.psmc /home/pyraloidea/LBSNPanalysis/zhongqundongtai/APT1/APT1.psmc /home/pyraloidea/LBSNPanalysis/zhongqundongtai/AZH10/AZH10.psmc > /home/pyraloidea/LBSNPanalysis/zhongqundongtai/HB/combine.psmc 
/home/pyraloidea/miniconda3/envs/psmc/bin/psmc_plot.pl -g 0.5 -u 3.5e-9 -M 'ADK11,AHS21,APT1,AZH10' combine combine.psmc

###fastsimcoal2
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_south$ bcftools view -S north_southtiqu.txt /home/pyraloidea/Fastsimcoal2/1_31_recode_with_ids.vcf - Ov > north_south_sfs.vcf
(base) pyraloidea@shpc-4067-instance-b6HRjrZJ:~$ conda activate my_python
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~$ cd /home/pyraloidea/Fastsimcoal2/north_northeast
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast$ /home/pyraloidea/Fastsimcoal2/easySFS/easySFS.py -i north_northeast_sfs.vcf -p north_northeast_sfs.txt --preview > proj_flag
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_south$ /home/pyraloidea/Fastsimcoal2/easySFS/easySFS.py -i file.vcf -p pop.txt -a --proj=3,5 -o easySFS/
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast$ chmod +x /home/pyraloidea/fsc27_linux64/fsc27093
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/no_geneflow$ for i in {1..50}; do mkdir run$i; cp no_geneflow.tpl no_geneflow.est no_geneflow_jointMAFpop1_0.obs run$i"/"; cd run$i; /home/pyraloidea/fsc27_linux64/fsc27093 -t no_geneflow.tpl -e no_geneflow.est -m -0 -C 10 -n 100000 -L 40 -s 0 -M -q -c 5; cd ..; done
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/best$ /home/pyraloidea/Fastsimcoal2/north_northeast/fsc-selectbestrun.sh
(Ryuyan) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/best/ongoing_geneflow_M12$ /home/pyraloidea/Fastsimcoal2/xinjiang_north/calculateAIC.sh ongoing_geneflow_M12
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/best/ongoing_geneflow_M12$ /home/pyraloidea/fsc27_linux64/fsc27093 -i ongoing_geneflow_M12_boot.par -n 100 -j -m -s0 -x -I -q
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/1_no_geneflow/bestrun/no_geneflow_boot$ for i in {1..100}; do cp no_geneflow_boot.est no_geneflow_boot.tpl no_geneflow_boot_jointMAFpop1_0.obs no_geneflow_boot.pv no_geneflow_boot$i"/"; cd no_geneflow_boot$i; /home/pyraloidea/fsc27_linux64/fsc27093 -t no_geneflow_boot.tpl -e no_geneflow_boot.est -m -0 -C 10 -n 100000 -L 40 -s 0 -M -q -c 5; cd ..; done
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/1_no_geneflow/bestrun/no_geneflow_boot$ cat no_geneflow_boot1/no_geneflow_boot/no_geneflow_boot.pv
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/1_no_geneflow/bestrun/no_geneflow_boot$ echo
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/1_no_geneflow/bestrun/no_geneflow_boot$ for i in {2..100}; do tail -n 1 no_geneflow_boot$i/no_geneflow_boot/no_geneflow_boot.pv; echo; done
(my_python) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/Fastsimcoal2/north_northeast/best/ongoing_geneflow_M12/ongoing_geneflow_M12_boot$ bash pv.sh > pv.txt
python
import pandas as pd
import numpy as np
pv_info = pd.read_csv('/home/pyraloidea/Fastsimcoal2/north_northeast/best/ongoing_geneflow_M12/ongoing_geneflow_M12_boot/pv.txt', sep='\t')
pv_info.quantile([0.05,0.95]).to_csv('/home/pyraloidea/Fastsimcoal2/north_northeast/1_no_geneflow/bestrun/no_geneflow_boot/result2.txt',sep='\t')

########################################################################################################################
# Part 2 Sympatric Mantel and partial Mantel regression tests by R package 'vegan'
########################################################################################################################
library(vegan)
library(dplyr)
setwd('D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source')
PSIG<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/1PSIG_data/PSIG_matrix.txt", header = TRUE) %>% as.matrix
Fst<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/2Fst_data/Fst_matrix.txt", header = TRUE) %>% as.matrix
env<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/3env_data/allenv_matrix.txt", header = TRUE) %>% as.matrix
geo<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/4geo_distance_data/geo_distance_matrix.txt", header = TRUE) %>% as.matrix
weidu<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/5weidu_data/weidu_matrix.txt", header = TRUE) %>% as.matrix
jingdu<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/6jingdu_data/jingdu_matrix.txt", header = TRUE) %>% as.matrix
wendu<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/7wendu_data/wendu_matrix.txt", header = TRUE) %>% as.matrix
jiangshui<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/8jiangshui_data/jiangshui_matrix.txt", header = TRUE) %>% as.matrix
haiba<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_source/9haiba_data/haiba_matrix.txt", header = TRUE) %>% as.matrix
### The correlation between genetic distance (Fst) and environmental/geographic variables
#Mantel test
mantel(env, Fst)
mantel(geo, Fst)
mantel(weidu, Fst)
mantel(jingdu, Fst)
mantel(wendu, Fst)
mantel(jiangshui, Fst)
mantel(haiba, Fst)
### The correlation between shared introgression(PSIG) and environmental/geographic variables
#mantel test
mantel(env, PSIG)
mantel(geo, PSIG)
mantel(weidu, PSIG)
mantel(jingdu, PSIG)
mantel(wendu, PSIG)
mantel(jiangshui, PSIG)
mantel(haiba, PSIG)
########################################################################################################################
# Part 2 Allopatric Mantel and partial Mantel regression tests by R package 'vegan'
########################################################################################################################
library(vegan)
library(dplyr)
setwd('D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric')
PSIG<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/1PSIG_data/PSIG_matrix.txt", header = TRUE) %>% as.matrix
Fst<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/2Fst_data/Fst_matrix.txt", header = TRUE) %>% as.matrix
env<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/3env_data/env_matrix.txt", header = TRUE) %>% as.matrix
geo<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/4geo_data/geo_matrix.txt", header = TRUE) %>% as.matrix
weidu<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/5weidu_data/weidu_matrix.txt", header = TRUE) %>% as.matrix
jingdu<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/6jingdu_data/jingdu_matrix1.txt", header = TRUE) %>% as.matrix
wendu<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/7wendu_data/wendu_matrix.txt", header = TRUE) %>% as.matrix
jiangshui<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/8jiangshui_data/jiangshui_matrix1.txt", header = TRUE) %>% as.matrix
haiba<- read.table("D:/RRRRRR/Workplace/SNPanalysis/correlation/mantel_test_allopatric/9haiba_data/haiba_matrix.txt", header = TRUE) %>% as.matrix
mantel(env, PSIG)
mantel(geo, PSIG)
mantel(weidu, PSIG)
mantel(jingdu, PSIG, method = "pearson")
mantel(wendu, PSIG)
mantel(jiangshui, PSIG)
mantel(haiba, PSIG)
mantel.partial(env, PSIG, geo)
mantel.partial(geo, PSIG, env)
mantel.partial(weidu, PSIG, geo)
mantel.partial(jingdu, PSIG, geo)
mantel.partial(wendu, PSIG, geo)
mantel.partial(jiangshui, PSIG, geo)
mantel.partial(haiba, PSIG, geo)
mantel.partial(env, PSIG, Fst)
mantel.partial(geo, PSIG, Fst)
mantel.partial(weidu, PSIG, Fst)
mantel.partial(jingdu, PSIG, Fst)
mantel.partial(wendu, PSIG, Fst)
mantel.partial(jiangshui, PSIG, Fst)
mantel.partial(haiba, PSIG, Fst)
mantel(env, Fst)
mantel(geo, Fst)
mantel(weidu, Fst)
mantel(jingdu, Fst)
mantel(wendu, Fst,)
mantel(jiangshui, Fst)
mantel(haiba, Fst)
mantel.partial(env, Fst, geo)
mantel.partial(geo, Fst, env)
mantel.partial(weidu, Fst, geo)
mantel.partial(jingdu, Fst, geo)
mantel.partial(wendu, Fst, geo)
mantel.partial(jiangshui, Fst, geo)
mantel.partial(haiba, Fst, geo)

########################################################################################################################
# Part 3 Identifying local adaptation and adaptive introgression by linux software 'baypass'
########################################################################################################################
vcftools --vcf ACB.vcf  --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 --recode --recode-INFO-all --stdout | gzip -c > north_chr1_31.vcf.gz
bgzip -d ACB_chr2_11.vcf.gz
bcftools view -S samplelist_north.txt /home/pyraloidea/LBSNPanalysis/Evo.recodeidtihuan.vcf -Ov > north.vcf
awk 'BEGIN{OFS="\t"} /^#/ {print $0; next} {$3=$1"_"$2; print $0}' ACB_chr2_11.vcf > ACB_chr2_11_with_ids.vcf
vcftools --vcf /home/pyraloidea/4adaptive_analysis/ACB_chr2_11/ACB_chr2_11_with_ids.vcf --plink --out ACB_chr2_11_with_ids.plink
plink --file ACB_chr12_21_with_ids.plink --make-bed --chr-set 27 --noweb --out ACB_chr12_21_with_ids
awk '{print $2,"\t",$1}' /home/pyraloidea/4adaptive_analysis/latitude.txt > ACB_chr12_21_with_ids_metadata_POPIND.txt
PLINK=ACB_chr12_21_with_ids.plink.ped
cut -f 3- $PLINK > x.delete
paste ACB_chr12_21_with_ids_metadata_POPIND.txt x.delete > ACB_chr12_21_with_ids.plink.ped
rm x.delete
plink --file ACB_chr12_21_with_ids.plink --chr-set 27 --freq counts --family --out ACB_chr12_21_with_ids
tail -n +2 ACB_chr12_21_with_ids.frq.strat | awk '{ $9 = $8 - $7 } 1' | awk '{print $7,$9}' | tr "\n" " " | sed 's/ /\n/6; P; D' > ACB_chr12_21_with_ids_baypass.txt
/home/pyraloidea/outlier_analysis/analysis/baypass/baypass_public-v2.4/sources/g_baypass -npop 28 -gfile ./ACB_chr1_with_ids_baypass.txt -outprefix ACB_chr1_with_ids_baypass -nthreads 2
R
library("ape")
setwd("/home/pyraloidea/adaptive_analysis/latitude_chr22_31")
source("/home/pyraloidea/outlier_analysis/programs/BayPass/baypass_utils.R")
omega <- as.matrix(read.table("ACB_chr22_31_with_ids_baypass_mat_omega.out"))
pi.beta.coef <- read.table("ACB_chr22_31_with_ids_baypass_summary_beta_params.out", header = TRUE)
bta14.data <- geno2YN("ACB_chr22_31_with_ids_baypass.txt")
simu.bta <- simulate.baypass(omega.mat = omega, nsnp = 5000, sample.size = bta14.data$NN, beta.pi = pi.beta.coef$Mean, pi.maf = 0, suffix = "btapods")
q()
/home/pyraloidea/outlier_analysis/analysis/baypass/baypass_public-v2.4/sources/g_baypass -npop 28 -gfile G.btapods -outprefix G.btapods -nthreads 2
R
setwd("/home/pyraloidea/adaptive_analysis/latitude_chr1")
source("/home/pyraloidea/outlier_analysis/programs/BayPass/baypass_utils.R")
library("ape")
library("corrplot")
pod.xtx <- read.table("G.btapods_summary_pi_xtx.out", header = T)
pod.thresh <- quantile(pod.xtx$M_XtX ,probs = 0.99)
pod.thresh
q()
cat ACB_chr1_with_ids_baypass_summary_pi_xtx.out | awk '$4>4.655721' > baypass_outliers.txt
grep -v "^#" /home/pyraloidea/adaptive_analysis/ACB_chr1/ACB_chr1_with_ids.vcf | cut -f1-3 | awk '{print $0"\t"NR}' > ACB_chr1_with_ids_SNPs.txt
awk 'FNR==NR{a[$1];next} (($4) in a)' baypass_outliers.txt ACB_chr1_with_ids_SNPs.txt | cut -f3 > baypass_outlierSNPIDs.txt
wc -l baypass_outlierSNPIDs.txt
R
setwd("/home/pyraloidea/adaptive_analysis/sympatric_jiangshui/precipitation1_chr1")
metadata <- read.table("/home/pyraloidea/adaptive_analysis/latitude.txt", sep = "\t", header = FALSE)
str(metadata)
pop_metadata <- aggregate(V3 ~ V2, data = metadata, mean)
pop_metadata[, 2]
write(pop_metadata[, 2], "pop_mean_latitude.txt")
q()
/home/pyraloidea/outlier_analysis/analysis/baypass/baypass_public-v2.4/sources/g_baypass -npop 28 -gfile ACB_chr1_with_ids_baypass.txt -efile pop_mean_latitude.txt -scalecov -auxmodel -nthreads 2 -omegafile ACB_chr1_with_ids_baypass_mat_omega.out -outprefix ACB_chr1_with_ids_baypass_latitude
R
setwd("/home/pyraloidea/adaptive_analysis/latitude_chr1")
covaux.snp.res.mass <- read.table("ACB_chr1_with_ids_baypass_latitude_summary_betai.out", header = T)
covaux.snp.xtx.mass <- read.table("ACB_chr1_with_ids_baypass_latitude_summary_pi_xtx.out", header = T)
pdf("Baypass_plots.pdf")
layout(matrix(1:3,3,1))
plot(covaux.snp.res.mass$BF.dB.,xlab="Mass",ylab="BFmc (in dB)")
abline(h=20, col="red")
plot(covaux.snp.res.mass$M_Beta,xlab="SNP",ylab=expression(beta~"coefficient"))
plot(covaux.snp.xtx.mass$M_XtX, xlab="SNP",ylab="XtX corrected for SMS")
dev.off()
null device 
1 
q()
cat ACB_chr1_with_ids_baypass_precipitation1_summary_betai.out | awk '$6>15' > ACB_chr1_baypass_precipitation1_BF20.txt
wc -l ACB_chr1_baypass_precipitation1_BF20.txt
awk 'FNR==NR{a[$2];next} (($4) in a)' ACB_chr1_baypass_precipitation1_BF20.txt ACB_chr1_with_ids_SNPs.txt | cut -f3 > baypass_precipitation1_outlierSNPIDs.txt
comm -12 <(sort baypass_precipitation1_outlierSNPIDs.txt) <(sort baypass_outlierSNPIDs.txt) > double_outliers.txt
wc -l double_outliers.txt

##########################################################################################################################################################
# Part 4  Linkage disequilibrium, population recombination rates and select sweep ananlysis for sympatric populations
##########################################################################################################################################################
vcftools --vcf ACB_with_ids.vcf --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 
--chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 
--recode --recode-INFO-all --stdout | gzip -c > ACB_with_ids.vcf.gz
vcftools --vcf /home/pyraloidea/4adaptive_analysis/ACB_chr2_11/ACB_chr2_11_with_ids.vcf --plink --out ACB_chr2_11_with_ids.plink
plink --vcf ACB_with_ids_1_31.vcf --make-bed --out plink_ACB_with_ids_1_31 --chr-set 27
plink --file ACB_with_ids.plink --indep-pairwise 20000 500 0.2 --out pruned --chr-set 27
plink --file ACB_with_ids.plink --extract pruned.prune.in --chr-set 27 --make-bed --out final_pruned_data
vcftools --vcf ACB_with_ids.vcf --snps pruned.prune.in --recode --recode-INFO-all --out filter_ACB_with_ids
PopLDdecay -InVCF /home/pyraloidea/LD_PR_CI_analysis/programs/filter_ACB_with_ids.recode.vcf -OutStat /home/pyraloidea/PopLDdecay/ACBLDdecay
perl Plot_OnePop.pl -inFile /home/pyraloidea/PopLDdecay/ACBLDdecay.stat.gz --output /home/pyraloidea/PopLDdecay/ACB -bin1 10 -bin2 100
java -Xmx27307m -jar beagle.18May20.d20.jar gt=ACB_with_id.vcf nthreads=8 out=pahse_ACB_with_id
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31;
> do
> vcftools --gzvcf pahse_ACB_with_ids.vcf.gz --recode --recode-INFO-all --chr ${i} --out pahse_ACB_with_ids.chr${i};
> done
for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29;
do 
vcftools --vcf JBC.chr${k}.recode.vcf --plink --out chr${k}.MT;
done
for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; 
do 
awk 'BEGIN{OFS=" "} {print 1,".",$4/1000000,$4}' chr${k}.MT.map > chr${k}.MT.map.distance;
done
for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; 
do 
/home/pyraloidea/xuanzexinhao/selscan-2.0.3/bin/linux/selscan-2.0.3 --ihs --vcf pahse_ACB_with_ids.chr${k}.recode.vcf --map chr${k}.MT.map.distance --out chr${k}.iHS;
done
for k in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31;
do 
awk  '{print '${k}',$2,$3,$4,$5,$6}' chr${k}.iHS.ihs.out > Chr${k}.ihs.out;
sed -i 's/ /\t/g' Chr${k}.ihs.out;      
done
for k in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 30 31;
do 
/home/software/selscan/bin/linux/norm --ihs --files Chr${k}.ihs.out --bp-win --winsize 50000;
done
for k in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 30 31;
do 
awk '{print '${k}',$1,$2,$4}' Chr${k}.ihs.out.100bins.norm.50kb.windows > Chr${k}.chart.ihs.out.50kb.windows;
cat ./*.chart.ihs.out.50kb.windows > all.ihs.out.50kb.windows;
done
sort -k 4n,4  all.ihs.out.50kb.windows > all.ihs.out.50kb.windows.sort
library(rehh)
hh <- data2haplohh(hap_file = "pahse_ACB_with_ids.chr26.recode.vcf",polarize_vcf = FALSE, vcf_reader = "data.table")
chicke_iHS = scan_hh(hh, polarized = FALSE)
iHS = ihh2ihs(chicke_iHS, freqbin = 1)
write.table(iHS, file = "28iHS.txt", sep = " ", quote = FALSE)
library(FastEPRR)
FastEPRR_VCF_step1(vcfFilePath="/home/pyraloidea/xuanzexinhao/pahse_ACB_with_ids.chr30.recode.vcf", winLength="50",srcOutputFilePath="/home/pyraloidea/xuanzexinhao/30_step1/chr30_step1")
FastEPRR_VCF_step2(srcFolderPath="/home/pyraloidea/xuanzexinhao/30_step1/",jobNumber=1,currJob=1,DXOutputFolderPath="/home/pyraloidea/xuanzexinhao/30_step2/")
FastEPRR_VCF_step3(srcFolderPath="/home/pyraloidea/xuanzexinhao/30_step1/", DXFolderPath="/home/pyraloidea/xuanzexinhao/30_step2/", finalOutputFolderPath="/home/pyraloidea/xuanzexinhao/30_step3")

##########################################################################################################################################################
# Part 5  Linkage disequilibrium, population recombination rates and select sweep ananlysis for allopatric populations
##########################################################################################################################################################
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; 
do  
vcftools --vcf /home/pyraloidea/4adaptive_analysis/north_vcf/north_with_ids.vcf --recode --recode-INFO-all --chr ${i} --out north_with_ids.chr${i}; 
done
sed -i 's/\t/,/g' chr_grid.txt
for line in `cat chr_grid.txt`;
do
chr=${line%%,*};
grid=${line##*,};
/home/pyraloidea/zhongqunxuanze/SweeD/sweed/SweeD -name ${chr}.50kb -input north_with_ids.chr1.recode.vcf -grid ${grid};
done
for k in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31; do 
awk '{print '${k}',$1,$2,$3,$4,$5}' SweeD_Report.${k}.50kb > SweeD_Report.chr.${k}.50kb; 
#sed -i 's/ /\t/g' SweeD_Report.chr.${k}.50kb; 
sed -i '1,3d' SweeD_Report.chr.*.50kb; 
cat ./SweeD_Report.chr.*.50kb > all.CLR.50k.txt; 
done

##########################################################################################################################################################
# Part 6 Identifying the numbers of adaptive introgressed SNPs and genes associated with different environmental factors and their functions.
##########################################################################################################################################################
vcftools --vcf north_with_ids.vcf --snps outliersnps_id.txt --recode --out outliersnps
/home/pyraloidea/ANNOVAR/annovar/convert2annovar.pl -format vcf4old outliersnps.recode.vcf > northoutliersnps.annovar
vcftools --vcf ACB.vcf  --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr 23 --chr 24 --chr 25 --chr 26 --chr 27 --chr 28 --chr 29 --chr 30 --chr 31 --recode --recode-INFO-all --stdout | gzip -c > north_chr1_31.vcf.gz
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/ANNOVAR$ gffread /home/pyraloidea/LBSNPanalysis/HF.gff3 -T -o HF.gtf(note:envs)
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/ANNOVAR$ conda install -c bioconda ucsc-gtfToGenePred
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/ANNOVAR$ gtfToGenePred -genePredExt HF.gtf HF_refGene.txt
/home/pyraloidea/ANNOVAR/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile /home/pyraloidea/LBSNPanalysis/HFgenome.fa --out HF_refGeneMrna.fa HF_refGene.txt
/home/pyraloidea/ANNOVAR/annovar/convert2annovar.pl -format vcf4old outliersnps.recode.vcf > northoutliersnps.annovar
/home/pyraloidea/ANNOVAR/annovar/table_annovar.pl /home/pyraloidea/4adaptive_analysis/north_vcf/north_chr1_31.vcf.gz /home/pyraloidea/ANNOVAR -buildver HF --outfile ./snp 
-remove -protocol refGene -operation g -nastring . -vcfinput
wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ mkdir 1-SnpEff
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ cp /home/pyraloidea/snpEFF/snpEff/snpEff.config /home/pyraloidea/snpEFF/snpEff_1/
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ echo "# my own genome, version Ostrinia_furnacalis
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ Ostrinia_furnacalis.genome:Ostrinia_furnacalis" >> snpEff.config
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ /home/pyraloidea/LBSNPanalysis/HFgenome.fa /home/pyraloidea/snpEFF/snpEff_1/data/Ostrinia_furnacalis/sequences.fa
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ /home/pyraloidea/LBSNPanalysis/HF.gff3 /home/pyraloidea/snpEFF/snpEff_1/data/Ostrinia_furnacalis/genes.gff
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ java -jar /home/daichunyan/software/snpEff/snpEff.jar build -c ./snpEff.config -gff3 -v Brassica_napus.ZS11.v0 -d -noCheckCds -noCheckProtein
java -jar ~/software/snpEff/snpEff.jar -c ~/software/Arabidopsis-SnpEff/snpEff.config -ud 1000 -csvStats north.csv -htmlStats north.html -o vcf Ostrinia_furnacalis north_chr1_31.vcf > /home/pyraloidea/snpEFF/north/north.anno.vcf
(jiyinjianshen) pyraloidea@shpc-4067-instance-b6HRjrZJ:~/snpEFF/snpEff_1$ java -jar /home/pyraloidea/snpEFF/snpEff/snpEff.jar -c ./snpEff.config -ud 1000 -csvStats north.csv -htmlStats north.html -o vcf Ostrinia_furnacalis north_chr1_31.vcf > north_anno.csv

#####################################################
# Part 7 identification of inversion by Delly v0.8.7
#####################################################
bwa index HFgenome.fa
samtools faidx HFgenome.fa
index=/home/pyraloidea/LBSNPanalysis/chr_inv/HFgenome.fa
data=/home/pyraloidea/LBSNPanalysis/chr_inv/1cleandata/
out=/home/pyraloidea/LBSNPanalysis/chr_inv/2align/
for i in $(cat sample.list); 
do 
bwa mem -t 20 -R "@RG\tID:${i}\tSM:${i}\tLB:reseq\tPL:Illumina" 
${index} 
${data}${i}"_1.fq.gz" 
${data}${i}"_2.fq.gz" 
-o ${out}/bwa/${i}".sam" 2>${out}/bwa/${i}".bwa.align.log"; 
done
for i in $(cat sample.list);
do
samtools view
-@ 20
-bS ${out}/bwa/${i}".sam"
-o ${out}/bwa/${i}".bam" 
done
for i in $(cat sample.list); 
do 
samtools sort -o ${out}/sort/${i}".sort.bam" ${out}/bwa/${i}".bam"; 
done
for i in $(cat sample.list);
do
samtools sort -o ${out}/sort/${i}".sort.bam" ${out}/bwa/${i}".bam";
done
for i in $(cat sample.list);
do
/home/pyraloidea/LBSNPanalysis/chr_inv/gatk-4.6.0.0/gatk MarkDuplicates REMOVE_DUPLICATES=false \
INPUT=${out}/sort/${i}.sort.bam \
OUTPUT=${out}/mark/${i}.sort.marked.bam \
METRICS_FILE=${out}/mark/${i}.markdup_metrics_list
done
for i in $(cat sample.list)
do
sambamba view -t 50 -h -f bam \
-F "mapping_quality >= 10 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${out}/mark/$i".sort.marked.bam" \
-o ${out}/filter/$i".sort.mark.filter.bam"
done
/home/pyraloidea/delly_test/delly_v1.2.6_linux_x86_64bit call -t INV -o APT1_delly.bcf -g HFgenome.fa /home/pyraloidea/LBSNPanalysis/chr_inv/2align/filter/APT1.sort.mark.filter.bam
bcftools view APT1_delly.bcf -O v -o APT1_delly.vcf
less APT1_delly.vcf|head -336 >APT1_delly_filter.vcf
less APT1_delly.vcf|awk '{if($7=="PASS" && $6 >=100)print}'>>APT1_delly_filter.vcf
bgzip APT1_delly_filter.vcf
bcftools index APT1_delly_filter.vcf.gz
bcftools isec -n+2 -c all -p merge APT1_delly_filter.vcf.gz APT2_delly_filter.vcf.gz APT3_delly_filter.vcf.gz APT4_delly_filter.vcf.gz APT5_delly_filter.vcf.gz
bcftools merge --force-samples -o APT_merge.vcf.gz -O z APT1_delly_filter.vcf.gz APT2_delly_filter.vcf.gz APT3_delly_filter.vcf.gz APT4_delly_filter.vcf.gz APT5_delly_filter.vcf.gz
bcftools index APT_merge.vcf.gz
bcftools view -T merge/sites.txt APT_merge.vcf.gz -O z -o APT.select.vcf.gz