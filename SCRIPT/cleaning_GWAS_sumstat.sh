# This script is related with the Vitamin D project
# To run GSMR on several traits, I've collected many GWAS summary sumstat in the litterature.
# To comply with GCTA format, I've cleaned these GAWAS sumstat: this is what this script does:


# A few rules
# - To go from Odd ratio to Beta, just use a log in awk







# ------ MULTIPLE SCLEROSIS ---------- #

# Data received per email, uq adress, 13 Jun 2018.
# Provides Odd rations.
# Alleles are ordered in a weird way.

# The Odd ratios does not make sense.
data <- read.table("Immunochip_FinalResults_LimitedDiscovery.txt", head=T, sep=" ")
library(tidyverse)
data %>% arrange(desc(OR)) %>% head(20)
# let's try to prune on allele frequency?
data %>% filter(Risk_Allele_Freq>0.05 & Risk_Allele_Freq<0.95) %>% arrange(desc(OR)) %>% head(20)
#OK this makes more sense.
# A few pvalues = 0 and must be converted to 10^-323?

#Odd ratio must be converted to Beta
echo "SNP A1 A2 freq b se p n" > MultipleSclerosis.ma
zcat Immunochip_FinalResults_LimitedDiscovery.txt  |  awk '{ if(NR>1 && $7>0.05 && $7<0.95){ print $4, $5, $6, $7, log($10), $9, 38589}}' >> MultipleSclerosis.ma






# ------ PARKINSON ---------- #
# Reformat Parkinson's Disease data from Constanza. I need to add the allele frequency of each SNP
ls /shares/compbio/Group-Wray/YanHoltz/DATA/ALLELE_FREQUENCY/*
cd /shares/compbio/Group-Wray/YanHoltz/DATA/GWAS/GWAS_SUMSTAT 
echo "SNP A1 A2 freq b se p n" > tmp
join <(cat META_ANALYSIS_10K23_beta_se_correctGC1_pdgene_sharing_280317.tbl | grep -v "^SNP" | sort -k 1,1) <(awk '{print $2, $6}' /shares/compbio/Group-Wray/YanHoltz/DATA/ALLELE_FREQUENCY/* | sort -k 1,1) | awk '{ print $1,$2,$3,$9,$4,$5,$6,100000}' >> tmp
mv tmp META_ANALYSIS_10K23_beta_se_correctGC1_pdgene_sharing_280317.ma







# ------ TANNING ---------- #

cd /shares/compbio/Group-Wray/YanHoltz/DATA/GWAS/GWAS_SUMSTAT/ORIGINAL_FILES

#frequency = /shares/compbio/Group-Wray/YanHoltz/DATA/ALLELE_FREQUENCY/HRC/*freq

# Now I can build the .ma file for GSMR
echo "SNP A1 A2 freq b se p n" > GWAS_tanning_UKBB.ma
join <(awk '{print $2, $5, $6}' /gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EURu_impHRC/ukbEURu_imp_chr*.bim | sort -k 1,1) <(awk '{print $2, $6}' /shares/compbio/Group-Wray/YanHoltz/DATA/ALLELE_FREQUENCY/HRC/*frq | sort) | join - <(zcat GWAS_tanning_UKBB.txt.gz | awk '{
if($12==0 && NR>1){ print $2, log($7), $8, 1.0e-323, $6 } ; if($12!=0 && NR>1){ print $2, log($7), $8, $12, $6 }}' | sort -k 1,1) >> GWAS_tanning_UKBB.ma





# ------ PROSTATE CANCER ---------- #

wget http://practical.icr.ac.uk/blog/wp-content/uploads/uploadedfiles/oncoarray/summaryData/results_onco_sample_euro_bycountry_nooverlap_imp_chr22_varname_se.fixed.txt.zip






# ------ IPAQ ---------- #

# Now I can build the .ma file for GSMR
cd /shares/compbio/Group-Wray/YanHoltz/DATA/GWAS/GWAS_SUMSTAT

echo "SNP A1 A2 freq b se p n" > ukbEUR_MAF1e-4_IPAQG.bolt.ma
join <(awk '{print $2, $5, $6}' /gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EUR_impHRC/ukbEUR_imp_chr*_v2_imp_QC_HRC.bim | head | sort -k 1,1) <(awk '{print $2, $6}' /shares/compbio/Group-Wray/YanHoltz/DATA/ALLELE_FREQUENCY/all_frequency*.frq | sort) | join - <(cat /shares/compbio/PCTG/a.xue/for_Yan/ukbEUR_MAF1e-4_IPAQG.bolt | awk '{if(NR>1) print $2, $7, $7/$8, $9, $6 }' | sort -k 1,1) >> ukbEUR_MAF1e-4_IPAQG.bolt.ma

































##GTEX Liver
2 causal associations have been detected:
  ```{r}
gtex <- read.table("0_DATA/smr_VitaminDXiaEtAl_GTEXLiver.smr", header=T)

gtex %>% filter(p_SMR < (0.05 / nrow(gtex)) )
```

Once more the 2 regions of the chromosome 11 are highlighted: around 15 Mb and around 71 Mb. This time, it's the genes SPON1 and RP11-66L16.2 that are found has having a causal effect on VitaminD. The gene SPON1 is located a bit before GWAS highlighted region (14.1 Mb vs 14.9 Mb). 

Note: the gene RP11-660L16.2 = ENSG00000254682 = AP002387.1

```{bash, eval=FALSE}
# Good folder
cd /shares/compbio/Group-Wray/YanHoltz/VITAMIND_XIA_ET_AL/5_SMR/ON_GTEX_LIVER

# Load the position of genes:
wget https://www.cog-genomics.org/static/bin/plink/glist-hg18

# Send smr plot for this loci (first of chromosome 11)
tmp_command="smr_Linux --bfile /gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EURu_HM3/ukbEURu_imp_chr11_v2_HM3_QC --gwas-summary /shares/compbio/Group-Wray/YanHoltz/VITAMIND_XIA_ET_AL/1_GWAS/GWAS_vitaminD_XiaEtAL.ma --beqtl-summary /gpfs/gpfs01/polaris/Q0286/GTExV7/Summary/besd/Liver --out myplot --plot --probe ENSG00000152268.8	 --probe-wind 500 --gene-list /shares/compbio/Group-Wray/YanHoltz/VITAMIND_XIA_ET_AL/5_SMR/ON_EQTLGEN_32K/glist-hg18"
qsubshcom "$tmp_command" 1 30G smr_plot_vitaminD_loc1 10:00:00 ""

# Send smr plot for this loci (Second of chromosome 11)
tmp_command="smr_Linux --bfile /gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EURu_HM3/ukbEURu_imp_chr11_v2_HM3_QC --gwas-summary /shares/compbio/Group-Wray/YanHoltz/VITAMIND_XIA_ET_AL/1_GWAS/GWAS_vitaminD_XiaEtAL.ma --beqtl-summary /gpfs/gpfs01/polaris/Q0286/GTExV7/Summary/besd/Liver --out myplot --plot --probe ENSG00000254682.1	 --probe-wind 500 --gene-list /shares/compbio/Group-Wray/YanHoltz/VITAMIND_XIA_ET_AL/5_SMR/ON_EQTLGEN_32K/glist-hg18"
qsubshcom "$tmp_command" 1 30G smr_plot_vitaminD_loc1 10:00:00 ""

# transfer locally
cd /Users/y.holtz/Dropbox/QBI/4_UK_BIOBANK_GWAS_PROJECT/VitaminD-GWAS/0_DATA
scp  y.holtz@delta.imb.uq.edu.au:/shares/compbio/Group-Wray/YanHoltz/VITAMIND_XIA_ET_AL/5_SMR/ON_GTEX_LIVER/plot/*  .
```

```{r, warning=FALSE, message=FALSE, fig.align="center", fig.width=12, fig.height=9}
# Make the plot
source("SCRIPT/plot_SMR.r") 
# Read the data file in R:
SMRData = ReadSMRData("0_DATA/myplot.ENSG00000152268.8.txt")
# Plot the SMR results in a genomic region centred around a probe:
SMRLocusPlot(data=SMRData, smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)
```





```{r, warning=FALSE, message=FALSE, fig.align="center", fig.width=12, fig.height=9}
# Make the plot
source("SCRIPT/plot_SMR.r") 
# Read the data file in R:
SMRData = ReadSMRData("0_DATA/myplot.ENSG00000254682.1.txt")
# Plot the SMR results in a genomic region centred around a probe:
SMRLocusPlot(data=SMRData, smr_thresh=8.4e-6, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)
```























########## 22 diseases GWAS summary QC ##########

## Read the file
PATH="/shares/compbio/Group-Wray/YanHoltz/DATA/GWAS/GWAS_SUMSTAT/TMP/"

library(data.table)

trait=read.table(paste(PATH,"ukbEUR_MAF1e-4_IPAQG.bolt",sep=""),header=T)

## Get sample size N
phen=fread("/shares/compbio/PCTG/a.xue/for_Yan/IPAQG_zscore_adj_by_sex_age.phen",header=T)
phen=as.data.frame(phen)

co=fread("/gpfs/gpfs01/polaris/Q0286/uqaxue/phen/covar_sex_age_10PCs.txt",header=T)
co=as.data.frame(co)

ind=Reduce(intersect,list(phen$IID,co$IID))
phen=phen[match(ind,phen$IID),]

N=sum(table(phen$IPAQG))
print(paste("The sample size is ",N),sep="")

## How many columns of P-value, keep the later one if there are two
if(ncol(trait)==12){trait=trait[,-12]}
trait=trait[,-4]

## Rename the colnames
colnames(trait)=c("SNP","CHR","BP","A1","A2","freq","N","b","se","P")

trait$N=N*(1-trait$N)
trait=trait[,c(2,1,3:10)]


## Clean
trait$N=as.integer(trait$N)

#### !!! Need to be careful if there is any NA or Inf in the data !!!! ####
print("Number of NAs is")
print(sum(rowSums(is.na(trait))>0))

print("Start to save summary statistics......")

## Save common SNPs with MAF >= 0.01
#write.table(trait[trait$freq>=0.01,],paste(PATH,"ukbEUR_trait_common.txt",sep=""),row.names=F,col.names=T,quote=F)

## Also save rare summary
#write.table(trait[trait$freq<0.01,],paste(PATH,"ukbEUR_trait_rare.txt",sep=""),row.names=F,col.names=T,quote=F)

## Save COJO & SMR format file

write.table(trait[trait$freq>=0.01,c("SNP","A1","A2","freq","b","se","P","N")],paste(PATH,"ukbEUR_trait_cojo.txt",sep=""),row.names=F,col.names=T,quote=F)
## Save LDSC format file
#hm3=read.table("/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/w_hm3.snplist",header=T)
#trait$Z=trait$b/trait$se
#snp=Reduce(intersect,list(hm3$SNP,trait$SNP))
#new=trait[match(snp,trait$SNP),]

#write.table(new[,c("SNP","A1","A2","Z","N","P")],paste(PATH,"trait_ldsc.txt",sep=""),row.names=F,col.names=T,quote=F)


## The End
print("trait QC completed!")








