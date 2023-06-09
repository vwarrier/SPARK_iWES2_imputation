# Pre-imputation

## Step 1: Delineate to European ancestry

```{R}
setwd("/mnt/b2/home4/arc/vw260/SPARK/SPARK_v3/")

ancestry = fread("SPARK.iWES_v2.ancestry.2023_01.tsv") #106781

head(ancestry)

EUR = subset(ancestry, superclass == "EUR")
fam = fread("SPARK.iWES_v2.genotyping_by_sequencing.fam")
fam2 = fam[fam$V2 %in% EUR$spid,] ##51869 genotyping ###26046 genotypebysequencing

write.table(fam2[,c("V1", "V2")], file = "EUR_keep_sequencing.txt", row.names = F, col.names = F, quote = F)

```


## Step 2: Pre-imputation quality control

```{bash}
./plink --bfile SPARK.iWES_v2.genotyping_by_sequencing --keep EUR_keep_sequencing.txt --make-bed --geno 0.05 --out ~/SPARK/SPARK_v3/Preimputation_sequencing/QC1output --threads 15
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_sequencing/QC1output --make-bed --me 0.05 0.1 --mind 0.02 --out ~/SPARK/SPARK_v3/Preimputation_sequencing/QC2output --threads 15
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_sequencing/QC2output --make-bed --out ~/SPARK/SPARK_v3/Preimputation_sequencing/QC3output --remove ~/SPARK/SPARK_v3/Preimputation_sequencing/QC2output.irem --threads 15
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_sequencing/QC3output --hwe 0.000001 --make-bed --threads 15 --out ~/SPARK/SPARK_v3/Preimputation_sequencing/QC4output 
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_sequencing/QC4output --het --out ~/SPARK/SPARK_v3/Preimputation_sequencing/check_het --threads 15
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_sequencing/QC4output --check-sex --out ~/SPARK/SPARK_v3/Preimputation_sequencing/check_sex --threads 15
```

```{R}
R
library(data.table)
het = fread("Preimputation_sequencing/check_het.het")
het$HET = (het$`N(NM)` - het$`O(HOM)`)/het$`N(NM)` #create heterozygosity stats
mean = mean(het$HET)
sd = sd(het$HET)
het$Z = (het$HET - mean)/sd #create Z scores of heterozygosity
hetoutlier = subset(het, abs(Z) > 3)
failedsample = hetoutlier[,c(1:2)]

sex_check = fread("Preimputation_sequencing/check_sex.sexcheck")
problem = subset(sex_check, STATUS == "PROBLEM")
problem2 = problem[,c(1:2)]

failedsample2 = rbind(failedsample, problem2)
failedsample2 = unique(failedsample2)

write.table(failedsample2, file = "Preimputation_sequencing/failedsample_het.txt", row.names = F, col.names = T, quote = F)
q()
n
```
```{bash}
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_sequencing/QC4output --make-bed --out ~/SPARK/SPARK_v3/Preimputation_sequencing/QC5output --remove ~/SPARK/SPARK_v3/Preimputation_sequencing/failedsample_het.txt --threads 15
```


## Step 3: Update for imputation 

```{R}
library(data.table)
bim = fread("Preimputation_sequencing/QC5output.bim")
bim$update <- gsub("GSA-", "", bim$V2)
write.table(bim[,c("V2", "update")], file = "Preimputation_genotype/updatenames.txt", row.names = F, col.names = F, quote = F)
q()
n
```
```{bash}
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC5output --update-name ~/SPARK/SPARK_v3/Preimputation_genotype/updatenames.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output 
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output --freq --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6outputQC6outputfreq

perl HRC-1000G-check-bim.pl -b  ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output.bim -f  ~/SPARK/SPARK_v3/Preimputation_genotype/QC6outputQC6outputfreq.frq -h -r PASS.Variantsbravo-dbsnp-all.tab

./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output --exclude Exclude-QC6output-HRC.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP1 --threads 20
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP1 --update-map Chromosome-QC6output-HRC.txt --update-chr --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP2 --threads 20
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP2 --update-map Position-QC6output-HRC.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP3 --threads 20
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP3 --flip Strand-Flip-QC6output-HRC.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP4 --threads 20
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP4 --reference-allele Force-Allele1-QC6output-HRC.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated --threads 20

rm ~/SPARK/SPARK_v3/Preimputation_genotype/TEMP*
```

## Step 4: Split into ubgroups and convert files for imputation

In R to split it into two subsets because TOPMED accepts only 25K max

```{R}
fam1 = fread("~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated.fam")
omega = sample_n(fam1, 23585, replace = FALSE)
alpha = fam1[!fam1$V2 %in% omega$V2,]

write.table(alpha, file = "~/SPARK/SPARK_v3/Preimputation_genotype/alpha_keep.txt", row.names = F, col.names = T, quote = F)
write.table(omega, file = "~/SPARK/SPARK_v3/Preimputation_genotype/omega_keep.txt", row.names = F, col.names = T, quote = F)
```
```{bash}
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated --keep ~/SPARK/SPARK_v3/Preimputation_genotype/alpha_keep.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated-alpha --threads 20 --output-chr chr26
./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated --keep ~/SPARK/SPARK_v3/Preimputation_genotype/omega_keep.txt --make-bed --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated-omega --threads 20 --output-chr chr26

for i in {1..23}; do ./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated-alpha --reference-allele Force-Allele1-QC6output-HRC.txt --chr ${i} --recode-vcf --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_alpha_chr${i} --output-chr chr26; done
for i in {1..23}; do vcf-sort ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_alpha_chr${i}.vcf | bgzip -c > ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_alpha_chr${i}.vcf.gz; done
for i in {1..23}; do rm ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_alpha_chr${i}.vcf | rm ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_alpha_chr${i}.log; done

for i in {1..23}; do ./plink --bfile ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output-updated-omega --reference-allele Force-Allele1-QC6output-HRC.txt --chr ${i} --recode-vcf --out ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_omega_chr${i} --output-chr chr26; done
for i in {1..23}; do vcf-sort ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_omega_chr${i}.vcf | bgzip -c > ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_omega_chr${i}.vcf.gz; done
for i in {1..23}; do rm ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_omega_chr${i}.vcf | rm ~/SPARK/SPARK_v3/Preimputation_genotype/QC6output_file_omega_chr${i}.log; done
```
