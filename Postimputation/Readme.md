
# Extract SNPs with MAF > 0.1% and r2 > 0.6
## Step1a
```{R}
library(data.table)
setwd("/mnt/beegfs/home4/arc/vw260/SPARK/SPARK_v3/Postimputation_genotyping/omega/")
c = NULL
for (i in 10:22){
  a = read.table(paste0("chr", i, ".info.gz"), header = T)
  a$Rsq = as.numeric(as.character(a$Rsq))
  b = subset(a, Rsq > 0.6)
  b = subset(b, ALT_Frq > 0.005 & ALT_Frq < 0.995)
  write.table(b[,1], file = paste0("chr", i, "extract.txt"), row.names = F, col.names = T, quote = F)
  c = rbind(c, b)
}

a = read.table("chrX.info.gz", header = T)
a$Rsq = as.numeric(as.character(a$Rsq))
b = subset(a, Rsq > 0.6)
b = subset(b, ALT_Frq > 0.005 & ALT_Frq < 0.995)
write.table(b[,1], file = "chrXextract.txt", row.names = F, col.names = T, quote = F)
c = rbind(c, b)

write.table(c, file = "AllextractSNPs.txt", row.names = F, col.names = T, quote = F)

```

## Step 1b
```{bash}
for i in {1..22}; do ./plink --vcf ./Postimputation_genotyping/omega/chr${i}.dose.vcf.gz --make-bed --out ./Postimputation_genotyping/omega/omega${i}  --extract ./Postimputation_genotyping/omega/chr${i}extract.txt --const-fid 0; done
```

#Step 2: Update SNPids and FAMnames, and recode to plink ped/map files for liftover
##Step 2a

try: https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/

```{bash}
wget https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
zgrep -v "^##" common_all_20180418.vcf.gz | cut -f1-3 > fileforrecoding.txt
awk '{print $1":"$2"\t"$3}' < fileforrecoding.txt > plinkrecodingfile.txt
```

## Step 2b

```{R}
setwd("/mnt/beegfs/home4/arc/vw260/SPARK/SPARK_v3/Postimputation_genotyping/omega")

library(data.table)
library(tidyr)


bim22 = fread("alpha22.bim")
bim22$CHROM_POS = paste(bim22$V1, bim22$V4, sep = ":")
bim22_update = bim22[, c("CHROM_POS", "V2")]
update = bim22_update

for (i in 1:21){
  bim = fread(paste0("alpha", i, ".bim"))
  bim$CHROM_POS = paste(bim$V1, bim$V4, sep = ":")
  bim_update = bim[, c("CHROM_POS", "V2")]
  update = rbind(update, bim_update)
}

plink_recoding = fread("/mnt/beegfs/home4/arc/vw260/SPARK/SPARK_v3/plinkrecodingfile.txt")
plink_recoding_file2 = merge(plink_recoding, update, by.x = "#CHROM:POS", by.y = "CHROM_POS")
plink_recoding_file = plink_recoding_file2[!duplicated(plink_recoding_file2$V2), ]

write.table(plink_recoding_file[,c("V2", "ID")], file = "~/SPARK/SPARK_v3/Postimputation_genotyping/alpha/plinkrecodingfile_final.txt", col.names = T, row.names = F, quote = F)


# Generate updated filenames
fileimputed = fread("omega1.fam")
fileimputed$oldFID = fileimputed$V1
fileimputed$oldIID = fileimputed$V2


fileimputed = fileimputed %>% separate(V2, into = c('FID', 'IID'), sep = 10)
fileimputed$FID = substr(fileimputed$FID,1,nchar(fileimputed$FID)-1)

write.table(fileimputed[,c("oldFID", "oldIID", "FID", "IID")], file = "~/SPARK/SPARK_v3/Postimputation_genotyping/omega/updatefamnames.txt", row.names = F, col.names = F, quote = F)
```

##Step 2c

```{bash}
for i in {1..22}; do ./plink --bfile ./Postimputation_genotyping/omega/omega${i}  --update-name ./Postimputation_genotyping/alpha/plinkrecodingfile_final.txt  --update-ids ./Postimputation_genotyping/omega/updatefamnames.txt --recode --out ./Postimputation_genotyping/omega/omega${i}_v4_try --threads 15; done
```



# Step 3: liftover to hg19 to keep it consistent with existing GWAS

```{bash}
for i in {1..22}; do python2 ~/ABCD/ABCDgenotype/Genotype_postimputation/TOPMED/Plink_files/liftOverPlink.py -m ./Postimputation_genotyping/omega/omega${i}_v4_try.map -p ./Postimputation_genotyping/omega/omega${i}_v4_try.ped  -o ./Postimputation_genotyping/omega/omega_chr${i}_hg19 -e ~/ABCD/ABCDgenotype/Genotype_postimputation/TOPMED/Plink_files/liftOver -c ~/ABCD/ABCDgenotype/Genotype_postimputation/TOPMED/Plink_files/hg38ToHg19.over.chain.gz; done

for i in {1..22}; do ./plink --file ./Postimputation_genotyping/omega/omega_chr${i}_hg19 --make-bed --out ~/SPARK/SPARK_v3/Postimputation_genotyping/omega/omega_chr${i}_hg19 --allow-extra-chr; done
```


# Step 4: post liftover cleaning
## Step 4a

```{R}
library(data.table)

setwd("/mnt/beegfs/home4/arc/vw260/SPARK/SPARK_v3/Postimputation_genotyping/omega")
for (i in 1:22){
  data1 = fread(paste0("alpha_chr",i,"_hg19.bim"))
  data2 = subset(data1, V1 == i) # keep only those with the right chromosome
  data3 = data2[data2$V2 %like% "rs", ] # keep only rsid
  data4 = data3[!duplicated(data3[, 2])]
  write.table(data4[,"V2"], file = paste0("alpha_snpstoextract_",i,"_hg19.txt"), row.names = F, col.names = F, quote = F)
}
```


## Step 4b
```{bash}
for i in {1..22}; do ./plink --bfile ./Postimputation_genotyping/alpha/alpha_chr${i}_hg19 --make-bed --out ./Postimputation_genotyping/alpha/alpha_chr${i}_hg19_cleaned  --extract ./Postimputation_genotyping/alpha/alpha_snpstoextract_${i}_hg19.txt --allow-extra-chr; done
```



