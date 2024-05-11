#该程序不需要运行
library(TwoSampleMR)

inputFile="GCST90199621_buildGRCh38.tsv.gz"     #输入文件
outFile="GCST90199621.exposure_data.csv"        #输出的结果文件
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/04.metData")     #设置工作目录

#读取输入文件, 并对输入文件进行格式转换
data<-read_exposure_data(filename=inputFile,
                         sep = "\t",
                         snp_col = "variant_id",
                         beta_col = "beta",
                         se_col = "standard_error",
                         effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele",
                         eaf_col = "effect_allele_frequency",
                         pval_col = "p_value",
                         chr_col="chromosome",
                         pos_col = "base_pair_location",
                         clump = F)

#根据pvalue<1e-05对数据进行过滤
exposure_dat<-subset(data, pval.exposure<1e-05)

#去除连锁不平衡的SNP
exposure_dat_clumped <- clump_data(exposure_dat, clump_kb=10000, clump_r2=0.001)
write.csv(exposure_dat_clumped, file=outFile, row.names=F)




