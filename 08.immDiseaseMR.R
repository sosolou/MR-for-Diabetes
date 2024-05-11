#该程序只有一处修改，疾病名称，然后运行即可
rm(list = ls()) #清除缓存
library(TwoSampleMR)



# 获取目录下所有文件夹的名称
directory_path <- "C:/Users/yuanfang/Desktop/ImmuneMedMR/07.metDiseaseMR"
folder_names <- list.dirs(directory_path, full.names = FALSE)
folder_names <- folder_names[folder_names != ""]
folder_names1<- unique(substr(folder_names, 1, 18))

exposureFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/02.immMR/1.imm.exposure.F10.csv"        #暴露数据文件
outcomeFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/02.immMR/2.outcome.csv"     #结局数据文件

#读取暴露数据输入文件
rt=read.csv(exposureFile, header=T, sep=",", check.names=F)

for (exposureID in folder_names1) {
  setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/08.immDiseaseMR")
  dir.create(exposureID, showWarnings = FALSE)

  #修改1：疾病名称
outcomeName="Gestational diabetes"  

folder_path <- paste0("C:/Users/yuanfang/Desktop/ImmuneMedMR/08.immDiseaseMR/",exposureID)
setwd(folder_path)

#提取这个免疫细胞的暴露数据
i=paste0("betaAll.", exposureID)
singleExposureFile=paste0("1.",i, ".exposure.csv")
exposure_set=rt[rt$id.exposure==exposureID,]
write.csv(exposure_set, file=singleExposureFile, row.names=F)
	
#读取这个免疫细胞的暴露数据
exposure_dat=read_exposure_data(filename=singleExposureFile,
                                sep = ",",
                                snp_col = "SNP",
                                beta_col = "beta.exposure",
                                se_col = "se.exposure",
                                pval_col = "pval.exposure",
                                effect_allele_col="effect_allele.exposure",
                                other_allele_col = "other_allele.exposure",
                                eaf_col = "eaf.exposure",
                                phenotype_col = "exposure",
                                id_col = "id.exposure",
                                samplesize_col = "samplesize.exposure",
                                chr_col="chr.exposure", pos_col = "pos.exposure",
                                clump=FALSE)

#读取结局数据
outcome_data=read_outcome_data(snps=exposure_dat$SNP,
		             filename=outcomeFile, sep = ",",
		             snp_col = "SNP",
		             beta_col = "beta.outcome",
		             se_col = "se.outcome",
		             effect_allele_col = "effect_allele.outcome",
		             other_allele_col = "other_allele.outcome",
		             pval_col = "pval.outcome",
		             eaf_col = "eaf.outcome")
	
#将暴露数据和结局数据合并
outcome_data$outcome=outcomeName
dat=harmonise_data(exposure_dat, outcome_data)
	
#输出用于孟德尔随机化的工具变量
outTab=dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file=paste0("2.",i, ".SNP.csv"), row.names=F)
	
#MR-PRESSO异常值检测(偏倚的SNP)
presso=run_mr_presso(dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0("3.",i, ".MR-PRESSO_Global.csv"))
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0("4.",i, ".MR-PRESSO_Outlier.csv"))
	
#孟德尔随机化分析
mrResult=mr(dat)

#筛选IVW方法pvalue小于0.05的结果
if(mrResult$pval[3]<0.05){
  #筛选五种方法OR方向一致的结果
  if(sum(mrResult$b>0)==nrow(mrResult) | sum(mrResult$b<0)==nrow(mrResult)){	
#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file=paste0("5.",i, ".MRresult.csv"), row.names=F)
	
#异质性检验
heterTab=mr_heterogeneity(dat)
write.csv(heterTab, file=paste0("6.",i, ".heterogeneity.csv"), row.names=F)
	
#多效性检验
pleioTab=mr_pleiotropy_test(dat)
write.csv(pleioTab, file=paste0("7.",i, ".pleiotropy.csv"), row.names=F)
	
#绘制散点图
pdf(file=paste0("8.",i, ".scatter_plot.pdf"), width=7, height=6.5)
p1=mr_scatter_plot(mrResult, dat)
print(p1)
dev.off()

#森林图
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file=paste0("9.",i, ".forest.pdf"), width=6.5, height=5)
p2=mr_forest_plot(res_single)
print(p2)
dev.off()
	
#漏斗图
pdf(file=paste0("10.",i, ".funnel_plot.pdf"), width=6.5, height=6)
p3=mr_funnel_plot(singlesnp_results = res_single)
print(p3)
dev.off()
	
#留一法敏感性分析
pdf(file=paste0("11.",i, ".leaveoneout.pdf"), width=6.5, height=5)
p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
print(p4)
dev.off()

rm(folder_path,i,singleExposureFile,exposure_set,exposure_dat,outcome_data,dat,outTab,
   presso,mrResult,mrTab,heterTab,pleioTab,res_single)	

  }
  
}

}
