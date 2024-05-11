#程序运行前需要把有意义的代谢物源文件下载到这个目录下
#本程序不需要修改
rm(list = ls()) #清除缓存
library(TwoSampleMR)
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/06.immMetMR")     #设置工作目录

#免疫细胞的暴露数据文件
exposureFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/02.immMR/1.imm.exposure.F10.csv"  
#免疫细胞过滤的结果文件
sigCellFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/03.immReverse/immune-disease.REVfilter.csv"       


#读取暴露数据输入文件
rt=read.csv(exposureFile, header=T, sep=",", check.names=F)

#读取免疫细胞过滤的结果文件
sigCell=read.csv(sigCellFile, header=T, sep=",", check.names=F)

#读取目录下所有代谢物的原始文件
files=dir()                                 #获取目录下所有文件
gzFiles=grep("\\.gz$", files, value=T)      #提取.gz结尾的文件

#对免疫细胞进行循环
for(cell in unique(sigCell$id.exposure)){
	#提取这个免疫细胞的暴露数据
	k=gsub("\\%|\\/|\\:", "_", cell)
	singleExposureFile=paste0(k, ".exposure.csv")
	exposure_set=rt[rt$id.exposure==cell,]
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
	
	#对代谢物的文件进行循环
	for(outcomeFile in gzFiles){
		outcomeName=gsub("_buildGRCh38.tsv.gz", "", outcomeFile)
		i=paste0("beta1.", cell, "_", outcomeName)
		#读取结局数据
		outcome_data=read_outcome_data(snps=exposure_dat$SNP,
			             filename=outcomeFile, sep = "\t",
			             snp_col = "variant_id",
			             beta_col = "beta",
			             se_col = "standard_error",
			             effect_allele_col = "effect_allele",
			             other_allele_col = "other_allele",
			             pval_col = "p_value",
			             eaf_col = "effect_allele_frequency")
		
		#将暴露数据和结局数据合并
		outcome_data$outcome=outcomeName
		dat=harmonise_data(exposure_dat, outcome_data)
		
		#MR-PRESSO异常值检测(偏倚的SNP)
		#presso=run_mr_presso(dat)
		#write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
		#write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
		
		#孟德尔随机化分析
		mrResult=mr(dat)
		mrTab=generate_odds_ratios(mrResult)
		
		#筛选IVW方法pvalue小于0.05的结果
		if(mrResult$pval[3]<0.05){
			#筛选五种方法OR方向一致的结果
			if(sum(mrTab$or>1)==nrow(mrTab) | sum(mrTab$or<1)==nrow(mrTab)){
				#输出用于孟德尔随机化的工具变量
				outTab=dat[dat$mr_keep=="TRUE",]
				write.csv(outTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
				
				#输出孟德尔随机化分析的结果
				write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)
				
				#异质性检验
				heterTab=mr_heterogeneity(dat)
				write.csv(heterTab, file=paste0(i, ".table.heterogeneity.csv"), row.names=F)
				
				#多效性检验
				pleioTab=mr_pleiotropy_test(dat)
				write.csv(pleioTab, file=paste0(i, ".table.pleiotropy.csv"), row.names=F)
				
				#绘制散点图
				pdf(file=paste0(i, ".scatter_plot.pdf"), width=7, height=6.5)
				p1=mr_scatter_plot(mrResult, dat)
				print(p1)
				dev.off()
				
				#森林图
				res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
				pdf(file=paste0(i, ".forest.pdf"), width=6.5, height=5)
				p2=mr_forest_plot(res_single)
				print(p2)
				dev.off()
				
				#漏斗图
				pdf(file=paste0(i, ".funnel_plot.pdf"), width=6.5, height=6)
				p3=mr_funnel_plot(singlesnp_results = res_single)
				print(p3)
				dev.off()
				
				#留一法敏感性分析
				pdf(file=paste0(i, ".leaveoneout.pdf"), width=6.5, height=5)
				p4=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
				print(p4)
				dev.off()
			}
		}
	}
}




