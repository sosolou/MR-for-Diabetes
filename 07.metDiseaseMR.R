#该程序只需修改1处，疾病的名称
#运行该程序前，需要把06.immMetMR文件夹中所有后缀为table.SNP复制到07.metDiseaseMR文件夹
library(TwoSampleMR)


# 设置工作目录为包含CSV文件的文件夹路径
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/07.metDiseaseMR")

# 获取文件夹中所有CSV文件的名称
file_list <- list.files(pattern = "\\.csv$")

# 遍历每个CSV文件名称，截取固定字符后创建新文件夹
for (file_name in file_list) {
  setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/07.metDiseaseMR")
  #新建文件夹
  new_folder_name <- substr(file_name, start = 7, stop = 37) 
  dir.create(new_folder_name, showWarnings = FALSE)
  #代谢物的ID
  exposureID=substr(file_name, start = 26, stop = 37)  #代谢物的ID   
  
  #修改1：结局中展示疾病的名字
  outcomeName="Gestational diabetes" 
  
  exposureFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/05.metMR/1.met.exposure.F10.csv"        #暴露数据文件
  outcomeFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/05.metMR/2.outcome.csv"  #结局数据文件
  
  
  #设置新的工作目录
  setwd(paste0("C:/Users/yuanfang/Desktop/ImmuneMedMR/07.metDiseaseMR/",new_folder_name))
  
  exposureFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/05.metMR/1.met.exposure.F10.csv"        #暴露数据文件
  outcomeFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/05.metMR/2.outcome.csv"  #结局数据文件
  
  
  #免疫细胞到代谢物的SNP文件
  beta1SNPfile=paste0("C:/Users/yuanfang/Desktop/ImmuneMedMR/06.immMetMR/",file_name)  
  
  
  
  #读取暴露数据输入文件
  rt=read.csv(exposureFile, header=T, sep=",", check.names=F)
  
  #提取目标代谢物的暴露数据
  i=paste0("beta2.", exposureID)
  singleExposureFile=paste0(i, ".exposure.csv")
  exposure_set=rt[rt$id.exposure==exposureID,]
  
  #去除免疫细胞到代谢物的SNP
  beta1SNP=read.csv(beta1SNPfile, header=T, sep=",", check.names=F)
  exposureOut=exposure_set[!(exposure_set$SNP %in% as.vector(beta1SNP[,"SNP"])),]
  write.csv(exposureOut, file=singleExposureFile, row.names=F)
  
  #读取这个代谢物的暴露数据
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
                                 filename="C:/Users/yuanfang/Desktop/ImmuneMedMR/05.metMR/2.outcome.csv", sep = ",",
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
  write.csv(outTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
  
  #MR-PRESSO异常值检测(偏倚的SNP)
  presso=run_mr_presso(dat)
  write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
  write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
  
  #孟德尔随机化分析
  mrResult=mr(dat)
  
 	#筛选IVW方法pvalue小于0.05的结果
		if(mrResult$pval[3]<0.05){
			#筛选五种方法OR方向一致的结果
			if(sum(mrResult$b>0)==nrow(mrResult) | sum(mrResult$b<0)==nrow(mrResult)){
			  #对结果进行OR值的计算
			  mrTab=generate_odds_ratios(mrResult)
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
			  
				rm(new_folder_name,exposureID,beta1SNPfile,rt,i,singleExposureFile,
				   exposure_set,beta1SNP,exposureOut,exposure_dat,outcome_data,dat,outTab,
				   presso,mrResult,heterTab,pleioTab,res_single)
				}
		  else
		  {
		    folder_path <- paste0("C:/Users/yuanfang/Desktop/ImmuneMedMR/07.metDiseaseMR/",new_folder_name)
		   unlink(folder_path, recursive = TRUE) 
		    
		    rm(new_folder_name,exposureID,beta1SNPfile,rt,i,singleExposureFile,
		       exposure_set,beta1SNP,exposureOut,exposure_dat,outcome_data,dat,outTab,
		       presso,mrResult,folder_path)
		  }
		}
  
}



