#本程序一共有三处修改
rm(list = ls()) #清除缓存
library("ieugwasr")
library("VariantAnnotation")
library("gwasglue")
library("TwoSampleMR")

#修改1：设置工作目录
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/02.immMR")  
#修改2：设置疾病名称
diseaseName="Gestational diabetes"     
#修改3：设置结局ID
diseaseID="finn-b-O15_PREG_DM"       


#############第一步：筛选暴露F>0的工具变量#############
data0 <- read.csv("C:/Users/yuanfang/Desktop/ImmuneMedMR/01.exposure/exposure_data.csv", header=T, sep=",", check.names=F)

#计算F检验值
data0$R2<-(2*data0$beta.exposure*data0$beta.exposure*data0$eaf.exposure*(1-data0$eaf.exposure)/(2*data0$beta.exposure*data0$beta.exposure*data0$eaf.exposure*(1-data0$eaf.exposure)+2*data0$se.exposure*data0$se.exposure*data0$samplesize.exposure*data0$eaf.exposure*(1-data0$eaf.exposure)))     #计算R2
data0$F<-data0$R2*(data0$samplesize.exposure-2)/(1-data0$R2)     #计算F检验值

#根据F值>10对数据进行过滤, 删除弱工具变量
data <- data0[as.numeric(data0$F)>10,]
write.csv(data, file="1.imm.exposure.F10.csv", row.names=F)
#############第一步：筛选暴露F>0的工具变量#############



#############第二步：合并暴露数据和结局数据，并做MR#############
#读取暴露数据
exposure=read_exposure_data(filename="1.imm.exposure.F10.csv",
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
outcome <- extract_outcome_data(snps=exposure$SNP, outcomes=diseaseID)
write.csv(outcome, file="2.outcome.csv", row.names=F)

#将暴露数据和结局数据合并
outcome$outcome <- diseaseName
dat <- harmonise_data(exposure, outcome)

#输出用于孟德尔随机化的工具变量
outTab <- dat[dat$mr_keep=="TRUE",]
write.csv(outTab, file="3.SNP.csv", row.names=F)

#孟德尔随机化分析
mrResult <- mr(dat)

#对结果进行OR值的计算
mrTab <- generate_odds_ratios(mrResult)
write.csv(mrTab, file="4.MRresult.csv", row.names=F)

#异质性检验
heterTab <- mr_heterogeneity(dat)
write.csv(heterTab, file="5.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab <- mr_pleiotropy_test(dat)
write.csv(pleioTab, file="6.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="7.scatter plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, dat)
dev.off()

#森林图
res_single <- mr_singlesnp(dat)      #得到每个工具变量对结局的影响
pdf(file="8.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="9.funnel plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="10.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()
#############第二步：合并暴露数据和结局数据，并做MR#############



#############第三步：根据P值、五种OR方法一致原则对MR的结果进行过滤#############
pvalFilter <- 0.01     #pvalue过滤条件(0.05   0.01   0.001)

#读取孟德尔随机化的结果文件
rt <- mrTab

#根据pvalue对孟德尔随机化分析的结果进行过滤
ivwRT <- rt[rt$method=="Inverse variance weighted",]
ivwRT <- ivwRT[ivwRT$pval<pvalFilter,]

#提取五种方法OR值方向一致的免疫细胞
ivw <- data.frame()
for(cell in unique(ivwRT$exposure)){
  cellData=rt[rt$exposure==cell,]
  if(sum(cellData$or>1)==nrow(cellData) | sum(cellData$or<1)==nrow(cellData)){
    ivw=rbind(ivw, ivwRT[ivwRT$exposure==cell,])
  }
}

#读取多效性的结果文件
pleRT <- pleioTab
#剔除pvalue小于0.05的免疫细胞
pleRT <- pleRT[pleRT$pval>0.05,]
cellLists <- as.vector(pleRT$exposure)
IVWfilter0 <- ivw[ivw$exposure %in% cellLists,]

#读取异质性检验,剔除pvalue小于0.05的免疫细胞
heter <- heterTab[heterTab$Q_pval>0.05,]
cellLists1 <- as.vector(heter$exposure)
IVWfilter <- IVWfilter0[IVWfilter0$exposure %in% cellLists1,]
length(IVWfilter$id.exposure)
write.csv(IVWfilter, file="11.immune-disease.IVWfilter.csv", row.names=F)
#############第三步：根据P值、五种OR方法一致原则对MR的结果进行过滤#############
#异质性检验
heterTab1 <- mr_heterogeneity(dat[dat$id.exposure %in% IVWfilter$id.exposure,])
write.csv(heterTab1, file="5.heterogeneity1.csv", row.names=F)

#多效性检验
pleioTab1 <- mr_pleiotropy_test(dat[dat$id.exposure %in% IVWfilter$id.exposure,])
write.csv(pleioTab1, file="6.pleiotropy1.csv", row.names=F)

#绘制散点图
pdf(file="7.scatter plot1.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult[mrResult$id.exposure %in% IVWfilter$id.exposure,], dat[dat$id.exposure %in% IVWfilter$id.exposure,])
dev.off()

#森林图
res_single <- mr_singlesnp(dat[dat$id.exposure %in% IVWfilter$id.exposure,])      #得到每个工具变量对结局的影响
pdf(file="8.forest1.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="9.funnel plot1.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single[res_single$id.exposure %in% IVWfilter$id.exposure,])
dev.off()

#留一法敏感性分析
pdf(file="10.leaveoneout1.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat[dat$id.exposure %in% IVWfilter$id.exposure,]))
dev.off()

save(list = ls(), file = "immune-disease MR.RData")

#load("immune-disease MR.RData")
