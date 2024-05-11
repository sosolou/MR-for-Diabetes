#本程序需修改两处
rm(list = ls()) #清除缓存
library(ieugwasr)
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)

setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/05.metMR")     #设置工作目录

#修改1：设置结局ID
outcomeID="finn-b-O15_PREG_DM"           
#修改2：设置疾病名称
outcomeName="Gestational diabetes"     


#######第一步：读取暴露数据，根据F值筛选与暴露强相关的工具变量###########
dat=read.csv("C:/Users/yuanfang/Desktop/ImmuneMedMR/04.metData/met.exposure_data.csv", header=T, sep=",", check.names=F)

#计算F检验值
dat$R2<-(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)/(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)+2*dat$se.exposure*dat$se.exposure*dat$samplesize.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)))     #计算R2
dat$F<-dat$R2*(dat$samplesize.exposure-2)/(1-dat$R2)     #计算F检验值

#根据F值>10对数据进行过滤, 删除弱工具变量
metData=dat[as.numeric(dat$F)>10,]
write.csv(metData, file="1.met.exposure.F10.csv", row.names=F)
#######第一步：读取暴露数据，根据F值筛选与暴露强相关的工具变量###########



#######第二步：合并暴露数据和结局数据###########
#读取暴露数据
exposure_dat=read_exposure_data(filename="1.met.exposure.F10.csv",
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
outcomeData=extract_outcome_data(snps=exposure_dat$SNP, outcomes=outcomeID)
write.csv(outcomeData, file="2.outcome.csv", row.names=F)

#将暴露数据和结局数据合并
outcomeData$outcome=outcomeName
data=harmonise_data(exposure_dat, outcomeData)
#######第二步：合并暴露数据和结局数据###########



#######第三步：进行MR分析###########
#输出用于孟德尔随机化的工具变量
outSNP=data[data$mr_keep=="TRUE",]
write.csv(outSNP, file="3.SNP.csv", row.names=F)

#孟德尔随机化分析
mrResult=mr(data)

#对结果进行OR值的计算
mrTab=generate_odds_ratios(mrResult)
write.csv(mrTab, file="4.MRresult.csv", row.names=F)

#异质性分析
heterTab=mr_heterogeneity(data)
write.csv(heterTab, file="5.heterogeneity.csv", row.names=F)

#多效性检验
pleioTab=mr_pleiotropy_test(data)
write.csv(pleioTab, file="6.pleiotropy.csv", row.names=F)

#绘制散点图
pdf(file="7.scatter_plot.pdf", width=7.5, height=7)
mr_scatter_plot(mrResult, data)
dev.off()

#森林图
res_single=mr_singlesnp(data)      #得到每个工具变量对结局的影响
pdf(file="8.forest.pdf", width=7, height=5.5)
mr_forest_plot(res_single)
dev.off()

#漏斗图
pdf(file="9.funnel_plot.pdf", width=7, height=6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="10.leaveoneout.pdf", width=7, height=5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(data))
dev.off()
#######第三步：进行MR分析###########


#############第四步：根据P值、五种OR方法一致原则对MR的结果进行过滤#############
pvalFilter <- 0.01      #pvalue过滤条件(0.05   0.01   0.001)

#读取孟德尔随机化的结果文件
rt <- mrTab

#根据pvalue对孟德尔随机化分析的结果进行过滤
ivwRT <- rt[rt$method=="Inverse variance weighted",]
ivwRT <- ivwRT[ivwRT$pval<pvalFilter,]

#提取五种方法OR值方向一致的代谢物
ivw=data.frame()
for(metabolite in unique(ivwRT$exposure)){
  metData=rt[rt$exposure==metabolite,]
  if(sum(metData$or>1)==nrow(metData) | sum(metData$or<1)==nrow(metData)){
    ivw=rbind(ivw, ivwRT[ivwRT$exposure==metabolite,])
  }
}

#读取多效性的结果文件
pleRT <- pleioTab
#剔除pvalue小于0.05的代谢物
pleRT=pleRT[pleRT$pval>0.05,]
metLists=as.vector(pleRT$exposure)
outTab0=ivw[ivw$exposure %in% metLists,]
row.names(outTab0)=outTab0[,1]


#读取异质性检验,剔除pvalue小于0.05的代谢物
heter <- heterTab[heterTab$Q_pval>0.05,]
metLists1 <- as.vector(heter$exposure)
outTab <- outTab0[outTab0$exposure %in% metLists1,]
length(outTab$id.exposure)


#在结果中添加代谢物的下载链接
metFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/04.metData/metExposureID.txt"
metRT=read.table(metFile, header=T, sep="\t", check.names=F, row.names=1,  comment.char="", quote="")
sameID=intersect(row.names(outTab), row.names(metRT))
outTab=cbind(outTab[sameID,,drop=F], metRT[sameID,"URL",drop=F])
write.csv(outTab, file="11.metabolite-disease.IVWfilter.csv", row.names=F)
#############第四步：根据P值、五种OR方法一致原则对MR的结果进行过滤#############



save(list = ls(), file = "met-disease MR.RData")

#load("met-disease MR.RData")



