###########需要修改一处##############
rm(list = ls()) #清除缓存
library(TwoSampleMR)
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/03.immReverse") 
     

#读取结局数据文件(IVW方法过滤的结果文件)
rt=read.csv("C:/Users/yuanfang/Desktop/ImmuneMedMR/02.immMR/11.immune-disease.IVWfilter.csv", header=T, sep=",", check.names=F)
row.names(rt)=rt$id.exposure

#获取暴露数据(以疾病为暴露)
dat=extract_instruments(outcomes = unique(rt$id.outcome), 
                        p1=5e-08, p2=5e-08, clump=TRUE,access_token = NULL)

#修改1：疾病sampleSize，如果有这个变量，就把下面这句删掉
dat$samplesize.exposure <- 16379684

#计算F检验值
dat$R2<-(2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)/
           (2*dat$beta.exposure*dat$beta.exposure*dat$eaf.exposure*
              (1-dat$eaf.exposure)+2*dat$se.exposure*dat$se.exposure*
              dat$samplesize.exposure*dat$eaf.exposure*(1-dat$eaf.exposure)))     #计算R2
dat$F<-dat$R2*(dat$samplesize.exposure-2)/(1-dat$R2)     #计算F检验值

#根据F值>10对数据进行过滤, 删除弱工具变量
exposure_dat=dat[as.numeric(dat$F)>10,]


#对结局数据进行循环(免疫细胞)
revPvalVec=c()
for(i in row.names(rt)){
  #提取结局数据
  outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes=i)
  
  #将暴露数据和结局数据合并
  outcome_dat$outcome=rt[i, "exposure"]
  dat <- harmonise_data(exposure_dat, outcome_dat)
  
  #MR-PRESSO异常值检测(偏倚的SNP)
  #presso=run_mr_presso(dat)
  #write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file=paste0(i, ".table.MR-PRESSO_Global.csv"))
  #write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file=paste0(i, ".table.MR-PRESSO_Outlier.csv"))
  
  #孟德尔随机化分析
  mrResult=mr(dat)
  mrTab=generate_odds_ratios(mrResult)
  #输出孟德尔随机化分析的结果
  write.csv(mrTab, file=paste0(i, ".table.MRresult.csv"), row.names=F)
  revPvalVec=c(revPvalVec, mrResult$pval[3])
  
  #对IVW方法pvalue小于0.05的结果可视化
  if(mrResult$pval[3]<0.05){
    #输出用于孟德尔随机化的工具变量
    outTab=dat[dat$mr_keep=="TRUE",]
    write.csv(outTab, file=paste0(i, ".table.SNP.csv"), row.names=F)
    
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

#输出反向孟德尔随机化分析的结果文件
outTab=cbind(rt, revPvale=revPvalVec)
outTab=outTab[as.numeric(outTab$revPvale)>0.05,]
write.csv(outTab, file="immune-disease.REVfilter.csv", row.names=F)

save(list = ls(), file = "immune-disease MRreverse.RData")

#load("immune-disease MRreverse.RData")



