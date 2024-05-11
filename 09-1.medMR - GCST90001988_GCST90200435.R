rm(list = ls()) #清除缓存
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/09.medMR/ebi-a-GCST90001988_GCST90200435")   

#三个输入文件
#修改1.免疫细胞到代谢物MR分析结果文件
beta1File="C:/Users/yuanfang/Desktop/ImmuneMedMR/06.immMetMR/beta1.ebi-a-GCST90001988_GCST90200435.table.MRresult.csv"     

#修改2.#代谢物到疾病MR分析结果文件
beta2File="C:/Users/yuanfang/Desktop/ImmuneMedMR/07.metDiseaseMR/ebi-a-GCST90001988_GCST90200435/beta2.GCST90200435.table.MRresult.csv"     

#修改3：免疫细胞到疾病MR分析结果文件
betaAllFile="C:/Users/yuanfang/Desktop/ImmuneMedMR/08.immDiseaseMR/ebi-a-GCST90001988/5.betaAll.ebi-a-GCST90001988.MRresult.csv"      


##########第一步：计算中介效应###########
#读取三个输入文件, 获取beta值
beta1RT=read.csv(beta1File, header=T, sep=",", check.names=F)
beta2RT=read.csv(beta2File, header=T, sep=",", check.names=F)
betaAllRT=read.csv(betaAllFile, header=T, sep=",", check.names=F)
b1=beta1RT[3,"b"]
se1=beta1RT[3,"se"]
b2=beta2RT[3,"b"]
se2=beta1RT[3,"se"]
bataAll=betaAllRT[3,"b"]

#计算中介效应(indirect effect)
beta12=b1 * b2
se12=sqrt((b1^2)*(se1^2)+(b2^2)*(se2^2))
ciLow=beta12-1.96*se12
ciHigh=beta12+1.96*se12
medEffect=paste0(signif(beta12, digits=3), "(", signif(ciLow, digits=3), ", ", signif(ciHigh, digits=3), ")")

#计算直接效应(direct effect)
betaDirect=bataAll-beta12
betaDirect

#中介效应所占的百分比
beta12Ratio=beta12/bataAll
ciLowRatio=ciLow/bataAll
ciHighRatio=ciHigh/bataAll
medRatioEffect=paste0(signif(beta12Ratio*100, digits=3), "%(", signif(ciLowRatio*100, digits=3), "%, ", signif(ciHighRatio*100, digits=3), "%)")

#计算pvalue
zValue=beta12/se12
pvalue=2*pnorm(q=abs(zValue), lower.tail=FALSE)

#输出结果的表格
outRT=data.frame(unique(beta1RT[,"exposure"]),
				  unique(beta2RT[,"exposure"]),
				  unique(betaAllRT[,"outcome"]),
                  medEffect, medRatioEffect, pvalue)
colnames(outRT)=c("Immune cell","Metabolite","outcome","Mediated effect", "Mediated proportion", "pvalue")
write.csv(outRT, file="MediatedEffect.csv", row.names=F)
##########第一步：计算中介效应###########



############第二步：绘制森林图############
#引用包
library(grid)
library(readr)
library(forestploter)

#读取孟德尔随机化分析的结果
files=c(beta1File,beta2File,betaAllFile)

data=data.frame()
for(i in files){
  rt=read.csv(i, header=T, sep=",", check.names=F)
  data=rbind(data, rt)
}
lineVec=cumsum(c(1,table(data[,c('exposure','outcome')])))

#对数据进行整理
data$' ' <- paste(rep(" ", 10), collapse = " ")
data$'OR(95% CI)'=ifelse(is.na(data$or), "", sprintf("%.3f (%.3f to %.3f)", data$or, data$or_lci95, data$or_uci95))
data$pval = ifelse(data$pval<0.001, "<0.001", sprintf("%.3f", data$pval))
data$exposure = ifelse(is.na(data$exposure), "", data$exposure)
data$nsnp = ifelse(is.na(data$nsnp), "", data$nsnp)
data2=data[,c('exposure','outcome')]
data[duplicated(data2),]$exposure=""
data[duplicated(data2),]$outcome=""
data$outcome[1] <- data$exposure[6]
#准备图形参数
tm <- forest_theme(base_size = 15,   #图形整体的大小
                   #可信区间的形状、线条类型、宽度、颜色、两端竖线高度
                   ci_pch = 16, ci_lty = 1, ci_lwd = 1.5, ci_col = "black", ci_Theight = 0.2, 
                   #参考线条的形状、宽度、颜色
                   refline_lty="dashed", refline_lwd=1, refline_col="grey20",
                   #x轴刻度字体的大小
                   xaxis_cex=0.8,
                   #脚注大小、颜色
                   footnote_cex = 0.6, footnote_col = "blue")

#绘制图形
plot <- forestploter::forest(data[, c("exposure","outcome","nsnp","method","pval"," ","OR(95% CI)")],
                             est = data$or,
                             lower = data$or_lci95,
                             upper = data$or_uci95,
                             ci_column = 6,     #可信区间所在的列
                             ref_line = 1,      #参考线条的位置
                             xlim = c(0.5, 1.5),    #X轴的范围
                             theme = tm,        #图形的参数
)

#修改图形中可信区间的颜色
boxcolor = c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#8491B4","#91D1C2","#DC0000","#7E6148")
boxcolor = boxcolor[as.numeric(as.factor(data$method))]
for(i in 1:nrow(data)){
  plot <- edit_plot(plot, col=6,row = i, which = "ci", gp = gpar(fill = boxcolor[i],fontsize=25)) # 改col，box的列
}

#设置pvalue的字体
pos_bold_pval = which(as.numeric(gsub('<',"",data$pval))<0.05)
if(length(pos_bold_pval)>0){
  for(i in pos_bold_pval){
    plot <- edit_plot(plot, col=5,row = i, which = "text", gp = gpar(fontface="bold"))  # 改col pvalue的列
  }
}

#在图形中增加线段
plot <- add_border(plot, part = "header", row =1,where = "top",gp = gpar(lwd =2))
plot <- add_border(plot, part = "header", row = lineVec, gp = gpar(lwd =1))
#设置字体大小, 并且将文字居中
plot <- edit_plot(plot, col=1:ncol(data),row = 1:nrow(data), which = "text", gp = gpar(fontsize=12))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),part="header",
                  x = unit(0.5, "npc"))
plot <- edit_plot(plot, col = 1:ncol(data), which = "text",hjust = unit(0.5, "npc"),
                  x = unit(0.5, "npc"))

#输出图形
pdf("forest.pdf", width=18, heigh=8)
print(plot)
dev.off()
