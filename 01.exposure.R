#该程序不需要运行
library(TwoSampleMR)

inputFile="exposureID.txt"      #暴露数据id文件
setwd("C:/Users/yuanfang/Desktop/ImmuneMedMR/01.exposure")     #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)

#对免疫细胞数据ID进行循环
outTab=data.frame()
for(id in rt$ID) {
	expoData=extract_instruments(id,
                                 p1 = 1e-5, p2 = 1e-5,
                                 clump = T,
                                 kb = 10000, r2 = 0.001)
    outTab=rbind(outTab, expoData)
}
write.csv(outTab, file="exposure_data.csv", row.names=F)




