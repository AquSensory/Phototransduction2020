library("DESeq2")
cts=read.table("CEL-Seq2_counting_report_highest75%var.txt",header=TRUE,row.names=1,sep="\t")
celltype<-c(rep("epithelial",5),rep("flask",3),rep("globular",4),
            rep("inner",5),rep("pigment",4),rep("spicule",3))
coldata<-data.frame(row.names=colnames(cts),condition=celltype)

###For specific contrast, here e cells vs not-e cells (this is more sensitive than a model including all samples together)###  
coldata$condition<-sub("pigment","ne",coldata$condition)
coldata$condition<-sub("flask","ne",coldata$condition)
coldata$condition<-sub("inner","ne",coldata$condition)
coldata$condition<-sub("spicule","ne",coldata$condition)
coldata$condition<-sub("globular","ne",coldata$condition)

###Create DESeq2 object, fit GLM, lfcshrink to shrink logfold ∆s for small cts###
dds<-DESeqDataSetFromMatrix(countData = cts, colData = coldata,design=~condition)
dds$condition<-factor(dds$condition,levels=c("epithelial","ne"))
dds<-DESeq(dds)
alpha<-0.1 
res<-results(dds,contrast=c("condition","epithelial","ne"),alpha=alpha) 
summary(res)
reslfc<-lfcShrink(dds,contrast=c("condition","epithelial","ne"),res=res,alpha=alpha)
summary(reslfc)
plotMA(reslfc,ylim=c(-3,3))
reslfcord<-reslfc[order(reslfc$padj),]
reslfcsig<-subset(reslfcord,padj<alpha)
write.csv(as.data.frame(reslfcsig), file="DEG_E_sig0.1.csv")
