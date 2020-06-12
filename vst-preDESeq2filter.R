library("DESeq2")
cts<-read.table("CEL-Seq2_counting_report.txt", header=TRUE,row.names=1, sep="\t")
cts<-cts[1:(nrow(cts)-5),]
cts<-cts[-(grep("ERCC*",rownames(cts))),]
cts<-cts[(rowSums(cts)>10),]
celltype<-c(rep("epithelial",5),rep("flask",3),rep("globular",4),
            rep("inner",5),rep("pigment",4),rep("spicule",3))
coldata<-data.frame(row.names=colnames(cts),condition=celltype)

### Create DESeq2 object, vst transform data###
dds<-DESeqDataSetFromMatrix(countData = cts, colData = coldata,design=~condition)
dds$condition<-factor(dds$condition,levels=c("epithelial","flask","globular",
                                             "inner","pigment","spicule"))
vsd<- varianceStabilizingTransformation(dds, blind=TRUE)
vsdMat<- assay(vsd)
write.table(vsdMat,"CEL-Seq2_counting_report_vsd.txt",sep="\t")

### Using transformed data, filter lowest 25% var genes to improve DESeq (Sha et al 2016: between 0.3 & 0.15 optimal range to inc DEG detection) ###
vsdMat<-as.data.frame(vsdMat)
vsdMat$var = apply(vsdMat, 1, var)
vsdMatd = vsdMat[vsdMat$var >= quantile(vsdMat$var, c(.25)), ]
vsdMatd$var <- NULL
write.table(vsdMatd,"CEL-Seq2_counting_report_vsd_highest75%var.txt",sep="\t")

