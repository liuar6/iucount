exitWithError=function(x){cat(x); quit()}
args=commandArgs(trailingOnly = T)
if (length(args)!=2) exitWithError("Exact two arguments should be provided.\n")

options(scipen=100)
df=read.table(args[1], sep="\t", stringsAsFactors=F,quote="")
df=df[df$V3=="exon", ]
exon_V5=sub('.*gene_id "([^"]*)".*', '\\1',df$V9)
exon_V6=sub('.*transcript_id "([^"]*)".*', '\\1',df$V9)
exon_V7=sub('.*gene_name "([^"]*)".*', '\\1',df$V9)
exon_V8=sub('.*transcript_name "([^"]*)".*', '\\1',df$V9)
exons=data.frame(V1=df$V1, V2=df$V4, V3=df$V5, V4=df$V7, V5=exon_V5, V6=exon_V6, V7=exon_V7, V8=exon_V8, stringsAsFactors=F)

V1=rep("", nrow(exons))
V2=rep(0, nrow(exons))
V3=rep(0, nrow(exons))
V4=rep("", nrow(exons))
V5=rep("", nrow(exons))
V6=rep("", nrow(exons))
V7=rep("", nrow(exons))
V8=rep("", nrow(exons))
tid=""
for (i in 1:nrow(exons)){
    if (tid!=exons[i, 6]){
        tid=exons[i, 6]
    } else{
        V1[i]=exons[i, 1]
        V4[i]=exons[i, 4]
        V5[i]=exons[i, 5]
        V6[i]=exons[i, 6]
        V7[i]=exons[i, 7]
        V8[i]=exons[i, 8]
        if (exons[i, 4]=="+") {V2[i]=lastEnd; V3[i]=exons[i, 2]-1}
        if (exons[i, 4]=="-") {V2[i]=exons[i, 3]; V3[i]=lastStart-1}
    }
    lastStart=exons[i, 2]
    lastEnd=exons[i, 3]
}
annotation=data.frame(V1, V2, V3, V4, V5, V6, V7, V8,stringAsFactors=F)
annotation=annotation[annotation$V4!="",]
#from here same file is created as extractintron.sh but run much slower
#some test
if (!all(annotation$V2-annotation$V3<0)) exitWithError("The gtf file is not valid.\n")
if (length(unique(annotation$V6))!=length(unique(paste0(annotation$V6,"_",annotation$V4)))) exitWithError("The gtf file is not valid.\n")
a=as.numeric(factor(annotation$V6, levels=unique(annotation$V6)))
if (!all(a[2:length(a)]-a[1:(length(a)-1)]>=0))exitWithError("The gtf file is not valid.\n")


#deduplicate
annotation$id=paste0(annotation$V1, "_", annotation$V2, "_", annotation$V3, "_", annotation$V4)
new_V5=tapply(annotation$V5, annotation$id, function(x) paste0(x[!duplicated(x)], collapse=","))
new_V6=tapply(annotation$V6, annotation$id, function(x) paste0(x[!duplicated(x)], collapse=","))
new_V7=tapply(annotation$V7, annotation$id, function(x) paste0(x[!duplicated(x)], collapse=","))
new_V8=tapply(annotation$V8, annotation$id, function(x) paste0(x[!duplicated(x)], collapse=","))
df1=data.frame(V5=new_V5, V6=new_V6, V7=new_V7, V8=new_V8, stringAsFactors=F)
df1$V1=unlist(lapply(strsplit(row.names(df1),"_"),function(x) x[[1]]))
df1$V2=unlist(lapply(strsplit(row.names(df1),"_"),function(x) x[[2]]))
df1$V3=unlist(lapply(strsplit(row.names(df1),"_"),function(x) x[[3]]))
df1$V4=unlist(lapply(strsplit(row.names(df1),"_"),function(x) x[[4]]))
row.names(df1)=1:nrow(df1)
df1=df1[paste0("V",1:8)]
write.table(df1,paste0(args[2]), quote=F, sep="\t", row.names=F,col.names=F)



