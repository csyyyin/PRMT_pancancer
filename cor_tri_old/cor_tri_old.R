
setwd("D:/PRMT/cor_tri_old")

library(ggplotify)

tumor_list = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

my_data = read.csv('D:/PRMT/cor_tri_old/PRMT.csv',header=T,check.names=F)

out_put_list <- c()
#trace('corrplot',edit=T)
for(i in tumor_list)
{
  res <- my_data[grep(i,my_data$CODE),]
  res <- res[grep("Tumor",res$Group),]
  res <- res[,-1:-2]
  rownames(res) = gsub('\n','',res[,length(res)])
  res <- res[,1:length(res)-1]
  library(corrplot)
  final <- cor(res)
  pdf(paste(i,".pdf",sep=''),width=length(colnames(final))/2,height=length(colnames(final))/2)
  
  temp = corrplot(final, type = "upper", order = 'alphabet', 
                  tl.col = "black", tl.srt = 90,tl.cex=2,method = "square",na.label=i)
  temp = temp$corr
  write.csv(temp,file=paste0(i,'.csv'))
  dev.off()
  corrplot(final, type = "upper", order = "alphabet", 
           tl.col = "black",tl.cex=0.5,method = "square",upper = "pie", tl.srt = 90)
  eval(parse(text = paste(i,'<- recordPlot()')))
  out_put_list <- append(out_put_list,c(i,NULL))
}

pdf("fig1.pdf",width=20,height=30)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=5,labels=out_put_list,label_size=40))',sep='')))

dev.off()
