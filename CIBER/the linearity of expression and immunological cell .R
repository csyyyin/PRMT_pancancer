setwd("~/projects/Tcell-pancancer")
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)

tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')

out_put_list <- c()
for (tumor_id in tumor_list)
{
  print(tumor_id)
  
  dat_input <- read.csv("~/diff/data/CIBER.csv",check.names = F,header=T)
  
  dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
  
  dat_input <- dat_input[grep("Tumor",dat_input$Group),]
  
  dat_input <- subset(dat_input,select=-c(Group,CODE))
  
  dat_input <- dat_input[!duplicated(dat_input$SampleName),]
  
  rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))
  
  dat_input <- subset(dat_input,select=-c(SampleName))
  
  colscluster = length(colnames(dat_input))/3
  
  
  data_all_tumor <- read.csv(paste("~/projects/Tcell/data/Tcell_tpm_exp.csv",sep=''),check.names = F)
  
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  
  if (length(colnames(dat_one_tumor))==0)
  {
    next
  }
  
  rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))
  
  dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]
  
  dat_one_tumor <- dat_one_tumor[grep("Tumor",dat_one_tumor$Group),]
  
  dat_one_tumor <- subset(dat_one_tumor,select=-c(Group,SampleName,CODE))
  
  one_tumor_sample <- unlist(rownames(dat_one_tumor))
  
  all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))
  
  dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]
  
  dat_im <- dat_input[match(all_name,rownames(dat_input)),]
  
  if(length(rownames(dat_gene))==0)
  {
    next
  }
  data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
  
  temp_data = data.frame()
  
  for(i in colnames(as.data.frame(data.corr$r)))
  {
    for(j in rownames(as.data.frame(data.corr$r)))
    {
      eval(parse(text=paste0("temp_row=data.frame(Row='",i,"',Col='",j,"',rvalue=","as.data.frame(data.corr$r)[,which(colnames(as.data.frame(data.corr$r))=='",i,"')][which(rownames(as.data.frame(data.corr$r))=='",j,"')]",",pvalue=","as.data.frame(data.corr$p)[,which(colnames(as.data.frame(data.corr$p))=='",i,"')][which(rownames(as.data.frame(data.corr$p))=='",j,"')])")))
      temp_data = rbind(temp_data,temp_row)
    }
  }
  pdf(paste('CIBER/corr/',tumor_id,".pdf",sep=''),width = 15,height = 8)
  p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = abs(-log10(pvalue)))) +
    geom_point() +  # 使用点表示数据
    scale_color_gradient2(low = '#00CC00', high = '#FF0000') +  # 调整颜色渐变,
    scale_size_continuous(range = c(1,10)) +  # 调整点的大小范围
    labs(title=tumor_id,x = "Gene", y = "", color = "R", size = "-log10(P)") +  # 添加轴标签
    ggtitle(tumor_id)+
    theme_bw()+
    theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
  print(p)
  dev.off()
  
  #拼图
  eval(parse(text = paste0(tumor_id,'= p')))
  out_put_list <- append(out_put_list,tumor_id)
}
##拼图
pdf("fig1.pdf",width = length(colnames(dat_input))*3.5,height = length(colnames(dat_one_tumor))*length(out_put_list)/20)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=5))',sep='')))

dev.off()
