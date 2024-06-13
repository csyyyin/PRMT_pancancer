
library("pheatmap")
library("jsonlite")

setwd(dir = "D:/PRMT/check")
temp = list.files(pattern="*.csv")


df = read.csv('D:/PRMT/check/p.csv',header=T,row.names=1)
df = replace(df,is.na(df),1)


paletteLength = 1000

getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

sig.mat <- matrix(sapply(as.matrix(df), getSig), nrow=nrow(as.matrix(df)))
str(sig.mat)


tran = function(temp)
{
  temp_p = c()
  for(i in temp)
  {
    if(i > 0)
    {
      if(i < 0.000000001)
      {
        i =  0.000000001
      }
      temp_p = append(temp_p,-log10(i))
    }
    else
    {
      if(i > -0.000000001)
      {
        i =  -0.000000001
      }
      temp_p = append(temp_p,log10(abs(i)))
    }
  }
  return(temp_p)
}

df_temp = sapply(df,tran)
rownames(df_temp) = rownames(df)
df = df_temp


#0-0.05
#myColor <- colorRampPalette(c( "red3", "white"))(paletteLength)
#myBreaks <- seq(0,0.05,length.out=paletteLength)

# #yzx
# myColor <- colorRampPalette(c( "#5B9C4B", "white","#5B9C4B"))(paletteLength)
# #gx
# myColor <- colorRampPalette(c( "#86C06C", "white","#86C06C"))(paletteLength)
# #normal tumor
# myColor <- colorRampPalette(c( "#B8DCA1", "white","#B8DCA1"))(paletteLength)
# CIBER
# myColor <- colorRampPalette(c( '#8197C6', "white",'##507AAF'))(paletteLength)

#reg
#myColor <- colorRampPalette(c( '#D78851', "white",'#507AAF'))(paletteLength)
#check
myColor <- colorRampPalette(c( '#9AD0AD', "white",'#9999FF'))(paletteLength)

myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))

pdf('heatmap.pdf',length(colnames(df))/2,length(rownames(df))/2)

xx <- pheatmap(df,
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",main='-log10(P)',number_color='black',border_color = "black", cluster_rows=F,cluster_cols=F, cellwidth = 15,cellheight = 15,display_numbers=sig.mat)
print(xx)
dev.off()

