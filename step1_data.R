rm(list = ls())
library(tidyverse)
library(readxl)
library(GEOquery)
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# 180081
gse_number = "GSE180081"
eSet <- getGEO(gse_number, destdir = '.', getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
pd_jia <- pData(eSet)
exp_jia <- read_excel("GSE180081_CathDx_HG38_realign_RPKMs81.xlsx")
dat <- exp_jia[,-1];dat$gene <- dat$`Gene Name`
symbol_all<- aggregate( . ~ gene, data = dat,max)
gene_name <- symbol_all$gene
exp_jia <- symbol_all[,-c(1,98)]
exp_jia <- as.data.frame(lapply(exp_jia,as.numeric))
rownames(exp_jia) <- gene_name;colnames(exp_jia) <- pd_jia$geo_accession
i <- c("geo_accession","characteristics_ch1","characteristics_ch1.4","characteristics_ch1.5",
       "characteristics_ch1.11","characteristics_ch1.13","characteristics_ch1.14",
       "characteristics_ch1.15","characteristics_ch1.16")
pd_jia <- pd_jia[,i]
colnames(pd_jia) <- c("ID","status","gender","age","bmi","smoking","htn","dyslipidemia","diabetes")
exp_jia <- fpkmToTpm(exp_jia)
q1 <- data.frame(str_split_fixed(pd_jia$status,":",n = 2))
pd_jia$status <- q1$X2
q2 <- data.frame(str_split_fixed(pd_jia$gender,":",n = 2))
pd_jia$gender <- q2$X2
q3 <- data.frame(str_split_fixed(pd_jia$age,":",n = 2))
pd_jia$age <- as.numeric(q3$X2)
q4 <- data.frame(str_split_fixed(pd_jia$bmi,":",n = 2))
pd_jia$bmi <- q4$X2
q5 <- data.frame(str_split_fixed(pd_jia$smoking,":",n = 2))
pd_jia$smoking <-  q5$X2
q6 <- data.frame(str_split_fixed(pd_jia$htn,":",n = 2))
pd_jia$htn <- q6$X2
q7 <- data.frame(str_split_fixed(pd_jia$dyslipidemia,":",n = 2))
pd_jia$dyslipidemia <- q7$X2
q8 <- data.frame(str_split_fixed(pd_jia$diabetes,":",n = 2))
pd_jia$diabetes <- q8$X2
Group_jia <- ifelse(pd_jia$status == " LOW", "Control", "CAD")
Group_jia = factor(Group_jia,levels = c("Control","CAD"))
table(Group_jia)



# 180082
gse_number = "GSE180082"
eSet <- getGEO(gse_number, destdir = '.', getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
pd_yi <- pData(eSet)
exp_yi <- read_excel("GSE180082_TruCAD_SeqLL_RPKMs82.xlsx")
dat <- exp_yi[,-1];dat$gene <- dat$`Gene Name`
symbol_all<- aggregate( . ~ gene, data = dat,max)
gene_name <- symbol_all$gene
exp_yi <- symbol_all[,-c(1,82)]
exp_yi <- as.data.frame(lapply(exp_yi,as.numeric))
rownames(exp_yi) <- gene_name;colnames(exp_yi) <- pd_yi$geo_accession
i <- c("geo_accession","cad status:ch1","gender:ch1")
pd_yi <- pd_yi[,i]
colnames(pd_yi) <- c("ID","status","gender")
exp_yi <- fpkmToTpm(exp_yi)
Group_yi <- ifelse(pd_yi$status == "LOW", "Control", "CAD")
Group_yi = factor(Group_yi,levels = c("Control","CAD"))
table(Group_yi)
pd_yi$event <- ifelse(Group_yi == "Control", 0, 1)
save(exp_jia,pd_jia,Group_jia,exp_yi,pd_yi,Group_yi, file = "step1_data.Rdata")


############################################################
### step2
rm(list = ls())
load("step1_data.Rdata")
library(limma)
library(edgeR)
exp <- exp_jia
dge <- DGEList(counts = exp)
dge <- calcNormFactors(dge)
design <- model.matrix(~Group_jia)
v <- voom(dge,design, normalize="quantile")
design <- model.matrix(~Group_jia)
fit <- lmFit(v, design)
fit= eBayes(fit)
DEG = topTable(fit, coef=2, n=Inf)
DEG = na.omit(DEG)
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -1);table(k1)
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > 1);table(k2)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
###   
library(ggplot2)
library(tinyarray)
library(edgeR)
exp[1:4,1:4]
# cpm 
dat = log2(cpm(exp)+1)
pca.plot = draw_pca(dat,Group_jia);pca.plot
pca.plot.3D = draw_pca(dat,Group_jia, style = "3D" );pca.plot.3D
cg = rownames(DEG)[DEG$change !="NOT"]
h = draw_heatmap(dat[cg,],Group_jia,n_cutoff = 2);h
v = draw_volcano(DEG,pkg = 3,logFC_cutoff = 1);v
library(patchwork)
h + v +plot_layout(guides = 'collect') &theme(legend.position = "none")
# 
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd","#f87669")),
                       labels = c("Control","CAD"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
m = Heatmap(t(scale(t(exp[cg,]))),name = " ",
            col = col_fun,
            top_annotation = top_annotation,
            column_split = Group_jia,
            show_heatmap_legend = T,
            border = F,
            show_column_names = F,
            show_row_names = F,
            column_title = NULL)
m
save(DEG, cg, file = "step2_deg.Rdata")
tiff(filename = "D:// /CAD-2//plot//1.1deg_pca.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//1.1deg_pca.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//1.1deg_pca(3D).tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//1.1deg_pca(3D).pdf",width = 10, height = 6.5, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//1.2deg_vol.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//1.2deg_vol.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//1.3deg_heat.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//1.3deg_heat.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()


#########################################################################
### step3
rm(list = ls())
load("step1_data.Rdata");rm(exp_yi);rm(pd_yi);rm(Group_yi)
exp <- exp_jia[Group_jia == "CAD"];pd <- pd_jia[colnames(exp),-2]
pd$gender <- ifelse(pd$gender == " M", 1, 0)
pd$age <- ifelse(pd$age > 60, 1, 0)
pd$bmi <- ifelse(pd$bmi == " L", 0, 1)
pd$smoking <- ifelse(pd$smoking == " N", 0, 1)
pd$htn <- ifelse(pd$htn == " N", 0, 1)
pd$dyslipidemia <- ifelse(pd$dyslipidemia == " N", 0, 1)
pd$diabetes <- ifelse(pd$diabetes == " N", 0, 1)
### 
femData <- exp
CD_gene <- c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","MTF1","GLS","CDKN2A")
CD_matrix <- femData[rownames(femData) %in% CD_gene, ]
library(ConsensusClusterPlus)
pac_matrix<-as.matrix(CD_matrix)
maxK= 6
results = ConsensusClusterPlus(pac_matrix,
                               maxK=maxK,
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               title= "SKCM_fe-relate",
                               clusterAlg="km",
                               distance="euclidean",
                               seed=123,
                               plot="pdf",
                               writeTable = T
)
pac_group_2<- as.data.frame(results[[2]][["consensusClass"]])
names(pac_group_2)<- "pac_group"
pac_group_2$ID<- rownames(pac_group_2)
pd <- merge(pac_group_2,pd,by="ID")
rownames(pd)<- pd$ID
Group <- ifelse(pd$pac_group == 1, "cluster_1", "cluster_2")
Group = factor(Group,levels = c("cluster_1", "cluster_2"))
table(Group)
#save(exp,pd,Group,file = "step3_CU.Rdata")


#########################################################################
### step4 CU-DEG
rm(list = ls())
load("step3_CU.Rdata")
library(limma)
library(edgeR)
exp = log2(exp+1)
dge <- DGEList(counts = exp)
dge <- calcNormFactors(dge)
design <- model.matrix(~Group)
v <- voom(dge,design, normalize="quantile")
design <- model.matrix(~Group)
fit <- lmFit(v, design)
fit= eBayes(fit)
DEG = topTable(fit, coef=2, n=Inf)
DEG = na.omit(DEG)
k1 = (DEG$P.Value < 0.05)&(DEG$logFC < -1);table(k1)
k2 = (DEG$P.Value < 0.05)&(DEG$logFC > 1);table(k2)
DEG$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(DEG$change)
head(DEG)
###   
library(ggplot2)
library(tinyarray)
exp[1:4,1:4]
# cpm 
dat = log2(cpm(exp)+1)
pca.plot = draw_pca(dat,Group);pca.plot
pca.plot.3D = draw_pca(dat,Group, style = "3D" );pca.plot.3D
cg = rownames(DEG)[DEG$change !="NOT"]
h = draw_heatmap(dat[cg,],Group,n_cutoff = 2);h
v = draw_volcano(DEG,pkg = 3,logFC_cutoff = 1);v
library(patchwork)
h + v +plot_layout(guides = 'collect') &theme(legend.position = "none")
# 
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#2fa1dd","#f87669")),
                       labels = c("cluster_1","cluster_2"),
                       labels_gp = gpar(col = "white", fontsize = 12)))
m = Heatmap(t(scale(t(exp[cg,]))),name = " ",
            col = col_fun,
            top_annotation = top_annotation,
            column_split = Group,
            show_heatmap_legend = T,
            border = F,
            show_column_names = F,
            show_row_names = F,
            column_title = NULL)
m
DEG_CU <- DEG; cg_Cu <- cg
save(DEG_CU, cg_Cu, file = "step4_degcu.Rdata")
tiff(filename = "D:// /CAD-2//plot//3.1deg_Cu_pca.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//3.1deg_Cu_pca.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//3.1deg_Cu_pca(3D).tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//3.1deg_Cu_pca(3D).pdf",width = 10, height = 6.5, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//3.2deg_Cu_vol.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//3.2deg_Cu_vol.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//3.3deg_Cu_heat.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//3.3deg_Cu_heat.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()


######################################################################
### step5   
rm(list = ls())
load("step2_deg.Rdata"); load("step3_CU.Rdata"); load("step4_degcu.Rdata")
cg1 = cg;cg2 = cg_Cu
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}

up = intersect(UP(DEG),UP(DEG_CU))
down = intersect(DOWN(DEG),DOWN(DEG_CU))
#degx_co <- c(up,down)
# 1.
#write.table(degx_co,
#            file="diffgene.txt",
#            row.names = F,
#            col.names = F,
#            quote = F)
# 
# 2.
#DEG$symbol <- rownames(DEG);DEG_CU$symbol <- rownames(DEG_CU)
#p1 = DEG[DEG$change != "stable",c("symbol","logFC")]
#p2 = DEG_CU[DEG_CU$change != "stable",c("symbol","logFC")]
#head(p1)
#write.table(p1,file = "deg_1.txt",sep = "\t", quote = F, row.names = F)
#write.table(p2,file = "deg_2.txt",sep = "\t", quote = F, row.names = F)
dat = log2(cpm(exp)+1)

#
up_genes = list(state = UP(DEG),
                cluster = UP(DEG_CU))

down_genes = list(state = DOWN(DEG),
                  cluster = DOWN(DEG_CU))

up.plot <- draw_venn(up_genes,"UPgene")
down.plot <- draw_venn(down_genes,"DOWNgene")

#

library(patchwork)
up.plot + down.plot
#save(up, down, file = "step5_CO_deg.Rdata")
tiff(filename = "D:// /CAD-2//plot//4.co_deg.tiff", res = 300, width = 2000, height = 1000,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//4.co_deg.pdf",width = 10, height = 5, onefile = FALSE)
dev.off()


#########################################################################
### step6    
rm(list = ls())
library(tidyverse)
load("step1_data.Rdata"); load("step5_CO_deg.Rdata")
exprSet <- exp_jia[c(up, down),]
pd <- pd_jia;
pd$event <- as.numeric(ifelse(str_detect(pd$status,"MID+"),
                              "1",
                              "0"))
pd$gender <- as.numeric(ifelse(pd$gender == 'M',0,1))
pd$bmi <- as.numeric(ifelse(pd$bmi == 'L',0,1))
pd$smoking <- as.numeric(ifelse(pd$smoking == 'N',0,1))
pd$htn <- as.numeric(ifelse(pd$htn == 'N',0,1))
pd$dyslipidemia <- as.numeric(ifelse(pd$dyslipidemia == 'N',0,1))
pd$diabetes <- as.numeric(ifelse(pd$diabetes == 'N',0,1))
x=as.data.frame(t(exprSet));x$ID <- rownames(x)
pdnew <- pd[,-c(2,10)]
x <- left_join(x,pdnew)
rownames(x) <- x$ID;x <- select(x,-'ID')
x <- as.matrix(x)
y=pd$event
library(glmnet)
#
#set.seed(1006)
#cvfit <- cv.glmnet(x=x, y=y, family = "binomial")
#save(cvfit, file = "cvfitR.Rdata")
load("cvfitR.Rdata")
plot(cvfit)

#
fit <- glmnet(x=x, y=y)
plot(fit,xvar = "lambda")
#model_lasso_min <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y,lambda=cvfit$lambda.min,family="binomial", nlambda=100, alpha=1)
#head(model_lasso_min$beta,20)
#choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
#length(choose_gene_min)
length(choose_gene_1se)
#save(choose_gene_min,file = "lasso_choose_gene_min.Rdata")
lasso.prob <- predict(cvfit, newx=x, s= cvfit$lambda.min)
re=cbind(event = pd$event ,lasso.prob)
re=as.data.frame(re)
colnames(re)=c('event','prob_min')
re$event <- ifelse(re$event == "0", "Control", "CAD")
re$event=as.factor(re$event)
head(re)
library(ggpubr) 
re1 <- re;colnames(re1) <- c("Group","LASSO_score");re1$Group <- ifelse(re1$Group == "CAD", "high", "low")
p = ggplot(dat = re1,aes(Group, LASSO_score))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means(method = "t.test")+
  theme_classic()

p

tiff(filename = "D:// /CAD-2//plot//6.1cv_fit.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//6.1cv_fit.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//6.2lambda.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//6.2lambda.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//6.3lasso_mod.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//6.3lasso_mod.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()

### logistic 
#rm(list = ls())
#load("step1.Rdata");load("lasso_choose_gene_min.Rdata")
log_exp <- exprSet[choose_gene_1se,]
log1<-glm(pd$event~.,family=binomial(link='logit'),data = as.data.frame(t(log_exp)))
summary(log1)
log2<-step(log1)
summary(log2)
pred<-predict(log2)
prob<-exp(pred)/(1+exp(pred))
re=cbind(pd$event ,prob)
head(re)
re=as.data.frame(re)
colnames(re)=c('event','prob_min')
re$event=as.factor(re$event)
library(pROC)
library(ggplot2)
m <- roc(pd$event, re$prob_min)



g <- ggroc(m,legacy.axes = T,size = 1,color = "#2fa1dd")
auc(m); ci.auc(m)
g + theme_minimal() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               colour = "grey", linetype = "dashed")+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of Train = ",format(round(as.numeric(auc(m)),2),nsmall = 2)),color = "#2fa1dd")
# 
table(pd_yi$status)
logctic_gene <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
exp_yi <- exp_yi[logctic_gene,rownames(pd_yi)]

lasso.prob.1 <- predict(log2, newdata = as.data.frame(t(exp_yi)))
lasso.prob.1<-exp(lasso.prob.1)/(1+exp(lasso.prob.1))
re.1=cbind(pd_yi$event ,lasso.prob.1)
head(re.1)
re.1=as.data.frame(re.1)
colnames(re.1)=c('event','prob_min')
re.1$event=as.factor(re.1$event)
m1 <- roc(pd_yi$event, re.1$prob_min)
auc(m1)
g1 <- ggroc(list(train = m,test = m1),legacy.axes = T,size = 1)

g1 + theme_minimal() +
  scale_color_manual(values = c("#f87669", "#2fa1dd"))+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               colour = "grey", linetype = "dashed")+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of Train = ",format(0.83,nsmall = 2)),color = "#f87669")+
  annotate("text",x = .75, y = .15,
           label = paste("AUC of Test = ",format(round(as.numeric(auc(m1)),2),nsmall = 2)),color = "#2fa1dd")

tiff(filename = "D:// /CAD-2//plot//6.4logistic.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//6.4logistic.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()



#####################################################################
### step7  immune
g <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
dat=pd_jia
dat$event <- as.numeric(ifelse(str_detect(dat$status,"MID+"),
                              "1",
                              "0"))
dat <- dat[,c("ID", "event")]
i=t(exprSet[g,])
i=i[rownames(pd_jia),]
dat=cbind(dat,i)
colnames(dat)
dat4 = dat
fp.3 <- predict(log2,type="response",newdata = dat4)##“link”, “response”, “terms”
medianTrainRisk=median(fp.3)
risk=as.vector(ifelse(fp.3>medianTrainRisk,"high","low"))
group_model <- ifelse(risk == "high", "high", "low")
group_model = factor(group_model,levels = c("high", "low"))
table(group_model)
save(group_model,dat4,file = "step7.1_immuGroup.Rdata")


#
rm(list = ls())
load("step1_data.Rdata");load("step7.1_immuGroup.Rdata")
library(GSVA)
library(tidyverse)
library(tinyarray)
geneSet <- read.csv("cellmarker.csv",header = F,sep = ",")
class(geneSet)
geneSet <- geneSet %>%
  column_to_rownames("V1")%>%t()
a <- geneSet
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
#save(l,file = "./gene_set.Rdata")
dat <- as.matrix(exp_jia)
dat <- dat[ ,colnames(dat) %in% pd_jia$ID]
ssgsea <- gsva(dat, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

draw_boxplot(ssgsea,group_model,color = c("#f87669","#2fa1dd"))
draw_pca(ssgsea,group_model,color = c( "#f87669","#2fa1dd"))
pca.plot.3D = draw_pca(ssgsea,group_model, style = "3D" );pca.plot.3D
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}
apply(ssgsea.1[,1:6], 2, range)
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#f87669","#2fa1dd")),
                       labels = c("high", "low"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

k <- Heatmap(t(scale(t(ssgsea.1))),name = " ",
             col = col_fun,
             top_annotation = top_annotation,
             column_split = group_model,
             show_heatmap_legend = T,
             column_gap = unit(0, "mm"),
             show_column_dend = F,
             border = F,
             show_column_names = F,
             show_row_names = T,
             column_title = NULL)
k
library(ggplotify)
k <- as.ggplot(k)#

#
library(Hmisc)
g <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
identical(colnames(ssgsea),colnames(dat))
nc = t(rbind(ssgsea,dat[g,]))
m = rcorr(nc)$r[1:nrow(ssgsea),(ncol(nc)-length(g)+1):ncol(nc)]
p = rcorr(nc)$P[1:nrow(ssgsea),(ncol(nc)-length(g)+1):ncol(nc)]
p[1:4,1:4]
library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
library(pheatmap)
pheatmap(t(m),
         display_numbers =t(tmp),
         angle_col ="45",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)

tiff(filename = "D:// /CAD-2//plot//7.1.1ssgsea_box.tiff", res = 300, width = 5000, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//7.1.1ssgsea_box.pdf",width = 18, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//7.1.2ssgsea_heat.tiff", res = 300, width = 3500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//7.1.2ssgsea_heat.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//7.1.3ssgsea_gene.tiff", res = 300, width = 3500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//7.1.3ssgsea_gene.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()

# estimate
library(estimate)
library(stringr)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="affymetrix")   ## platform
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='CAD'
scores=estimate(dat,pro)
head(scores)
library(tibble)
scores = rownames_to_column(as.data.frame(scores),"id")
scores$id = str_replace_all(scores$id,"\\.","-")
scores$Group = group_model
library(ggplot2)
library(ggpubr)
b1 = ggplot(dat = scores,aes(Group,StromalScore))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()
b2 = ggplot(dat = scores,aes(Group,ImmuneScore ))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()
b3 = ggplot(dat = scores,aes(Group,ESTIMATEScore ))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()
b4 = ggplot(dat = scores,aes(Group,TumorPurity ))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()
b1+b2+b3+plot_layout(nrow=1)
tiff(filename = "D:// /CAD-2//plot//7.2.1estimate.tiff", res = 300, width = 3500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//7.2.1estimate.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()


## cibersort
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
source('Cibersort.R')
LM22.file <- "LM22.txt"
#write.csv(exp_jia,file = "Train_exp.csv")
exp.file <- "Train_exp.csv"
tme <- CIBERSORT(LM22.file,exp.file,perm = 1000,QN = TRUE)
TCGA_TME.results <- tme
group_list <- group_model
table(group_list) 

## 3. 
# 3.1 
TME_data <- as.data.frame(TCGA_TME.results[,1:22])

TME_data$group <- group_list
TME_data$sample <- row.names(TME_data)

# 2.2 
TME_New = melt(TME_data)
## Using group, sample as id variables
colnames(TME_New)=c("Group","Sample","Celltype","Composition")  #
head(TME_New)
# 3.3 
plot_order = TME_New[TME_New$Group=="high",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)
## `summarise()` ungrouping output (override with `.groups` argument)
TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)

if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)

box_TME
save(group_model, ssgsea, scores, tme, file = "immu.Rdata")
tiff(filename = "D:// /CAD-2//plot//7.3Cibersort.tiff", res = 300, width = 5000, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//7.3Cibersort.pdf",width = 15, height = 8, onefile = FALSE)
dev.off()


#################################################################
### step8 
rm(list = ls())
load("step1_data.Rdata");load("immu.Rdata")
#immu_exp <- t(exp_jia);immu_exp <- as.data.frame(immu_exp)
#immu_exp$group <- group_model;immu_exp$group <- ifelse(immu_exp$group == "high", "Group1", "Group2")
#immu_exp <- t(immu_exp);immu_exp <- as.data.frame(immu_exp)
#immu_exp$Symbol <- rownames(immu_exp)
#cols <- colnames(immu_exp)
#new_cols <- c(cols[97], cols[1:(length(cols) - 1)])
#immu_exp <- immu_exp[, new_cols]
#write.table(immu_exp,file = "immu_exp.2.txt",quote = FALSE,sep = '\t',row.names = FALSE)

#a <- read.delim("immu_exp.2.txt")
#View(a);dim(a)


tide <- read.csv("immu_tide.csv")
tide$id = tide$Patient; o <- scores[,-c(2,3,4)]
tide <- merge(tide, o, by = "id");tide <- tide[,-2]
b1 = ggplot(dat = tide,aes(Group,TIDE))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()+
  labs(title = "TIDE")
b1


# immu
#immuCell.1 <- read.delim("467996cf-f91b-41c3-8f02-f5d69af87344.txt")
#immuCell.2 <- read.delim("467996cf-f91b-41c3-8f02-f5d69af87344_group.txt")
immuCell.3 <- read.delim("ImmuCellAI_icb_result.txt")
immuCell.3$id = immuCell.3$X; o <- scores[,-c(2,3,4,5)]
immuCell.3 <- merge(immuCell.3, o, by = "id");immuCell.3 <- immuCell.3[,-2]
immuCell.3 <- immuCell.3[,c("Score","Group")]
colnames(immuCell.3) <- c("treat_score", "Group")
immuCell_result = ggplot(dat = immuCell.3,aes(Group,treat_score))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()+
  labs(title = "ImmuCellAI")
immuCell_result
tiff(filename = "D:// /CAD-2//plot//8.1immucellai.tiff", res = 300, width = 1500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//8.1immucellai.pdf",width = 8,height = 8,onefile = FALSE)
dev.off()


#easier
library(oncoPredict);library(easier)
RNA_tpm <- exp_jia
RNA_counts <- exp_jia
genes_info <- easier:::reannotate_genes(cur_genes = rownames(RNA_tpm))

## symbol
non_na <- !is.na(genes_info$new_names)
RNA_tpm <- RNA_tpm[non_na, ]
genes_info <- genes_info[non_na, ]

## 
RNA_tpm <- RNA_tpm[-which(genes_info$new_names == "entry withdrawn"), ]
genes_info <- genes_info[-which(genes_info$new_names == "entry withdrawn"), ]

## 
newnames_dup <- unique(genes_info$new_names[duplicated(genes_info$new_names)])
newnames_dup_ind <- do.call(c, lapply(newnames_dup, function(X) which(genes_info$new_names == X)))
newnames_dup <- genes_info$new_names[newnames_dup_ind]

## 
tmp <- RNA_tpm[genes_info$old_names[genes_info$new_names %in% newnames_dup],]

## 
RNA_tpm <- RNA_tpm[-which(rownames(RNA_tpm) %in% rownames(tmp)),]
RNA_counts <- RNA_counts[rownames(RNA_tpm),]

## 
dup_genes <- genes_info$new_names[which(genes_info$new_names %in% newnames_dup)]
names(dup_genes) <- rownames(tmp)
if (anyDuplicated(newnames_dup)){
  tmp2 <- stats::aggregate(tmp, by = list(dup_genes), FUN = "mean")
  rownames(tmp2) <- tmp2$Group.1
  tmp2$Group.1 <- NULL
  # 
  RNA_tpm <- rbind(RNA_tpm, tmp2)
  RNA_counts <- rbind(RNA_counts, tmp2)
}

hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response_scores <- compute_scores_immune_response(RNA_tpm = RNA_tpm, 
                                                         selected_scores = hallmarks_of_immune_response)
pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts, remove_sig_genes_immune_response = TRUE)
cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm, cancer_type = "UCEC")
ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, cancer_type = "UCEC")
tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
predictions <- predict_immune_response(pathways = pathway_activities,
                                       immunecells = cell_fractions,
                                       tfs = tf_activities,
                                       lrpairs = lrpair_weights,
                                       ccpairs = ccpair_scores,
                                       cancer_type = "UCEC", 
                                       verbose = TRUE)
easier_derived_scores <- retrieve_easier_score(predictions_immune_response = predictions,
                                               easier_with_TMB = c("weighted_average", 
                                                                   "penalized_score"),
                                               weight_penalty = 0.5)
easier_derived_scores$id = rownames(easier_derived_scores); o <- scores[,-c(2,3,4)]
easier_derived_scores <- merge(easier_derived_scores, o, by = "id")
easier_polt = ggplot(dat = easier_derived_scores,aes(Group,easier_score))+
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()+
  labs(title = "Easier")
easier_polt
tiff(filename = "D:// /CAD-2//plot//8.2easier.tiff", res = 300, width = 1500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//8.2easier.pdf",width = 8,height = 8,onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//8.3treat.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
immuCell_result + easier_polt
dev.off()
pdf(file = "D:// /CAD-2//plot//8.3treat.pdf",width = 8,height = 8,onefile = FALSE)
immuCell_result + easier_polt
dev.off()



#############################################################################
### step9 

rm(list = ls())
library(clusterProfiler)
library(ggthemes)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
load("step5_CO_deg.Rdata");load("step1_data.Rdata");load("step2_deg.Rdata")
co_gene <- c(up,down)
#write.csv(co_gene, file = "Co_gene.txt")
# 1.GO 
DEG_up <- DEG[up, ]
DEG_down <- DEG[down, ]

DEG_up$symbol <- rownames(DEG_up)
DEG_down$symbol <- rownames(DEG_down)
s2e.1 <- bitr(DEG_up$symbol, 
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
dim(DEG_up);dim(DEG_down)
DEG_up <- inner_join(DEG_up,s2e.1,by=c("symbol"="SYMBOL"))
dim(DEG_up)
length(unique(DEG_up$symbol))
s2e.2 <- bitr(DEG_down$symbol, 
              fromType = "SYMBOL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
dim(DEG_down)
DEG_down <- inner_join(DEG_down,s2e.2,by=c("symbol"="SYMBOL"))
dim(DEG_down)
length(unique(DEG_down$symbol))

gene_up = DEG_up$ENTREZID
gene_down = DEG_down$ENTREZID
gene_diff = c(gene_up,gene_down)


#(2)
if(!file.exists("CAD_GO.Rdata")){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb= org.Hs.eg.db,
                  ont = "ALL",
                  readable = TRUE)
  ego_BP <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     readable = TRUE)
  ego_MF <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "MF",
                     readable = TRUE)
  ego_CC <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "CC",
                     readable = TRUE)
  #ont "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego, ego_BP, ego_MF, ego_CC, file = "CAD_GO.Rdata")
}
load("CAD_GO.Rdata")

#(3)
barplot(ego)
dotplot(ego)
dotplot(ego, split = "ONTOLOGY", font.size = 10, 
        showCategory = 5) + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
#geneList 
geneList = c(DEG_down$logFC, DEG_up$logFC) 
names(geneList)= c(DEG_down$ENTREZID, DEG_up$ENTREZID)

#(3)
#Gene-Concept Network
cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
cnetplot(ego, showCategory = 5,foldChange=geneList, circular = TRUE, colorEdge = TRUE)


Biobase::package.version("enrichplot")
emapplot(pairwise_termsim(ego)) #
#(4)
#goplot(ego)
goplot(ego_BP)
goplot(ego_CC)
goplot(ego_MF)
#(5)Heatmap-like functional classification
heatplot(ego,foldChange = geneList,showCategory = 8)
# 2.KEGG pathway analysis----

gene_up = DEG_up[ ,'ENTREZID'] 
gene_down = DEG_down[ ,'ENTREZID'] 
gene_diff = c(gene_up,gene_down)
#（2）
if(!file.exists("CAD_KEGG.Rdata")){
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa')
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa')
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa')
  save(kk.diff,kk.down,kk.up,file = "CAD_KEGG.Rdata")
}
load("CAD_KEGG.Rdata")

#(3)https://mp.weixin.qq.com/s/NglawJgVgrMJ0QfD-YRBQg
table(kk.diff@result$p.adjust<0.05)
table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)
up$kk@result = mutate(up$kk@result, change = "up")
down$kk@result = mutate(down$kk@result, change = "down")
kk = rbind(up$kk@result[1:8, ], down$kk@result[1:10, ])
ud_enrich = function(df) {
  df$pl = ifelse(df$change == "up", -log10(df$p.adjust), 
                 log10(df$p.adjust))
  df = arrange(df, change, pl)
  df$Description = factor(df$Description, levels = unique(df$Description), 
                          ordered = TRUE)
  tmp = with(df, labeling::extended(range(pl)[1], range(pl)[2], 
                                    m = 5))
  lm = tmp[c(1, length(tmp))]
  lm = c(floor(min(df$pl)), ceiling(max(df$pl)))
  ggplot(df, aes(x = Description, y = pl)) + geom_bar(stat = "identity", 
                                                      aes(fill = change), width = 0.7) + scale_fill_manual(values = color) + 
    coord_flip() + theme_light() + ylim(lm) + scale_x_discrete(labels = function(x) str_wrap(x, 
                                                                                             width = 30)) + scale_y_continuous(breaks = tmp, labels = abs(tmp)) + 
    theme(panel.border = element_blank())
}
color = c("#2fa1dd", "#f87669")
kp = ud_enrich(kk)
kp
#g_kegg +scale_y_continuous(labels = c(2,0,2,4,6))

## GSEA
library(clusterProfiler);library(fgsea)
alldiff <- DEG[order(DEG$logFC,decreasing = T),]
id <- alldiff$logFC
names(id) <- rownames(alldiff)
gmtfile <- "h.all.v7.5.1.symbols.gmt"
hallmark <- read.gmt(gmtfile)
hallmark$term <- gsub('HALLMARK_','',hallmark$term)
hallmark.list <- hallmark %>% split(.$term) %>% lapply( "[[", 2)
fgseaRes <- fgsea(pathways = hallmark.list,
                  stats = id,
                  minSize=1,
                  maxSize=10000)
sig <- fgseaRes[fgseaRes$pval<0.05,]
sig <- sig[order(sig$NES,decreasing = T),]
plotGseaTable(hallmark.list[sig$pathway],id, fgseaRes,gseaParam = 0.5)
# 2.
deg <- DEG[co_gene, ]
deg$symbol <- rownames(deg)
p = deg[ , c("symbol", "logFC")]
head(p)
write.table(p,
            file = "deg.txt",
            sep = "\t",
            quote = F,
            row.names = F)
# 

tiff(filename = "D:// /CAD-2//plot//5.1GO(bar).tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.1GO(bar).pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.1GO(dot).tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.1GO(dot).pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.1GO(dot_3).tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.1GO(dot_3).pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.2Concept.tiff", res = 300, width = 5000, height = 3500,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.2Concept.pdf",width = 15, height = 12, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.2Concept_cir.tiff", res = 300, width = 3500, height = 2500,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.2Concept_cir.pdf",width = 15, height = 12, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.3GO_pathway_BP.tiff", res = 300, width = 5000, height = 3500,compression = "lzw")
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.3GO_pathway_CC.tiff", res = 300, width = 4000, height = 2800,compression = "lzw")
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.3GO_pathway_MF.tiff", res = 300, width = 3500, height = 2500,compression = "lzw")
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.4Pathway_gene(heat).tiff", res = 300, width = 5500, height = 2800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.4Pathway_gene(heat).pdf",width = 18, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.5kegg.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.5kegg.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()
tiff(filename = "D:// /CAD-2//plot//5.6gsea.tiff", res = 300, width = 2500, height = 1800,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//5.6gsea.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()


#####################################################################
### step10 
rm(list = ls())
load("step1_data.Rdata");load(("immu.Rdata"))
exprSet <- as.matrix(exp_jia)
library(ggplot2);library(ggpubr)
o <- scores[,-c(2,3,4)]
dat.exp.1 = data.frame(gene = exprSet["HIST1H4E",],group_risk = o$Group)
colnames(dat.exp.1) <- c("HIST1H4E","risk_group")
p1 = ggplot(dat.exp.1,aes(risk_group, HIST1H4E)) +
  geom_boxplot(aes(fill = risk_group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.2 = data.frame(gene = exprSet["IL6ST",],group_risk = o$Group)
colnames(dat.exp.2) <- c("IL6ST","risk_group")
p2 = ggplot(dat.exp.2,aes(risk_group, IL6ST)) +
  geom_boxplot(aes(fill = risk_group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.3 = data.frame(gene = exprSet["RN7SKP45",],group_risk = o$Group)
colnames(dat.exp.3) <- c("RN7SKP45","risk_group")
p3 = ggplot(dat.exp.3,aes(risk_group, RN7SKP45)) +
  geom_boxplot(aes(fill = risk_group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.4 = data.frame(gene = exprSet["LST1",],group_risk = o$Group)
colnames(dat.exp.4) <- c("LST1","risk_group")
p4 = ggplot(dat.exp.4,aes(risk_group, LST1)) +
  geom_boxplot(aes(fill = risk_group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.5 = data.frame(gene = exprSet["SNORD50B",],group_risk = o$Group)
colnames(dat.exp.5) <- c("SNORD50B","risk_group")
p5 = ggplot(dat.exp.5,aes(risk_group, SNORD50B)) +
  geom_boxplot(aes(fill = risk_group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()
library(patchwork)
p1+p2+p3+p4+p5+plot_layout(nrow=1)
tiff(filename = "D:// /CAD-2//plot//9.exp.tiff", res = 300, width = 5000, height = 1500,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//9.exp.pdf",width = 18, height = 5, onefile = FALSE)
dev.off()

# GSE120774
gse_number = "GSE120774"
eSet <- getGEO(gse_number, destdir = '.', getGPL = F)
class(eSet)
length(eSet)
eSet = eSet[[1]]
pd_eat <- pData(eSet)
exp_eat <- exprs(eSet)
boxplot(exp_eat)
pd_eat <- pd_eat[,c("title", "geo_accession", "disease:ch1", "gender:ch1")]
e <- pd_eat$title[str_detect(pd_eat$title,"Epicardial fat")]
pd_eat <- filter(pd_eat, pd_eat$title %in% e)
pd_eat <- pd_eat[,-1]
Group_eat <- ifelse(pd_eat$`disease:ch1` == "None", "Control", "CAD")
Group_eat = factor(Group_eat,levels = c("Control","CAD"))
table(Group_eat)
pd_eat$event <- ifelse(Group_eat == "Control", 0, 1)
exp_eat <- exp_eat[,pd_eat$geo_accession]
library(tinyarray)
ids <- AnnoProbe::idmap('GPL6244')
exp_eat <- as.data.frame(exp_eat)
exp_eat <- mutate(exp_eat,probe_id=rownames(exp_eat))
ids = ids[!duplicated(ids$symbol),]
exp_eat <- inner_join(exp_eat,ids,by="probe_id")
nrow(exp_eat)
rownames(exp_eat) <- exp_eat$symbol
exp_eat <- select(exp_eat, -c("probe_id", "symbol"))
dim(exp_eat);dim(pd_eat);table(Group_eat)
save(exp_eat, pd_eat, Group_eat, file = "eat_test.Rdata")

exp_eat <- as.matrix(exp_eat)
dat.exp.1 = data.frame(gene = exp_eat["HIST1H4E",],group = Group_eat)
colnames(dat.exp.1) <- c("HIST1H4E","Group")
Q1 = ggplot(dat.exp.1,aes(Group, HIST1H4E)) +
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.2 = data.frame(gene = exp_eat["IL6ST",],group = Group_eat)
colnames(dat.exp.2) <- c("IL6ST","Group")
Q2 = ggplot(dat.exp.2,aes(Group, IL6ST)) +
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.3 = data.frame(gene = exp_eat["RN7SKP45",],group = Group_eat)
colnames(dat.exp.3) <- c("RN7SKP45","Group")
Q3 = ggplot(dat.exp.3,aes(Group, RN7SKP45)) +
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.4 = data.frame(gene = exp_eat["LST1",],group = Group_eat)
colnames(dat.exp.4) <- c("LST1","Group")
Q4 = ggplot(dat.exp.4,aes(Group, LST1)) +
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()

dat.exp.5 = data.frame(gene = exp_eat["SNORD50B",],group = Group_eat)
colnames(dat.exp.5) <- c("SNORD50B","Group")
Q5 = ggplot(dat.exp.5,aes(Group, SNORD50B)) +
  geom_boxplot(aes(fill = Group))+
  scale_fill_manual(values = c("#f87669","#2fa1dd"))+
  stat_compare_means()+
  theme_classic()
library(patchwork)
Q1+Q2+Q4+Q5

### 
g <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
library(corrplot);library(pheatmap)
M = cor(t(exprSet[g,]))
pheatmap(M)
library(paletteer)
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10)
corrplot(M, type="upper",
         method="pie",
         order="hclust", 
         col=my_color,
         tl.col="black", 
         tl.srt=45)
library(cowplot)
cor_plot <- recordPlot() 

a = exprSet[g,]
M = cor(t(a))
library(circlize)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
mat = M
df = mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(names_to  = "to",
               values_to  = "value",
               cols = 2:(nrow(mat)+1) #
  ) 
df = df[df$from != df$to,]
df = df[!duplicated(df$value),]

#
library(RColorBrewer)
col_sample = c(brewer.pal(n = 12,name = "Set3")
               #brewer.pal(n = 8,name = "Set2")
)
col_sample
border_color <- sample(col_sample,nrow(mat)) 
#border_color = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3") 

#
col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5)
# range(mat)
# col <- colorRamp2(range(mat), c("#FA0A0A", "#FFFFFF")) 

# ComplexHeatmapt
lgd = Legend(col_fun = col, title = "foo")
ComplexHeatmap:::width(lgd)
ComplexHeatmap:::height(lgd)

pushViewport(viewport(width = 0.9, height = 0.9))
grid.rect()  # border
draw(lgd, x = unit(5, "cm"), y = unit(3, "cm"), just = c("left", "bottom"))
# draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
# draw(lgd, x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
popViewport()


chordDiagram(x = df, 
             directional = 0, #
             grid.col = border_color, #
             col = col, #
             annotationTrack = c('grid', 'name', 'axis'), #
             annotationTrackHeight = c(0.05, 0.1),#
             transparency = 0.25#
             )

tiff(filename = "D:// /CAD-2//plot//9.2gene_cor.tiff", res = 300, width = 3000, height = 2000,compression = "lzw")
cor_plot
dev.off()
pdf(file = "D:// /CAD-2//plot//9.2gene_cor.pdf",width = 12, height = 8, onefile = FALSE)
cor_plot
dev.off()
tiff(filename = "D:// /CAD-2//plot//9.3gene_cor_cir.tiff", res = 300, width = 3000, height = 2000,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//9.3gene_cor_cir.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()

#####################################################################
### step11 
rm(list = ls())
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
#输入文件
load("step1_data.Rdata")
g <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
sgene="SNORD50B"       #
data=as.matrix(exp_jia)
#
group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))
#
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
logFC_t=1
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') ))
table(deg$g)
deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))
geneList = data_all_sort$logFC #
names(geneList) <- data_all_sort$ENTREZID #
head(geneList)
#
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
#write.table(af,file=paste0("2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)
#
num=5
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
#
num=5
SNORD50B <- gseaplot2(kk2,
          title = "SNORD50B",
          geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
#
gseaplot2(kk2,
          title = "Th1 and Th2 cell differentiation",  #
          "hsa00260", #
          color="red", #
          base_size = 20, #
          subplots = 1:3, 
          pvalue_table = T) # 
tiff(filename = "D:// /CAD-2//plot//10.5SNORD50B.tiff", res = 300, width = 3000, height = 2000,compression = "lzw")
dev.off()
pdf(file = "D:// /CAD-2//plot//10.5SNORD50B.pdf",width = 12, height = 8, onefile = FALSE)
dev.off()


######################################################################
### step12 ceRNA
rm(list = ls())
library(multiMiR)
load("step1_data.Rdata")
x <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
gene2mir <- get_multimir(org     = 'hsa',
                         target  = x,
                         table   = 'predicted',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
mit = gene2mir@data[gene2mir@data$database=="microcosm",];dim(mit)
table(mit$target_symbol)
#mit = mit[stringr::str_detect(mit$experiment,
#                              "Degradome sequencing"),];dim(mit)
miRNAs = unique(mit$mature_mirna_id)
# 
starbase = data.table::fread("starBaseV3_hg19_CLIP-seq_lncRNA_all.txt");dim(starbase)
load("anno.Rdata")
lnc_anno$gene_id = stringr::str_remove(lnc_anno$gene_id,"\\.\\d")
p1 = starbase$geneName %in% lnc_anno$gene_name;table(p1)
p2 = starbase$geneID %in% lnc_anno$gene_id;table(p2)
starbase = starbase[p2,];dim(starbase)
lnc_mi = starbase[starbase$miRNAname %in% miRNAs,]
colnames(lnc_mi)
p2 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>4;table(p2)
#degradome evidence 
p3 = lnc_mi$pancancerNum >10  &lnc_mi$clipExpNum>4 & lnc_mi$degraExpNum >0;table(p3)
lnc_mi$geneName[p2]
lnc_mi$geneName[p3]
ez1 = mit[,3:4]
ez2 = lnc_mi[p2,c(2,4)]
colnames(ez1) = colnames(ez2)

library(dplyr)
ez = rbind(ez1,ez2);dim(ez)
ez = distinct(ez,miRNAname,geneName);dim(ez)
tp = data.frame(nodes = c(ez1$miRNAname,
                          ez2$miRNAname,
                          ez1$geneName,
                          ez2$geneName),
                type = rep(c("mi","pc","lnc"),
                           times = c(nrow(ez1)+nrow(ez2),
                                     nrow(ez1),
                                     nrow(ez2))
                )
)
dim(tp);head(tp)
tp = distinct(tp,nodes,.keep_all = T)
table(tp$type)
write.table(ez,file = "ez.txt",quote = F,sep = "\t",row.names = F)
write.table(tp,file = "tp.txt",quote = F,sep = "\t",row.names = F)


############################################################################
### step13 
rm(list = ls())
load("step1_data.Rdata");rm(exp_yi);rm(pd_yi);rm(pd_jia);rm(Group_yi)
exp_jia <- exp_jia[, Group_jia == "CAD"]
g <- c("HIST1H4E","IL6ST","RN7SKP45","LST1","SNORD50B")
immunoinhibitor <- c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R",
                     "CTLA4","HAVCR2","IDO1","IL10","IL10RB","KDR","KIR2DL1",
                     "KIR2DL3","LAG3","LGALS9","PDCD1","PDCD1LG2","PVRL2",
                     "TGFB1","TGFBR1","TIGIT","VTCN1")
immunostimulator <- c("C10orf54","CD27","CD276","CD28","CD40","CD40LG","CD48",
                      "CD70","CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2",
                      "ICOS","ICOSLG","IL2RA","IL6","IL6R","KLRC1","KLRK1","LTA",
                      "MICB","NT5E","PVR","RAET1E","TMEM173","TMIGD2","TNFRSF13B",
                      "TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25",
                      "TNFRSF4","TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14",
                      "TNFSF15","TNFSF18","TNFSF4","TNFSF9","ULBP1")
MHC <- c("B2M","HLA-A","HLA-B","HLA-C","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1",
         "HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-E",
         "HLA-F","HLA-G","TAP1","TAP2","TAPBP")
exprSet <- exp_jia
library(corrplot) 
library(vegan) 
library(ggcor)
library(ggplot2)
library(tidyverse)
library(ambient)
### 
varespec <- t(exprSet[g,])
varechem <- t(exprSet[immunoinhibitor,])
mantel <- mantel_test(varespec, varechem,
                      spec.select = list(HIST1H4E = 1,
                                         IL6ST = 2,
                                         RN7SKP45 = 3,
                                         LST1 = 4,
                                         SNORD50B = 5)) %>%  
  mutate(rd = cut(r, breaks = c(-Inf, 0.4, 0.6, Inf), 
                  labels = c("< 0.4", "0.4 - 0.6", ">= 0.6")),# 
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#
a <- quickcor(varechem, type = "upper") +#
  scale_fill_gradient2(low="#2fa1dd", high="#f87669", mid="white")+
  geom_star()+ #
  anno_link(aes(colour = pd, size = rd), data = mantel) +#
  scale_size_manual(values = c(0.5, 1, 2))+ 
  scale_fill_gradientn(values = seq(0,1,0.2),
                       colors = c("#1543B3", "#2fa1dd", "#f5f4f4", "#D04B35", "#920A07")) +
  scale_color_manual(values = c( "#f87669", "#82D035","#2fa1dd"))+
  guides(size = guide_legend(title = "Mantel's r",#
                             order = 2), 
         colour = guide_legend(title = "Mantel's p",  
                               order = 3), 
         fill = guide_colorbar(title = "Pearson's r", order = 4)) 
tiff(filename = "D:// /CAD-2//plot//13.1inhibitor.tiff", res = 300, width = 4000, height = 3000,compression = "lzw")
a
dev.off()
pdf(file = "D:// /CAD-2//plot//13.1inhibitor.pdf",width = 12, height = 8, onefile = FALSE)
a
dev.off()

varespec <- t(exprSet[g,])
varechem.1 <- t(exprSet[immunostimulator,])
mantel <- mantel_test(varespec, varechem.1,
                      spec.select = list(HIST1H4E = 1,
                                         IL6ST = 2,
                                         RN7SKP45 = 3,
                                         LST1 = 4,
                                         SNORD50B = 5)) %>%  
  mutate(rd = cut(r, breaks = c(-Inf, 0.4, 0.6, Inf), 
                  labels = c("< 0.4", "0.4 - 0.6", ">= 0.6")),#
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#
b <- quickcor(varechem.1, type = "upper") +#
  scale_fill_gradient2(low="#2fa1dd", high="#f87669", mid="white")+
  geom_star()+ #
  anno_link(aes(colour = pd, size = rd), data = mantel) +#
  scale_size_manual(values = c(0.5, 1, 2))+ 
  scale_fill_gradientn(values = seq(0,1,0.2),
                       colors = c("#1543B3", "#2fa1dd", "#f5f4f4", "#D04B35", "#920A07")) +
  scale_color_manual(values = c( "#f87669", "#82D035","#2fa1dd"))+
  guides(size = guide_legend(title = "Mantel's r",#
                             order = 2), 
         colour = guide_legend(title = "Mantel's p",  
                               order = 3), 
         fill = guide_colorbar(title = "Pearson's r", order = 4)) 
tiff(filename = "D:// /CAD-2//plot//13.2stimulator.tiff", res = 300, width = 4000, height = 3000,compression = "lzw")
b
dev.off()
pdf(file = "D:// /CAD-2//plot//13.2stimulator.pdf",width = 15, height = 10, onefile = FALSE)
b
dev.off()

### MHC
varespec <- t(exprSet[g,])
varechem.2 <- t(exprSet[MHC,])
mantel <- mantel_test(varespec, varechem.2,
                      spec.select = list(HIST1H4E = 1,
                                         IL6ST = 2,
                                         RN7SKP45 = 3,
                                         LST1 = 4,
                                         SNORD50B = 5)) %>%  
  mutate(rd = cut(r, breaks = c(-Inf, 0.4, 0.6, Inf), 
                  labels = c("< 0.4", "0.4 - 0.6", ">= 0.6")),#
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf), 
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))#
c <- quickcor(varechem.2, type = "upper") +#
  scale_fill_gradient2(low="#2fa1dd", high="#f87669", mid="white")+
  geom_star()+ #
  anno_link(aes(colour = pd, size = rd), data = mantel) +#
  scale_size_manual(values = c(0.5, 1, 2))+ 
  scale_fill_gradientn(values = seq(0,1,0.2),
                       colors = c("#1543B3", "#2fa1dd", "#f5f4f4", "#D04B35", "#920A07")) +
  scale_color_manual(values = c( "#f87669", "#82D035","#2fa1dd"))+
  guides(size = guide_legend(title = "Mantel's r",#
                             order = 2), 
         colour = guide_legend(title = "Mantel's p",  
                               order = 3), 
         fill = guide_colorbar(title = "Pearson's r", order = 4)) 
tiff(filename = "D:// /CAD-2//plot//13.3MHC.tiff", res = 300, width = 4000, height = 3000,compression = "lzw")
c
dev.off()
pdf(file = "D:// /CAD-2//plot//13.3MHC.pdf",width = 12, height = 8, onefile = FALSE)
c
dev.off()


data_exp <- as.data.frame(t(exprSet))
library(ggstatsplot)
p1 <- ggscatterstats(
  data = data_exp,
  x = BTLA,
  y = HIST1H4E,
  point.args = list(color = "#82D035",size = 2),
  smooth.line.args = list(color = "#10100A",method = "lm"),
  xlab = "BTLA",
  ylab = "HIST1H4E",
  title = "The correlation between BTLA and HISTIH4E",
  messages = FALSE,
  xsidehistogram.args = list(fill = "#2fa1dd",color = "#10100A"),
  ysidehistogram.args = list(fill = "#f87669",color = "#10100A"),
  marginal.type = "density"
)
p1

p2 <- ggscatterstats(
  data = data_exp,
  x = TMEM173,
  y = LST1,
  point.args = list(color = "#82D035",size = 2),
  smooth.line.args = list(color = "#10100A",method = "lm"),
  xlab = "TMEM173",
  ylab = "LST1",
  title = "The correlation between TMEM173 and LST1",
  messages = FALSE,
  xsidehistogram.args = list(fill = "#2fa1dd",color = "#10100A"),
  ysidehistogram.args = list(fill = "#f87669",color = "#10100A"),
  marginal.type = "density"
)
p2

p3 <- ggscatterstats(
  data = data_exp,
  x = TAPBP,
  y = LST1,
  point.args = list(color = "#82D035",size = 2),
  smooth.line.args = list(color = "#10100A",method = "lm"),
  xlab = "TAPBP",
  ylab = "LST1",
  title = "The correlation between TAPBP and LST1",
  messages = FALSE,
  xsidehistogram.args = list(fill = "#2fa1dd",color = "#10100A"),
  ysidehistogram.args = list(fill = "#f87669",color = "#10100A"),
  marginal.type = "density"
)
p3

p4 <- ggscatterstats(
  data = data_exp,
  x = TMEM173,
  y = HIST1H4E,
  point.args = list(color = "#82D035",size = 2),
  smooth.line.args = list(color = "#10100A",method = "lm"),
  xlab = "TMEM173",
  ylab = "HIST1H4E",
  title = "The correlation between TMEM173 and HISTIH4E",
  messages = FALSE,
  xsidehistogram.args = list(fill = "#2fa1dd",color = "#10100A"),
  ysidehistogram.args = list(fill = "#f87669",color = "#10100A"),
  marginal.type = "density"
)
p4
p1 + p4 + p2 + p3
tiff(filename = "D:// /CAD-2//plot//13.4cor.tiff", res = 300, width = 4000, height = 3000,compression = "lzw")
p1 + p4 + p2 + p3
dev.off()
pdf(file = "D:// /CAD-2//plot//13.4cor.pdf",width = 12, height = 8, onefile = FALSE)
p1 + p4 + p2 + p3
dev.off()
