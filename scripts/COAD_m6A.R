if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('z:/projects/codes/mg_base.R')
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}
mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
                     quick=T,mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#
  #,observation=pbc[2,] #
  #
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #
  #              ,showP = T #
  #              ,droplines = F#
  #,colors = mg_colors[1:3] #
  #,rank="decreasing") #
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=as.data.frame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

getC_index=function(riskscore,os,status){
  inds=which(!is.na(riskscore)&!is.na(os)&!is.na(status))
  riskscore=riskscore[inds]
  os=os[inds]
  status=status[inds]
  c1=survcomp::concordance.index(x=riskscore, surv.time=os, surv.event=status,
                                 method="noether")
  #c2=concordance.index(x=riskscore[order(rnorm(length(riskscore)))], surv.time=os, surv.event=status,
  #                     method="noether")
  #p=min(cindex.comp(c1, c2)$p.value,cindex.comp(c2, c1)$p.value)  
  return(c1)
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2024)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}

#
genecode=read.delim('GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)

#Sidra-LUMC AC-ICAM队列  PMID：37202560#########
###
silu.data=read.delim('origin_datas/coad_silu_2022/data_mrna_seq_expression.txt',check.names = F)
silu.data=silu.data[!duplicated(silu.data$Hugo_Symbol),]
rownames(silu.data)=silu.data$Hugo_Symbol
silu.data=silu.data[,-1]
silu.data[1:5,1:5]
range(silu.data)
dim(silu.data)



####
silu.patient=read.delim('origin_datas/coad_silu_2022/data_clinical_patient.txt',check.names = F)
silu.patient=silu.patient[-c(1:4),]
head(silu.patient)

silu.sample=read.delim('origin_datas/coad_silu_2022/data_clinical_sample.txt',check.names = F)
silu.sample=silu.sample[-c(1:4),]
head(silu.sample)

silu.cli=merge(silu.patient,silu.sample,by.x='#Patient Id',by.y='#Patient Identifier')
head(silu.cli)
colnames(silu.cli)
silu.cli=data.frame(Samples=silu.cli$`Sample Identifier`,
                    Age=silu.cli$`Age At Diagnosis`,Gender=silu.cli$Sex,
                    OS=str_split_fixed(silu.cli$`Overall Survival Status`,':',2)[,1],
                    OS.time=silu.cli$`Overall Survival (Months)`,
                    DSS=str_split_fixed(silu.cli$`Disease-specific Survival Status`,':',2)[,1],
                    DSS.time=silu.cli$`Months of disease-specific survival`,
                    PFS=str_split_fixed(silu.cli$`Progression Free Status`,':',2)[,1],
                    PFS.time=silu.cli$`Progress Free Survival (Months)`,
                    silu.cli[,c(8:11,31,39,40,45)])
rownames(silu.cli)=silu.cli$Samples
silu.cli=crbind2DataFrame(silu.cli)
silu.cli$Path.Tumor.Stage=paste0('T',silu.cli$Path.Tumor.Stage)
silu.cli$Path.Nodes.Stage=paste0('N',silu.cli$Path.Nodes.Stage)
silu.cli$Path.Metastasis.Stage=paste0('M',silu.cli$Path.Metastasis.Stage)
table(silu.cli$Path.Metastasis.Stage)
silu.cli$Path.Metastasis.Stage[silu.cli$Path.Metastasis.Stage=='MNA']<-NA
silu.cli$AJCC.Path.Stage=paste0('Stage',silu.cli$AJCC.Path.Stage)
head(silu.cli)
silu.cli$DSS.time
###
silu.cli=silu.cli%>%drop_na(DSS.time)
fivenum(silu.cli$Age)
silu.cli$Age1=ifelse(silu.cli$Age>69,'>69','<=69')
dim(silu.cli)

silu.exp=silu.data[rownames(silu.data)%in%mrna_genecode$SYMBOL,silu.cli$Samples]
dim(silu.exp)# 17099   320
range(silu.exp)



#01. ########
dir.create('results/01.m6A_subtype')
m6A.regulators=read.xlsx('origin_datas/m6A_regulators_PMID_35646049.xlsx')
m6A.regulators=m6A.regulators$Regulator
length(m6A.regulators)
#22

m6A.regulators.cox=cox_batch(dat = silu.exp[m6A.regulators,silu.cli$Samples],
                             time = silu.cli$DSS.time,event = silu.cli$DSS)
m6A.regulators.cox=na.omit(m6A.regulators.cox)
table(m6A.regulators.cox$p.value<0.05)
m6A.regulators.cox[m6A.regulators.cox$p.value<0.05,]
write.csv(m6A.regulators.cox,'results/01.m6A_subtype/m6A.regulators.cox.csv')

pdf('results/01.m6A_subtype/m6A.regulators.forest.pdf',height = 10,width = 6,onefile = F)
bioForest(rt = m6A.regulators.cox,col=c('#80B1D3','#FDB462'))
dev.off()

library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[1]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[1]
consen_gene=rownames(m6A.regulators.cox[m6A.regulators.cox$p.value<0.05,])
length(consen_gene)
silu_consen_data=as.matrix(silu.exp[intersect(consen_gene,rownames(silu.exp)),])
silu_consen_data=t(scale(t(silu_consen_data),scale = T))   #11,3
silu_consen_data=as.matrix(silu_consen_data)
dim(silu_consen_data)
silu_clust_subtype <- ConsensusClusterPlus(silu_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "silu_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = F
                                           , seed = 123456)
k=3
subtype.cols=c("#1B9E77","#D95F02","#7570B3")
silu.subtype <- data.frame(Samples = names(silu_clust_subtype[[k]]$consensusClass),
                           Subtype=silu_clust_subtype[[k]]$consensusClass)
silu.subtype$Subtype=paste0('S',silu.subtype$Subtype)
table(silu.subtype$Subtype)
silu.subtype.cli=merge(silu.subtype,silu.cli,by='Samples')

silu.subtype.km=ggsurvplot(fit=survfit(Surv(DSS.time, DSS) ~ Subtype,data = silu.subtype.cli),
                           data=silu.subtype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                           ggtheme=custom_theme(),palette = subtype.cols,
                           linetype = c("solid", "dashed","strata")[1],
                           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                           legend.title = "Subtype",legend.labs=c('S1','S2','S3'))
silu.subtype.km=mg_merge_plot(silu.subtype.km$plot,silu.subtype.km$table,
                              ncol = 1,nrow = 2,heights = c(2.5,1),align = 'v')
silu.subtype.km
ggsave('results/01.m6A_subtype/silu.subtype.km.pdf',silu.subtype.km,height = 5,width = 5)

pdf('results/01.m6A_subtype/silu.subtype.expr.pdf',height = 5,width = 5)
Heatmap(as.matrix(t(scale(t(silu.exp[consen_gene,silu.subtype.cli$Samples])))),
        name = "Expr",
        column_split = silu.subtype.cli$Subtype,
        column_title_gp = gpar(fill =subtype.cols),
        cluster_rows = T, cluster_columns = F,
        cluster_row_slices = F, cluster_column_slices=T,
        show_row_dend = F, show_column_dend = F,
        show_row_names = T, show_column_names = F,
        col = circlize::colorRamp2(c(-3, 0, 3), c('#8787FFFF', 'white', '#FF8888FF')))
dev.off()

library(ggbiplot)
silu.subtype.pca <- prcomp(t(silu.exp[consen_gene,silu.subtype.cli$Samples]), scale=T)
pdf('results/01.m6A_subtype/silu.subtype.PCA.pdf',height = 5,width = 5)
ggbiplot(silu.subtype.pca, scale=1, groups = silu.subtype.cli$Subtype,
         ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values = subtype.cols) + 
  theme_bw() +
  theme(legend.direction = 'horizontal', legend.position = 'top', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab('PCA1') + ylab('PCA2')+ggtitle('silu')
dev.off()

my_mutibarplot=function(df,xlab='group',leg.title='',cols=brewer.pal(6,"Set2")){
  prop.pval=round(-log10(chisq.test(df)$p.value),2)
  df.prop=prop.table(df,margin=2)
  df.prop=reshape2::melt(df.prop)
  colnames(df.prop)<-c("type","group","Percentage")
  df.prop$Percentage<-round(df.prop$Percentage,digits=2)
  p=ggplot(df.prop,aes(x=group,y=Percentage,fill=type))+
    geom_bar(position = "fill",stat="identity")+
    scale_fill_manual(values = cols)+
    xlab(xlab)+labs(fill = leg.title,title = 'Chi-Squared Test',subtitle  =  paste0('-log10(pvalue) = ',prop.pval))+
    theme_bw()+theme(text=element_text(family = 'Times'),#legend.position = legend.position,
                     panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p
  return(p)
}

p1=my_mutibarplot(table(silu.subtype.cli$Path.Tumor.Stage,silu.subtype.cli$Subtype),
               xlab = 'Subtype',leg.title = 'Tumor Stage')
p2=my_mutibarplot(table(silu.subtype.cli$Path.Nodes.Stage,silu.subtype.cli$Subtype),
               xlab = 'Subtype',leg.title = 'Nodes Stage')
p3=my_mutibarplot(table(silu.subtype.cli$Path.Metastasis.Stage,silu.subtype.cli$Subtype),
               xlab = 'Subtype',leg.title = 'Metastasis Stage')
p4=my_mutibarplot(table(silu.subtype.cli$AJCC.Path.Stage,silu.subtype.cli$Subtype),
               xlab = 'Subtype',leg.title = 'AJCC Stage')

pdf('results/01.m6A_subtype/silu.subtype.clinical.pdf',height = 4,width = 16)
mg_merge_plot(p1,p2,p3,p4,ncol=4,labels = 'F')
dev.off()

#02.#############
dir.create('results/02.subtype_GSEA')
silu.subtype.use1=silu.subtype
silu.subtype.use1$Subtype1=ifelse(silu.subtype.use1$Subtype=='S1','S1','Other')
silu.subtype.use1$Subtype2=ifelse(silu.subtype.use1$Subtype=='S2','S2','Other')
silu.subtype.use1$Subtype3=ifelse(silu.subtype.use1$Subtype=='S3','S3','Other')

silu.geneList1=getGeneFC(gene.exp=silu.exp[,silu.subtype.use1$Samples],
                         group=silu.subtype.use1$Subtype1,
                         ulab='S1',dlab='Other')
silu.geneList2=getGeneFC(gene.exp=silu.exp[,silu.subtype.use1$Samples],
                         group=silu.subtype.use1$Subtype2,
                         ulab='S2',dlab='Other')
silu.geneList3=getGeneFC(gene.exp=silu.exp[,silu.subtype.use1$Samples],
                         group=silu.subtype.use1$Subtype3,
                         ulab='S3',dlab='Other')

set.seed(123456)
h.all.gmt<-read.gmt("h.all.v7.5.1.entrez.gmt")
silu.hallmark.gsea1<-GSEA(silu.geneList1,TERM2GENE = h.all.gmt,seed=T)
silu.hallmark.gsea2<-GSEA(silu.geneList2,TERM2GENE = h.all.gmt,seed=T)
silu.hallmark.gsea3<-GSEA(silu.geneList3,TERM2GENE = h.all.gmt,seed=T)

library(enrichplot)
library(ggplot2)

silu.hallmark.gsea.res1=silu.hallmark.gsea1@result
silu.hallmark.gsea.res2=silu.hallmark.gsea2@result
silu.hallmark.gsea.res3=silu.hallmark.gsea3@result

rownames(silu.hallmark.gsea.res1)=gsub("HALLMARK_","",rownames(silu.hallmark.gsea.res1))
rownames(silu.hallmark.gsea.res2)=gsub("HALLMARK_","",rownames(silu.hallmark.gsea.res2))
rownames(silu.hallmark.gsea.res3)=gsub("HALLMARK_","",rownames(silu.hallmark.gsea.res3))

silu.hallmark.union=Reduce(union,list(rownames(silu.hallmark.gsea.res1),
                                      rownames(silu.hallmark.gsea.res2),
                                      rownames(silu.hallmark.gsea.res3)))
length(silu.hallmark.union)#33

silu.hallmark.heatmap.dat=matrix(0,nrow = 3,ncol = length(silu.hallmark.union))
rownames(silu.hallmark.heatmap.dat)=c('S1 vs OTHER', 'S2 vs OTHER','S3 vs OTHER')
colnames(silu.hallmark.heatmap.dat)=silu.hallmark.union

silu.hallmark.heatmap.dat[1,match(rownames(silu.hallmark.gsea.res1),colnames(silu.hallmark.heatmap.dat))]=silu.hallmark.gsea.res1$NES
silu.hallmark.heatmap.dat[2,match(rownames(silu.hallmark.gsea.res2),colnames(silu.hallmark.heatmap.dat))]=silu.hallmark.gsea.res2$NES
silu.hallmark.heatmap.dat[3,match(rownames(silu.hallmark.gsea.res3),colnames(silu.hallmark.heatmap.dat))]=silu.hallmark.gsea.res3$NES
range(silu.hallmark.heatmap.dat)


pdf('results/02.subtype_GSEA/Fig2A.pdf',height = 10,width = 6)
Heatmap(as.matrix(t(silu.hallmark.heatmap.dat)), name = "NES"
        , col = circlize::colorRamp2(c(-3, 0, 3), c('#66A61E', 'white', '#E7298A'))
        , border = TRUE
        , show_column_names = T
        , show_column_dend = F
        , show_row_dend =F
        , cluster_columns=T
        , cluster_rows=T
        , rect_gp = gpar(col = "white", lwd = 1)
        , row_names_gp = gpar(fontsize = 10))
dev.off()

head(silu.subtype.cli)
silu.subtype.cli %>%
  ggplot(aes(x=Subtype,y=ICR.Score,fill=Subtype))+
  geom_violin()+  
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  stat_compare_means(aes(group=Subtype), label = "p.format", method = 'kruskal.test')+
  scale_fill_manual(values =subtype.cols)+theme_bw()+ylab('ICR Score')+
  theme(text = element_text(family = 'Times',size = 15),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave('results/02.subtype_GSEA/Fig2B.pdf',height = 5,width = 6)


# library(progeny)
# silu.pathway.activ=progeny(as.matrix(silu.exp),scale = T)
# dim(silu.pathway.activ)
# range(silu.pathway.activ)
# head(silu.pathway.activ)
load("results/02.subtype_GSEA/silu.pathway.activ.RData")
silu.pathway.activ[1:5,1:5]

pathway.df1=silu.pathway.activ[silu.subtype.cli$Samples,]
pathway.df1=as.data.frame(pathway.df1)

pathway.df1$Subtype=silu.subtype.cli$Subtype
pathway.df1=melt(pathway.df1)
head(pathway.df1)
table(pathway.df1$variable)
pathway.df1 %>%
  ggplot(aes(x=variable,y=value,fill=Subtype))+
  geom_boxplot()+stat_compare_means(aes(group=Subtype), label = "p.signif", method = 'kruskal.test')+
  scale_fill_manual(values =subtype.cols)+
  xlab('')+ylab('Score')+ggtitle('PROGENy')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggsave('results/02.subtype_GSEA/Fig2c.pdf',height = 5,width = 6)


#03.##########
dir.create('results/03.m6A_model')
pre.genes=consen_gene
length(pre.genes)
#9

silu_model_data <- cbind(silu.cli[, c("DSS.time", "DSS")],
                         t(silu.exp[pre.genes, silu.cli$Samples]))
colnames(silu_model_data) <- gsub('-', '_', colnames(silu_model_data))

##LASSO####
library(glmnet)
silu.lasso=get_riskscore.lasso(dat = silu_model_data[,-c(1:2)],
                               os = silu_model_data$DSS,
                               os.time = silu_model_data$DSS.time )
length(silu.lasso$lasso.gene)
silu.lasso$plot


######
fmla <- as.formula(paste0("Surv(DSS.time, DSS) ~"
                          ,paste0(silu.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(silu_model_data))
cox=step(cox)

lan <- coef(cox)
lan
length(lan)
paste0(round(lan, 3), '*', names(lan),collapse = '+')
#"0.559*METTL3+-1.327*YTHDC2+0.125*IGF2BP3"

gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
gene.coef.fig=ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#FDB462", "#6154AB")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,2),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,2), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),legend.position = 'top',
                    panel.grid.major = element_blank(),panel.grid.minor = element_blank())
gene.coef.fig

library(survminer)
gene.forest=ggforest(cox, data = silu_model_data, 
                     main = "Hazardratio", fontsize =1.0, 
                     noDigits = 2)
gene.forest






##silu#########
risktype.col=c("#E6AB02","#66A61E")
risk.silu=as.numeric(lan%*%as.matrix(t(silu_model_data[silu.subtype.cli$Samples,names(lan)])))
silu.risktype.cli=data.frame(silu.subtype.cli,Riskscore=risk.silu)
silu.risktype.cli$Risktype=ifelse(silu.risktype.cli$Riskscore>median(risk.silu),'High','Low')
head(silu.risktype.cli)

silu.roc.DSS=ggplotTimeROC(time = silu.risktype.cli$DSS.time,
                           status = silu.risktype.cli$DSS,
                           score = silu.risktype.cli$Riskscore,mks = c(1:5))
silu.roc.DSS

silu.km.DSS=ggsurvplot(fit=survfit( Surv(DSS.time/12, DSS) ~ Risktype,
                                    data = silu.risktype.cli),
                       data=silu.risktype.cli,
                       conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                       surv.median.line = 'hv',title='AC-ICAM',
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ggtheme = custom_theme(),
                       legend = 'top', # 
                       legend.title = "Risktype",legend.labs=c('High','Low'))
silu.km.DSS=mg_merge_plot(silu.km.DSS$plot,silu.km.DSS$table,nrow=2,heights = c(2.5,1),align = 'v')
silu.km.DSS


##GSE33113########
# gset <- getGEO("GSE33113", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# saveRDS(gset,file = 'origin_datas/GEO/GSE33113.rds')
GSE33113=readRDS('origin_datas/GEO/GSE33113.rds')
GSE33113.cli=pData(GSE33113)
head(GSE33113.cli)
GSE33113.cli=data.frame(Samples=GSE33113.cli$geo_accession,
                        DSS=GSE33113.cli$`meta or recurrence within 3 years:ch1`,
                        DSS.time=GSE33113.cli$`time to meta or recurrence:ch1`)
GSE33113.cli$DSS
GSE33113.cli=GSE33113.cli%>%drop_na(DSS)
GSE33113.cli$DSS=ifelse(GSE33113.cli$DSS=='yes',1,0)
rownames(GSE33113.cli)=GSE33113.cli$Samples

load("origin_datas/GEO/GSE33113.exp.RData")
range(GSE33113.exp,na.rm = T)
GSE33113.exp=log2(GSE33113.exp+1)
GSE33113.exp=GSE33113.exp[,GSE33113.cli$Samples]

###
GSE33113_model_data <- cbind(GSE33113.cli[, c("DSS.time", "DSS")],
                             t(GSE33113.exp[pre.genes, GSE33113.cli$Samples]))
colnames(GSE33113_model_data) <- gsub('-', '_', colnames(GSE33113_model_data))


risk.GSE33113=as.numeric(lan%*%as.matrix(t(GSE33113_model_data[GSE33113.cli$Samples,names(lan)])))
GSE33113.risktype.cli=data.frame(GSE33113.cli,Riskscore=risk.GSE33113)
GSE33113.risktype.cli$Risktype=ifelse(GSE33113.risktype.cli$Riskscore>median(risk.GSE33113),'High','Low')

GSE33113.risktype.cli=crbind2DataFrame(GSE33113.risktype.cli)
GSE33113.roc=ggplotTimeROC(time = GSE33113.risktype.cli$DSS.time,
                           status = GSE33113.risktype.cli$DSS,
                           score = GSE33113.risktype.cli$Riskscore,mks = c(1:5))
GSE33113.roc

GSE33113.km.DSS=ggsurvplot(fit=survfit( Surv(DSS.time/365, DSS) ~ Risktype,
                                        data = GSE33113.risktype.cli),
                           data=GSE33113.risktype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                           surv.median.line = 'hv',title='GSE33113',
                           linetype = c("solid", "dashed","strata")[1],
                           palette = risktype.col,ggtheme = custom_theme(),
                           legend ='top', # 
                           legend.title = "Risktype",legend.labs=c('High','Low'))
GSE33113.km.DSS=mg_merge_plot(GSE33113.km.DSS$plot,GSE33113.km.DSS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE33113.km.DSS


pdf('results/03.m6A_model/Fig3_1.pdf',height = 4,width = 16)
mg_merge_plot(silu.lasso$plot,gene.forest,ncol=2,widths = c(1.2,1),labels = c('','C'))
dev.off()

pdf('results/03.m6A_model/Fig3_2.pdf',height = 8.5,width = 10)
mg_merge_plot(silu.roc.DSS,silu.km.DSS,GSE33113.roc,GSE33113.km.DSS,ncol=2,nrow=2,labels = LETTERS[4:7])
dev.off()

library(tinyarray)
pdf('results/03.m6A_model/Fig3_3.pdf',height = 4.25,width = 6,onefile = F)
draw_heatmap(silu.exp[c('YTHDC2','METTL3','IGF2BP3'),order(silu.risktype.cli$Riskscore)],
             group_list = factor(silu.risktype.cli$Risktype[order(silu.risktype.cli$Riskscore)]),
             cluster_cols = F,cluster_rows=F,
             show_rownames = T,
             legend = T,
             annotation_legend = T,
             #n_cutoff = 3,
             color_an = risktype.col)
dev.off()

pdf('results/03.m6A_model/Fig3_4.pdf',height = 4.25,width = 6,onefile = F)
draw_heatmap(GSE33113.exp[c('YTHDC2','METTL3','IGF2BP3'),order(GSE33113.risktype.cli$Riskscore)],
             group_list = factor(GSE33113.risktype.cli$Risktype[order(GSE33113.risktype.cli$Riskscore)]),
             cluster_cols = F,cluster_rows=F,
             border_color = NA,
             show_rownames = T,legend = T,
             annotation_legend = T,
             #n_cutoff = 3,
             color_an = risktype.col)
dev.off()


#04.############
dir.create('results/04.nomogram')
#########
silu.risktype.cli.use=silu.risktype.cli[order(silu.risktype.cli$Riskscore),]
colnames(silu.risktype.cli.use)
cli.colors=list()
for(i in c(18,4,11:14)){
  var=silu.risktype.cli.use[i]
  var.clas=names(table(silu.risktype.cli.use[,i]))
  var.color=brewer.pal(6,"Set2")[1:length(var.clas)]
  names(var.color)=var.clas
  cli.colors=c(cli.colors,list(var.color))
}
names(cli.colors)=colnames(silu.risktype.cli.use)[c(18,4,11:14)]

library(ComplexHeatmap)

median(risk.silu)
silu.riskscore.barplot=columnAnnotation(Riskscore = anno_barplot(silu.risktype.cli.use$Riskscore,
                                                                 baseline =median(risk.silu),
                                                                 bar_width=1,
                                                                 gp=gpar(fill=ifelse(silu.risktype.cli.use$Riskscore>median(risk.silu),
                                                                                     risktype.col[1],risktype.col[2])))
                                        ,annotation_name_side ='left'
                                        ,height =unit(4,'inches'))
draw(silu.riskscore.barplot)

silu.cli.heatmap=columnAnnotation(df = silu.risktype.cli.use[,c(18,4,11:14)]
                                  ,annotation_name_side='left'
                                  ,annotation_height =unit(4,'inches')
                                  ,col = cli.colors)

ht_list=silu.riskscore.barplot %v% silu.cli.heatmap
pdf('results/04.nomogram/Fig4a.pdf',height = 8,width = 12,onefile = F)
ht_list
dev.off()

chisq.test(table(silu.risktype.cli$Risktype,silu.risktype.cli$Age1))
#p-value = 0.01843
chisq.test(table(silu.risktype.cli$Risktype,silu.risktype.cli$Gender))
#p-value = 0.7368
chisq.test(table(silu.risktype.cli$Risktype,silu.risktype.cli$Path.Tumor.Stage))
#p-value = 0.00197
chisq.test(table(silu.risktype.cli$Risktype,silu.risktype.cli$Path.Nodes.Stage))
#p-value = 0.0003965
chisq.test(table(silu.risktype.cli$Risktype,silu.risktype.cli$Path.Metastasis.Stage))
#p-value = 0.06853
chisq.test(table(silu.risktype.cli$Risktype,silu.risktype.cli$AJCC.Path.Stage))
#p-value = 0.001426

###########
head(silu.risktype.cli)
silu_cox_datas=silu.risktype.cli
colnames(silu_cox_datas)
silu_cox_datas=crbind2DataFrame(silu_cox_datas)

table(silu_cox_datas$AJCC.Path.Stage)
silu_cox_datas$AJCC.Path.Stage[silu_cox_datas$AJCC.Path.Stage=='Stage1'|silu_cox_datas$AJCC.Path.Stage=='Stage2']<-'Stage1+2'
silu_cox_datas$AJCC.Path.Stage[silu_cox_datas$AJCC.Path.Stage=='Stage3'|silu_cox_datas$AJCC.Path.Stage=='Stage4']<-'Stage3+4'

table(silu_cox_datas$Path.Tumor.Stage)
silu_cox_datas$Path.Tumor.Stage[silu_cox_datas$Path.Tumor.Stage=='T1'|silu_cox_datas$Path.Tumor.Stage=='T2']<-'T1+T2'
silu_cox_datas$Path.Tumor.Stage[silu_cox_datas$Path.Tumor.Stage=='T3'|silu_cox_datas$Path.Tumor.Stage=='T4']<-'T3+T4'

table(silu_cox_datas$Path.Nodes.Stage)
silu_cox_datas$Path.Nodes.Stage[silu_cox_datas$Path.Nodes.Stage=='N1'|silu_cox_datas$Path.Nodes.Stage=='N2']<-'N1+N2'

table(silu_cox_datas$Path.Metastasis.Stage)




#Age
Age_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~Age,
                             data=silu_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

##Sex
Gender_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~Gender,
                                data=silu_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

colnames(silu_cox_datas)
# #stage
# Stage_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~AJCC.Path.Stage,
#                                data=silu_cox_datas))
# Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
#                                 HR = round(Stage_sig_cox[[7]][,2],3),
#                                 lower.95 = round(Stage_sig_cox[[8]][,3],3),
#                                 upper.95 = round(Stage_sig_cox[[8]][,4],3),
#                                 p.value=round(Stage_sig_cox[[7]][,5],3))
# Stage_sig_cox_dat

#t_stage
t_stage_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~Path.Tumor.Stage,
                                 data=silu_cox_datas))
t_stage_sig_cox_dat <- data.frame(Names=rownames(t_stage_sig_cox[[8]]),
                                  HR = round(t_stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(t_stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(t_stage_sig_cox[[8]][,4],3),
                                  p.value=round(t_stage_sig_cox[[7]][,5],3))
t_stage_sig_cox_dat

#n_stage
n_stage_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~Path.Nodes.Stage,
                                 data=silu_cox_datas))
n_stage_sig_cox_dat <- data.frame(Names=rownames(n_stage_sig_cox[[8]]),
                                  HR = round(n_stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(n_stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(n_stage_sig_cox[[8]][,4],3),
                                  p.value=round(n_stage_sig_cox[[7]][,5],3))
n_stage_sig_cox_dat

#m_stage
m_stage_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~Path.Metastasis.Stage,
                                 data=silu_cox_datas))
m_stage_sig_cox_dat <- data.frame(Names=rownames(m_stage_sig_cox[[8]]),
                                  HR = round(m_stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(m_stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(m_stage_sig_cox[[8]][,4],3),
                                  p.value=round(m_stage_sig_cox[[7]][,5],3))
m_stage_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(DSS.time, DSS)~Riskscore,
                                 data=silu_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     # Stage_sig_cox_dat,
                     t_stage_sig_cox_dat,
                     n_stage_sig_cox_dat,
                     m_stage_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <-c('Age','Gender','Tumor.Stage','Nodes.Stage','Metastasis.Stage','Riskscore')
data.sig$Features=rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/04.nomogram/Univariate.pdf',height = 4,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='#4169E1',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()

#########
head(silu.risktype.cli)
silu_cox_datas=silu.risktype.cli
colnames(silu_cox_datas)
silu_cox_datas=crbind2DataFrame(silu_cox_datas)
#write.csv(silu_cox_datas,'results/04.nomogram/silu_cox_datas.csv')

muti_sig_cox <- summary(coxph(formula=Surv(DSS.time, DSS)~Path.Tumor.Stage+Path.Nodes.Stage+
                                Path.Metastasis.Stage+Riskscore, 
                              data=silu_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
# write.csv(data.muti,'results/04.nomogram/Multivariate_table.csv')
# rownames(data.muti) <- c('Age','Gender','AJCC.Stage','Tumor.Stage','Nodes.Stage','Metastasis.Stage','Riskscore')
# data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/04.nomogram/Multivariate.pdf',height = 4,width =6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='#4169E1',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()



#########

pdf('results/04.nomogram/nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=silu_cox_datas$Riskscore,
                                Path.Nodes.Stage=silu_cox_datas$Path.Nodes.Stage,
                                Path.Metastasis.Stage=silu_cox_datas$Path.Metastasis.Stage),
                     os = silu_cox_datas$DSS.time,
                     status = silu_cox_datas$DSS,
                     mks = c(1,3,5)
)
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1,3,5))

#05.#########
dir.create('results/05.Risktype.TME')

silu.immu.ssgaea=immu_ssgsea(silu.exp)
silu.immu.ssgaea[1:5,1:5]
tme.df1=silu.immu.ssgaea[silu.risktype.cli$Samples,]#silu.cibersort[silu.risktype.cli$Samples,1:22]
tme.df1=as.data.frame(tme.df1)
tme.df1$Risktype=silu.risktype.cli$Risktype
tme.df1=melt(tme.df1)
head(tme.df1)
table(tme.df1$variable)
tme.df1 %>%
  ggplot(aes(x=variable,y=value,fill=Risktype))+
  geom_boxplot()+stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('')+ylab('Score')+ggtitle('ssGSEA')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggsave('results/05.Risktype.TME/Fig5a.pdf',height = 5,width = 12)



load('results/05.Risktype.TME/silu_durg_ic50_res.RData')
library(ggcorrplot)
library(psych)
IC50_RS_cor <- corr.test(x =silu.risktype.cli$Riskscore,
                         y = silu_durg_ic50_res[silu.risktype.cli$Samples,],
                         method = "spearman",adjust = "BH",ci = F)


IC50_RS_cor_res=data.frame(drugs=colnames(silu_durg_ic50_res))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.3)
IC50_RS_cor_res=IC50_RS_cor_res[IC50_RS_cor_res$p.adj<0.05 &abs(IC50_RS_cor_res$cor)>0.3,]
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)
library(rcartocolor)
ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_continuous(type = "gradient")+
  scale_color_gradient(low = "#6B58EEFF", high = "#E6AB02")+
  geom_segment(aes(yend=drugs,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Drugs')+theme_bw()+
  theme(text = element_text(family = 'Times',size=14))
ggsave('results/05.Risktype.TME/Fig5b.pdf',height = 5,width = 6)



#######GSEA
silu.geneList=getGeneFC(gene.exp=silu.exp[,silu.risktype.cli$Samples],group=silu.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')

h.all.gmt<-read.gmt("h.all.v2023.1.Hs.entrez.gmt")
set.seed(147)
silu.hallmark.gsea<-GSEA(silu.geneList,TERM2GENE = h.all.gmt,seed=T)
library(enrichplot)
library(ggplot2)
library(GseaVis)
gsea.dotplot=dotplotGsea(data = silu.hallmark.gsea,topn = 10,order.by = 'NES')
pdf('results/05.Risktype.TME/risktype_GSEA.pdf',height = 5,width = 8)
gsea.dotplot$plot+theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))+
  ggtitle('High Risk vs Low Risk')
dev.off()
write.csv(silu.hallmark.gsea@result,file = 'results/05.Risktype.TME/risktype_GSEA_res.csv')


save(silu.exp,file = 'results/silu.exp.RData')
save.image(file = 'project.RData')
