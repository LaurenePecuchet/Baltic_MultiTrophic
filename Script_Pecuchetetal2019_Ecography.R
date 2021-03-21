##################################################
###     Multiple Factor Analysis 16.10.19      ###
### Script for Baltic Sea Multi-trophic study  ###
### Pecuchet et al. 2019 Ecography             ###

rm(list=ls())

library(reshape2)
library(FactoMineR)
library(missMDA)
library(ggplot2)
library(gridExtra)
library(factoextra)

options(na.action = "na.omit")
colMax <- function(data) sapply(data, max, na.rm = TRUE) # get maximum value

### Datasets ####

# community weighted mean (CWM) traits
CWMt<-read.table("Multitrophic_CWM_traits_BalticSea.csv") ### Time-series of CWM traits across areas and organisms

summary(CWMt)
unique(CWMt$Trait)


### Multiple Factor Analysis across 3 areas ####
##  Corresponding to Figure 3-5 from Pecuchet et al. 2019, Ecograhy

### Kattegat ####

CWMt_Kat<-subset(CWMt, Area=="Kattegat"&Variable=="abundance")

# Keep yrs with data for at least 3 out of the 4 groups
table(unique(CWMt_Kat[,c("Taxa","Year")]))
CWMt_Kat<-subset(CWMt_Kat, Year%in%subset(as.data.frame(table(unique(CWMt_Kat[,c("Taxa","Year")])$Year)), Freq>2)$Var1)
mat<-dcast(CWMt_Kat, Year~ Taxa+ Trait, value.var="Value")

mat<-mat[,colMax(mat)>c(0.25)] ### Keep traits that are expressed at least during one year on more than 25% of the individuals
mat<-mat[,colMeans(mat,na.rm=T)>c(0.10)] ### Keep traits that are expressed on average during the time-series in at least 10% of the individuals

colnames(mat)<-sub('.*\\_', '', colnames(mat))

Kat_CWM=as.matrix(scale(mat[,-1])) ### Standardize data to compare across surveys
rownames(Kat_CWM)<-mat$Year
vyears.kat<-mat$Year
colnames(Kat_CWM)
b<-8
f<-6
p<-6
z<-8

# Estimate missing CWM traits data using regularised iterative MFA algorithm
res.impute_Kat <-imputeMFA(as.data.frame(Kat_CWM),group = c(b,f,p,z))

# MFA, groups decomposed by organism types
res.mfa.kat <-MFA(res.impute_Kat$completeObs,group = c(b,f,p,z),name.group = c("Benthos","Fish","Phytoplankton","Zooplankton"))

####### Get the Years BIPLOT ##
ind <- get_mfa_ind(res.mfa.kat)
ind

MFA_Trends<-rbind(data.frame(Years=mat$Year, Group="All",ind$coord),
                  data.frame(Years=mat$Year, Group="Benthos",res.mfa.kat$separate.analyses$Benthos$ind$coord),
                  data.frame(Years=mat$Year, Group="Fish",res.mfa.kat$separate.analyses$Fish$ind$coord),
                  data.frame(Years=mat$Year, Group="Phytoplankton",res.mfa.kat$separate.analyses$Fish$ind$coord),
                  data.frame(Years=mat$Year, Group="Zooplankton",res.mfa.kat$separate.analyses$Zooplankton$ind$coord))

## CWM traits loadings (x=Dimension 1 ; y=Dimension 2)
imp.kat<-data.frame(lab=rownames(res.mfa.kat$quanti.var$coord),x=res.mfa.kat$quanti.var$coord[,1], y=res.mfa.kat$quanti.var$coord[,2],
                    cor1=res.mfa.kat$quanti.var$cor[,1], cor2=res.mfa.kat$quanti.var$cor[,2], contrib1=res.mfa.kat$quanti.var$contrib[,1],
                    contrib2=res.mfa.kat$quanti.var$contrib[,2],contrib3=res.mfa.kat$quanti.var$contrib[,3], 
                    col1=c(rep("#A37C27",b),rep("#563838",f),rep("#6A8A82",p),rep("#A7414A",z)))

fviz_screeplot(res.mfa.kat)

colnames(res.mfa.kat$ind$coord)<-c("Dim1","Dim2","Dim3","Dim4","Dim5")

##### Plot the MFA results ##
pdim1.Kat<-ggplot()+
  geom_line(aes(vyears.kat,res.mfa.kat$ind$coord[,1]), size=2)+
  theme_bw()+
  labs(x = "Year", y=paste("Dim 1 (",round(res.mfa.kat$eig[1,2],1),"%)", sep = ""))+
  theme(plot.margin=unit(c(0,0.2,0.05,0.2),"cm"))

pdim2.Kat<-ggplot()+
  geom_line(aes(vyears.kat,res.mfa.kat$ind$coord[,2]), size=2)+
  theme_bw()+
  labs(x = "Year", y=paste("Dim 2 (",round(res.mfa.kat$eig[2,2],1),"%)", sep = ""))+
  theme(plot.margin=unit(c(0,0.2,0.05,0.2),"cm"))


contrib1<-facto_summarize(res.mfa.kat, element = "group", result = "contrib", axes=1)
contrib1$name <- factor(contrib1$name , levels = c("Phytoplankton","Zooplankton","Benthos","Fish"))
pcontrib1<-ggplot(contrib1, aes(name , contrib))+
  geom_bar(stat="identity", position = "identity",fill=c("#A37C27","#563838","#6A8A82","#A7414A"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim1")+
  theme(axis.text.x=element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

contrib2<-facto_summarize(res.mfa.kat, element = "group", result = "contrib", axes=2)
contrib2$name <- factor(contrib1$name , levels = c("Phytoplankton","Zooplankton","Benthos","Fish"))
pcontrib2<-ggplot(contrib2, aes(name , contrib))+
  geom_bar(stat="identity", position = "identity",fill=c("#A37C27","#563838","#6A8A82","#A7414A"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim2")+
  theme(axis.text.x=element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

### CWM traits contribution
imp.kat2<-imp.kat[order(abs(imp.kat$x), decreasing = T),][1:14,]
imp.kat2$lab <- factor(imp.kat2$lab, levels = imp.kat2$lab[order(imp.kat2$x, decreasing = F)])

p_loading1tr<-ggplot(imp.kat2, aes(lab, x))+
  geom_bar(stat="identity", position = "identity",fill=imp.kat2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 1")

imp.kat2<-imp.kat[order(abs(imp.kat$y), decreasing = T),][1:10,]
imp.kat2$lab <- factor(imp.kat2$lab, levels = imp.kat2$lab[order(imp.kat2$y, decreasing = F)])

p_loading2tr<-ggplot(imp.kat2, aes(lab, y))+
  geom_bar(stat="identity", position = "identity",fill=imp.kat2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 2")

### Assemble plots ##
pdim1.Kat2<-pdim1.Kat+annotation_custom(ggplotGrob(pcontrib1), xmin = 1982, xmax = 1994, 
                                        ymin = -4, ymax = -0.5)+ggtitle("(a) Multi-trophic CWM trait dynamic - Dim 1")#+theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

pdim2.Kat2<-pdim2.Kat+annotation_custom(ggplotGrob(pcontrib2), xmin = 1994, xmax = 2007, 
                                        ymin = -2.6, ymax = -0.4)+ggtitle("(b) Multi-trophic CWM trait dynamic - Dim 2")#+theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

p_loading1tr<-p_loading1tr+ggtitle("(c) CWM traits")+theme(plot.title = element_text(size = 10, face = "bold"))
p_loading2tr<-p_loading2tr+ggtitle("(e) CWM traits")+theme(plot.title = element_text(size = 10, face = "bold"))


## Kattegat MFA on CWM traits (without environment)
tiff("Kattegat.tiff", height = 15, width = 20, units = 'cm', compression = "lzw", res = 600)
grid.arrange(pdim1.Kat2, pdim2.Kat2,p_loading1tr,p_loading2tr, ncol=2, heights=c(2,3))
dev.off()



### Central Baltic Sea ####
CWMt_CBS<-subset(CWMt, Area=="EBS"&Variable=="abundance")
# Years with at least 3 out of 4 group
CWMt_CBS<-subset(CWMt_CBS, Year%in%subset(as.data.frame(table(unique(CWMt_CBS[,c("Taxa","Year")])$Year)), Freq>2)$Var1)
mat<-dcast(CWMt_CBS, Year~ Taxa+ Trait, value.var="Value")

mat<-mat[,colMax(mat)>c(0.25)] ### Keep traits that are on average at least 10% (espeially for benthos and zooplankton)
mat<-mat[,colMeans(mat,na.rm=T)>c(0.10)]

colnames(mat)<-sub('.*\\_', '', colnames(mat))
CBS_CWM=as.matrix(scale(mat[,-1])) ### Standardize data to compare across surveys
rownames(CBS_CWM)<-mat$Year
vyears.CBS<-mat$Year
colnames(CBS_CWM)
b<-9
p<-5
f<-6
z<-6

res.impute_CBS <-imputeMFA(as.data.frame(CBS_CWM),group = c(b,p,f,z))

res.mfa.CBS <-MFA(res.impute_CBS$completeObs,group = c(b,p,f,z),name.group = c("Benthos","Phytoplankton","Fish","Zooplankton"))

####### Get the Years BIPLOT ##
ind <- get_mfa_ind(res.mfa.CBS)
ind

MFA_Trends<-rbind(data.frame(Years=mat$Year, Group="All",ind$coord),
                  data.frame(Years=mat$Year, Group="Benthos",res.mfa.CBS$separate.analyses$Benthos$ind$coord),
                  data.frame(Years=mat$Year, Group="Fish",res.mfa.CBS$separate.analyses$Fish$ind$coord),
                  data.frame(Years=mat$Year, Group="Phytoplankton",res.mfa.CBS$separate.analyses$Fish$ind$coord),
                  data.frame(Years=mat$Year, Group="Zooplankton",res.mfa.CBS$separate.analyses$Zooplankton$ind$coord))

## Without adding species as sup
imp.CBS<-data.frame(lab=rownames(res.mfa.CBS$quanti.var$coord),x=res.mfa.CBS$quanti.var$coord[,1], y=res.mfa.CBS$quanti.var$coord[,2],
                    cor1=res.mfa.CBS$quanti.var$cor[,1], cor2=res.mfa.CBS$quanti.var$cor[,2], contrib1=res.mfa.CBS$quanti.var$contrib[,1],
                    contrib2=res.mfa.CBS$quanti.var$contrib[,2],contrib3=res.mfa.CBS$quanti.var$contrib[,3], 
                    col1=c(rep("#A37C27",b),rep("#6A8A82",p),rep("#563838",f),rep("#A7414A",z))
)

## Fit ##
fviz_screeplot(res.mfa.CBS) # % explained variance by each dimensions

colnames(res.mfa.CBS$ind$coord)<-c("Dim1","Dim2","Dim3","Dim4","Dim5")

##### Plot the results ##
pdim1.CBS<-ggplot()+
  geom_line(aes(vyears.CBS,res.mfa.CBS$ind$coord[,1]), size=2)+
  theme_bw()+
  labs(x = "Year", y=paste("Dim 1 (",round(res.mfa.CBS$eig[1,2],1),"%)", sep = ""))+
  theme(plot.margin=unit(c(0,0.2,0.05,0.2),"cm"))

pdim2.CBS<-ggplot()+
  geom_line(aes(vyears.CBS,res.mfa.CBS$ind$coord[,2]), size=2)+
  theme_bw()+
  labs(x = "Year", y=paste("Dim 2 (",round(res.mfa.CBS$eig[2,2],1),"%)", sep = ""))+
  theme(plot.margin=unit(c(0,0.2,0.05,0.2),"cm"))

contrib1<-facto_summarize(res.mfa.CBS, element = "group", result = "contrib", axes=1)
contrib1$name <- factor(contrib1$name , levels = c("Phytoplankton","Zooplankton","Benthos","Fish"))
pcontrib1<-ggplot(contrib1, aes(name , contrib))+
  geom_bar(stat="identity", position = "identity",fill=c("#A37C27","#6A8A82","#563838","#A7414A"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim1")+
  theme(axis.text.x=element_blank(),
        #panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

contrib2<-facto_summarize(res.mfa.CBS, element = "group", result = "contrib", axes=2)
contrib2$name <- factor(contrib1$name , levels = c("Phytoplankton","Zooplankton","Benthos","Fish"))
pcontrib2<-ggplot(contrib2, aes(name , contrib))+
  geom_bar(stat="identity", position = "identity",fill=c("#A37C27","#6A8A82","#563838","#A7414A"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim2")+
  theme(axis.text.x=element_blank(),
        #panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

## CWM traits contribution
imp.CBS2<-imp.CBS[order(abs(imp.CBS$x), decreasing = T),][1:14,]
imp.CBS2$lab <- factor(imp.CBS2$lab, levels = imp.CBS2$lab[order(imp.CBS2$x, decreasing = F)])

p_loading1tr<-ggplot(imp.CBS2, aes(lab, x))+
  geom_bar(stat="identity", position = "identity",fill=imp.CBS2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 1")

imp.CBS2<-imp.CBS[order(abs(imp.CBS$y), decreasing = T),][1:13,]
imp.CBS2$lab <- factor(imp.CBS2$lab, levels = imp.CBS2$lab[order(imp.CBS2$y, decreasing = F)])

p_loading2tr<-ggplot(imp.CBS2, aes(lab, y))+
  geom_bar(stat="identity", position = "identity",fill=imp.CBS2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 2")


## Assemble plots
pdim1.CBS2<-pdim1.CBS+annotation_custom(ggplotGrob(pcontrib1), xmin = 1987, xmax = 1998, 
                                        ymin = 0.8, ymax = 4)+ggtitle("(a) Multi-trophic CWM trait dynamic - Dim 1")#+theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

pdim2.CBS2<-pdim2.CBS+annotation_custom(ggplotGrob(pcontrib2), xmin = 1984, xmax = 1995, 
                                        ymin = -0.05, ymax = 2.35)+ggtitle("(b) Multi-trophic CWM trait dynamic - Dim 2")#+theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

p_loading1tr<-p_loading1tr+ggtitle("(c) CWM traits")+theme(plot.title = element_text(size = 10, face = "bold"))
p_loading2tr<-p_loading2tr+ggtitle("(e) CWM traits")+theme(plot.title = element_text(size = 10, face = "bold"))


## Central Baltic Sea - MFA on CWM traits (without environment)
tiff("CBS.tiff", height = 15, width = 20, units = 'cm', compression = "lzw", res = 600)
grid.arrange(pdim1.CBS2, pdim2.CBS2,p_loading1tr,p_loading2tr, ncol=2, heights=c(2,3))
dev.off()


### Gulf of Riga ####
GoR<-subset(CWMt, Area=="GoR"&Variable=="abundance")
table(GoR$Taxa, GoR$Year) ## Number of taxa per year (when = 0, there is no data for this year*taxa)
table(unique(GoR[,c("Taxa","Year")]))
GoR<-subset(GoR, Year%in%subset(as.data.frame(table(unique(GoR[,c("Taxa","Year")])$Year)), Freq>2)$Var1)

mat<-dcast(GoR, Year~ Taxa+ Trait, value.var="Value")
mat<-mat[,colMax(mat)>c(0.25)] ### Keep traits that are expressed at least 1 year by more than 25% of the individuals in the community
mat<-mat[,colMeans(mat,na.rm=T)>c(0.10)] ### Keep traits that are on average at least 10% (espeially for benthos and zooplankton)

colnames(mat)<-sub('.*\\_', '', colnames(mat))

GoR_CWM=as.matrix(scale(mat[,-1])) ### Standardize data to compare across surveys
rownames(GoR_CWM)<-mat$Year
vYears<-mat$Year
colnames(GoR_CWM)

b<-7
f<-7
p<-6
z<-7

# Estimate missing CWM traits data using regularised iterative MFA algorithm
res.impute_GoR <-imputeMFA(as.data.frame(GoR_CWM),group = c(b,f,p,z))
# MFA, groups decomposed by organism types
res.mfa <-MFA(res.impute_GoR$completeObs,group = c(b,f,p,z),name.group = c("Benthos","Fish","Phytoplankton","Zooplankton"))

####### Get the Years BIPLOT ##
ind <- get_mfa_ind(res.mfa)
ind

MFA_Trends<-rbind(data.frame(Years=mat$Year, Group="All",ind$coord),
                  data.frame(Years=mat$Year, Group="Benthos",res.mfa$separate.analyses$Benthos$ind$coord),
                  data.frame(Years=mat$Year, Group="Fish",res.mfa$separate.analyses$Fish$ind$coord),
                  data.frame(Years=mat$Year, Group="Phytoplankton",res.mfa$separate.analyses$Fish$ind$coord),
                  data.frame(Years=mat$Year, Group="Zooplankton",res.mfa$separate.analyses$Zooplankton$ind$coord))

## Without adding species as sup
imp<-data.frame(lab=rownames(res.mfa$quanti.var$coord),x=res.mfa$quanti.var$coord[,1], y=res.mfa$quanti.var$coord[,2],
                cor1=res.mfa$quanti.var$cor[,1], cor2=res.mfa$quanti.var$cor[,2], contrib1=res.mfa$quanti.var$contrib[,1],
                contrib2=res.mfa$quanti.var$contrib[,2],contrib3=res.mfa$quanti.var$contrib[,3],
                col1=c(rep("#A37C27",b),rep("#563838",f),rep("#6A8A82",p),rep("#A7414A",z)))

fviz_screeplot(res.mfa) # % explained by each dimensions

colnames(res.mfa$ind$coord)<-c("Dim1","Dim2","Dim3","Dim4","Dim5")

##### Plot the results ##
pdim1.GoR<-ggplot()+
  geom_line(aes(vYears,res.mfa$ind$coord[,1]), size=2)+
  theme_bw()+
  labs(x = "Year", y=paste("Dim 1 (",round(res.mfa$eig[1,2],1),"%)", sep = ""))+
  theme(plot.margin=unit(c(0,0.2,0.05,0.2),"cm"))

pdim2.GoR<-ggplot()+
  geom_line(aes(vYears,res.mfa$ind$coord[,2]), size=2)+
  theme_bw()+
  labs(x = "Year", y=paste("Dim 2 (",round(res.mfa$eig[2,2],1),"%)", sep = ""))+
  theme(plot.margin=unit(c(0,0.2,0.05,0.2),"cm"))

contrib1<-facto_summarize(res.mfa, element = "group", result = "contrib", axes=1)
contrib1$name <- factor(contrib1$name , levels = c("Phytoplankton","Zooplankton","Benthos","Fish"))
pcontrib1<-ggplot(contrib1, aes(name , contrib))+
  geom_bar(stat="identity", position = "identity",fill=c("#A37C27","#563838","#6A8A82","#A7414A"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim1")+
  theme(axis.text.x=element_blank(),
        #panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

contrib2<-facto_summarize(res.mfa, element = "group", result = "contrib", axes=2)
contrib2$name <- factor(contrib1$name , levels = c("Phytoplankton","Zooplankton","Benthos","Fish"))
pcontrib2<-ggplot(contrib2, aes(name , contrib))+
  geom_bar(stat="identity", position = "identity",fill=c("#A37C27","#563838","#6A8A82","#A7414A"))+
  geom_hline(yintercept = 25, linetype = 2,color = "red")+
  theme_bw()+xlab("")+ylab("")+
  ggtitle("Contribution (%) to Dim2")+
  theme(axis.text.x=element_blank(),
        #panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 8, face = "bold"))

### CWM traits contribution
imp.GoR2<-imp[order(abs(imp$x), decreasing = T),][1:18,]
imp.GoR2$lab <- factor(imp.GoR2$lab, levels = imp.GoR2$lab[order(imp.GoR2$x, decreasing = F)])

p_loading1tr<-ggplot(imp.GoR2, aes(lab, x))+
  geom_bar(stat="identity", position = "identity",fill=imp.GoR2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 1")

imp.GoR2<-imp[order(abs(imp$y), decreasing = T),][1:10,]
imp.GoR2$lab <- factor(imp.GoR2$lab, levels = imp.GoR2$lab[order(imp.GoR2$y, decreasing = F)])

p_loading2tr<-ggplot(imp.GoR2, aes(lab, y))+
  geom_bar(stat="identity", position = "identity",fill=imp.GoR2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 2")


# Exemple with Environment data ####
GoR_env<-read.table("GoR_Environment.csv")
GoR_env<-subset(GoR_env, Year%in%rownames(GoR_CWM)) #Keep the same years than in the CWM traits data
GoR_env<-as.matrix(scale(GoR_env[,-c(1)])) # Scale variables

GoR_env<-cbind(GoR_env,res.impute_GoR$completeObs)
colnames(GoR_env)

n.GoR<-12 #number of environmental variables

res.mfa.env.GoR <- MFA(GoR_env,
                       group = c(n.GoR,b,f,p,z),
                       name.group = c("Environment","Benthos","Fish","Phytoplankton","Zooplankton"),
                       num.group.sup = c(1)) # Environment as supplementary group, i.e. not influencing the MFA results

res.mfa.env.GoR$quanti.var.sup$coord # Loadings of the environment

imp.env.GoR<-data.frame(lab=c(rownames(res.mfa.env.GoR$quanti.var$coord),rownames(res.mfa.env.GoR$quanti.var.sup$coord)),
                        x=c(res.mfa.env.GoR$quanti.var$coord[,1],res.mfa.env.GoR$quanti.var.sup$coord[,1]), 
                        y=c(res.mfa.env.GoR$quanti.var$coord[,2],res.mfa.env.GoR$quanti.var.sup$coord[,2]), 
                        cor1=c(res.mfa.env.GoR$quanti.var$cor[,1],res.mfa.env.GoR$quanti.var.sup$cor[,1]),
                        cor2=c(res.mfa.env.GoR$quanti.var$cor[,2],res.mfa.env.GoR$quanti.var.sup$cor[,2]), 
                        contrib1=c(res.mfa.env.GoR$quanti.var$contrib[,1],rep(0,12)),
                        contrib2=c(res.mfa.env.GoR$quanti.var$contrib[,2],rep(0,12)), 
                        col1=c(rep("#A37C27",b),rep("#563838",f),rep("#A7414A",p),rep("#6A8A82",z),rep("black",10),rep("grey45",1),rep("black",1)))

imp.env.GoR2<-imp.env.GoR[(sum(c(b,f,p,z))+1):(sum(c(b,f,p,z))+n.GoR),]
imp.env.GoR2<-imp.env.GoR2[order(abs(imp.env.GoR2$x), decreasing = T),][1:10,]
imp.env.GoR2$lab <- factor(imp.env.GoR2$lab, levels = imp.env.GoR2$lab[order(imp.env.GoR2$x, decreasing = F)])

p_loading1env<-ggplot(imp.env.GoR2, aes(lab, x))+
  geom_bar(stat="identity", position = "identity",fill=imp.env.GoR2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 1")

imp.env.GoR2<-imp.env.GoR[(sum(c(b,f,p,z))+1):(sum(c(b,f,p,z))+n.GoR),]
imp.env.GoR2<-imp.env.GoR2[order(abs(imp.env.GoR2$y), decreasing = T),][1:10,]
imp.env.GoR2$lab <- factor(imp.env.GoR2$lab, levels = imp.env.GoR2$lab[order(imp.env.GoR2$y, decreasing = F)])

p_loading2env<-ggplot(imp.env.GoR2, aes(lab, y))+
  geom_bar(stat="identity", position = "identity",fill=imp.env.GoR2$col1)+
  geom_hline(yintercept = 0.5, linetype = 2,color = "red")+
  geom_hline(yintercept = -0.5, linetype = 2,color = "red")+
  theme_bw()+
  coord_flip()+
  labs(x = "", y="Loadings in Dim 2")


### Assemble plots
pdim1.GoR2<-pdim1.GoR+annotation_custom(ggplotGrob(pcontrib1), xmin = 2000, xmax = 2014, 
                                        ymin = 0.3, ymax = 3)+ggtitle("(a) Multi-trophic CWM trait dynamic - Dim 1")#+theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

pdim2.GoR2<-pdim2.GoR+annotation_custom(ggplotGrob(pcontrib2), xmin = 2000, xmax = 2014, 
                                        ymin = 0.8, ymax = 3.7)+ggtitle("(b) Multi-trophic CWM trait dynamic - Dim 2")#+theme(plot.margin=unit(c(0,0.5,0,0.5),"cm"))

p_loading1tr<-p_loading1tr+ggtitle("(c) CWM traits")+theme(plot.title = element_text(size = 10, face = "bold"))
p_loading1env<-p_loading1env+ggtitle("(d) Environment")+theme(plot.title = element_text(size = 10, face = "bold"))
p_loading2tr<-p_loading2tr+ggtitle("(e) CWM traits")+theme(plot.title = element_text(size = 10, face = "bold"))
p_loading2env<-p_loading2env+ggtitle("(f) Environment")+theme(plot.title = element_text(size = 10, face = "bold"))

mat<-rbind(c(1,1,2,2),c(3,4,5,6))

## GoR MFA on CWM traits with environment as supplementary variables
tiff("GoR.tiff", height = 15, width = 25, units = 'cm', compression = "lzw", res = 400)
grid.arrange(pdim1.GoR2, pdim2.GoR2,p_loading1tr, p_loading1env,p_loading2tr, p_loading2env, ncol=4,layout_matrix=mat, heights=c(9,7))
dev.off()


