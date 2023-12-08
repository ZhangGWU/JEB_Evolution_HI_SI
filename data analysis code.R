#################################################################################
############  JEB Evolution of habitat isolation and sexual isolation analysis by Linyi Zhang 2023#################
###################################################################################
########################### package to load ##########################
library(ggplot2)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(bbmle)
library(plyr)
library(patchwork)
library(grid)
library(cowplot)
library(stringr)
library(car)
library(ggpattern)
library(vegan)
library(brms)
library(rstanarm)
library(rstan)
library(ggplot2)
library(gridExtra)
library(brms)
library(lme4)
####################################################################################
######################### Habitat isolation analysis ###############################
####################################################################################

################# read host preference data ############
host.preference<-read.csv("Bt host preference new.csv")
######### curate host preference data ##########
head(host.preference)
host.preference$host.plant<-"Qv"
host.preference$host.plant[host.preference$HP=="Qg"]<-"Qg"
host.preference$Block<-paste(host.preference$Year,host.preference$Block,sep="_")
host.preference$species<-"B. treatae"
host.preference$species[host.preference$HP=="Qv-T"]<-"B. kinyesi"
host.preference$species[host.preference$HP=="Qg"]<-"B. fossoria"
head(host.preference)
#######################################################################################################
#### define a function to calculate the amount of time spend on each host plant ####
H.pre1<-function(data,n){
  for(i in 1:length(data[,1])){
    data$Qv.Times[i]<-length(which(data[i,(14+n):34]=="Qv"))
    data$Qg.Times[i]<-length(which(data[i,(14+n):34]=="Qg"))
  }
  data
}  

host.pre<-H.pre1(host.preference,1)
#### calculate the preference = #times on native host/(times on native host+times on non-native host) ####
host.pre$preference[host.pre$host.plant=="Qv"]<-
  host.pre$Qv.Times[host.pre$host.plant=="Qv"]/(host.pre$Qv.Times[host.pre$host.plant=="Qv"]+host.pre$Qg.Times[host.pre$host.plant=="Qv"])

host.pre$preference[host.pre$host.plant=="Qg"]<-
  host.pre$Qg.Times[host.pre$host.plant=="Qg"]/(host.pre$Qv.Times[host.pre$host.plant=="Qg"]+host.pre$Qg.Times[host.pre$host.plant=="Qg"])

host.pre$Total.times<-host.pre$Qv.Times+host.pre$Qg.Times
#### calculate the latency = time takes to land on one host plant ####
for(i in 1:length(host.pre[,1])){
  host.pre$latency1[i]<-which(host.pre[i,15:34]!="C")[1]
}

for(i in 1:length(host.pre[,1])){
  host.pre$latency2[i]<-which(host.pre[i,15:34]==host.pre[i,"host.plant"])[1]
}
###############################################################################
########## overview of the distribution of host preference ####################
ggplot(data=host.pre,aes(x=preference))+facet_grid(vars(Sex),vars(HP))+geom_histogram()
ggplot(data=host.pre,aes(x=preference))+facet_grid(vars(Method),vars(HP))+geom_histogram()
################################################################################
#################### host preference, only females, dataset curation #######################
host.pre$HP<-factor(host.pre$HP,levels=c("Qv-T","Qv-F","Qg"))
host.pre$species[host.pre$HP=="Qv-T"]<-"B. kinseyi"
host.pre$species[host.pre$HP=="Qv-F"]<-"B. treatae"
host.pre$species[host.pre$HP=="Qg"]<-"B. fossoria"
host.pre$Year<-as.factor(host.pre$Year)
host.pre$count<-1

host.pre.new<-host.pre[host.pre$Total.times!=0,]
host.pre.new$preference[host.pre.new$preference==0]<-host.pre.new$preference[host.pre.new$preference==0]+0.0000001
host.pre.new$preference[host.pre.new$preference==1]<-host.pre.new$preference[host.pre.new$preference==1]-0.0000001
host.pre.new<-host.pre.new[!is.na(host.pre.new$Source),]
host.pre.new$pre.times[host.pre.new$host.plant=="Qv"]<-host.pre.new$Qv.Times[host.pre.new$host.plant=="Qv"]
host.pre.new$pre.times[host.pre.new$host.plant=="Qg"]<-host.pre.new$Qg.Times[host.pre.new$host.plant=="Qg"]

head(host.pre.new)
ddply(host.pre.new,"Source",summarize,N=sum(count))
host.pre.new<-host.pre.new[host.pre.new$Source!="Lake Griffin",] ### remove lake griffin individuals due to very low sample size ##
hp.sum<-ddply(host.pre.new,c("HP","Source"),summarize,N=sum(count))
###### ### write.csv(hp.sum,"hp.sum.csv") #####
##########################################################################################
# compare the strength of host preference between different species #
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi" &
                                 host.pre.new$Source=="Rice"],mu=0.5,alternative = 'greater')
### mean = 0.6161 ####
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi" &
                                 host.pre.new$Source=="Golden Meadow"],mu=0.5,alternative = 'greater')
### mean = 0.6597 ####
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi" &
                                 host.pre.new$Source=="Picayune "],mu=0.5,alternative = 'greater')
### mean = 0.6273 ####
t.test(host.pre.new$preference[host.pre.new$species=="B. kinseyi"],mu=0.5,alternative = 'greater')
### mean = 0.6289 ####
hp.pop<-c(0.6161,0.6597,0.6273)
hp.pop.name<-c("Rice","GM","PY")
hp.pop<-data.frame(hp.pop.name,hp.pop)
hp.mean<-0.6289
hp.species<-"B. kinseyi"
hp.mean<-data.frame(hp.species,hp.mean)
hp.pop$species<-"B. kinseyi"
##########################################################################################
#### summarize the host preference data across different wasp species and sites #####
hp.sum<-ddply(host.pre.new,c("species","Source"),summarize,N=sum(count),pref=mean(preference),
              pre.SE=sd(preference)/sqrt(N))

hp.sum$species<-factor(hp.sum$species,levels=c("B. kinseyi","B. fossoria","B. treatae"))
hp.sum<-hp.sum[order(hp.sum$species),]
hp.sum$site<-c("GM","PY","RU","ABS","DCK","LL","CC","KRE")
ddply(host.pre.new,c("species"),summarize,N=sum(count),pref=mean(preference),
      pre.SE=sd(preference)/sqrt(N))

ddply(hp.sum,c("species"),summarize,N=sum(N))
##########################################################################################
########## testing whether host preference differ from 0.5 ################
t.test(host.pre$preference[host.pre$HP=="Qv-T"& host.pre$Total.times!=0],mu=0.5,alternative = 'greater')
t.test(host.pre$preference[host.pre$HP=="Qv-F"& host.pre$Total.times!=0],mu=0.5,alternative = 'greater')
t.test(host.pre$preference[host.pre$HP=="Qg" & host.pre$Total.times!=0],mu=0.5,alternative = 'greater')
###########################################################################
#### compare host preference among three species using GLMMM ####
options(digits = 5)
host.pre.new$species<-as.factor(host.pre.new$species)
host.pre.new$species<-factor(host.pre.new$species,levels=c("B. kinseyi","B. fossoria","B. treatae"))
host.pre.new$preference[host.pre.new$preference==1]<-1-0.01
host.pre.new$preference[host.pre.new$preference==0]<-0.01

hp<-glmmTMB(pre.times/Total.times~Method*species+(1|Source),family=betabinomial(link = "logit"),weights=Total.times,
            data=host.pre.new)
summary(hp)

pairs(emmeans(hp,"species"))
ddply(hp.sum,c("species"),summarize,TN=sum(N))

##### generating Figure 2A ###############
hp.sum$N<-ddply(hp.sum,c("species"),summarize,TN=sum(N))$TN
hp.sum$N<-ddply(hp.sum,c("species"),summarize,TN=sum(N))$TN

hp.sum[9,"species"]<-"B. treatae"
hp.sum[9,"Source"]<-"ylank"
hp.sum[9,"N"]<-NA
hp.sum[9,"pref"]<-NA
hp.sum[9,"pref.SE"]<-NA
host.pre.new[427,]<-host.pre.new[426,]
host.pre.new[427,"HP"]<-"Qv-F"
host.pre.new[427,"species"]<-"B. treatae"
host.pre.new[427,"Source"]<-"ylank"
host.pre.new[427,38:46]<-NA
host.pre.new[427,"preference"]<--0.3
host.pre.new$Source<-factor(host.pre.new$Source,levels=hp.sum$Source)
hp.sum$Source<-as.factor(hp.sum$Source)
hp.sum$Source<-factor(hp.sum$Source,levels=hp.sum$Source)
HP.fig<-ggplot()+
  geom_hline(yintercept = 0.5, linetype="dashed",size=0.5)+
  geom_bar(data = hp.sum,aes(y=pref,x=species,fill=Source,color=species),width=0.6,color="black",stat="identity", position=position_dodge(width=0.95))+
  geom_errorbar(data = hp.sum,aes(ymin=pref-pre.SE,ymax=pref+pre.SE,x=species,fill=Source),stat="identity", width=.2,
                position=position_dodge(width=0.95))+
  geom_segment(aes(x=c(0.6,1.6,2.58),y=c(0.75,0.78,0.9),xend=c(1.4,2.4,3.1),yend=c(0.75,0.78,0.9)))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1.05),breaks=seq(0,1,0.25))+
  scale_fill_manual(values=c("mediumspringgreen","mediumspringgreen","mediumspringgreen",
                             "lightblue","lightblue","lightblue","burlywood1","burlywood1","burlywood1"))+
  scale_color_manual(values=c("darkgreen","blue3","orange4"))+
  ylab("Host preference")+xlab("")+
  annotate("text",x=c(1,2,2.9),y=c(0.82,0.86,0.96),label=c("a","a","b"),size=20)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(colour="black"),
        axis.text.y = element_text(size=24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=30,margin=margin(0,20,0,0)),
        plot.margin=unit(c(1,0,3,0),"lines"),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")

HP.fig
HP.fig1<-ggdraw(add_sub(HP.fig,c("GM","PY","RU","ABS","DCK","LL","CC","KRE"), 
                        x=c(0.08,0.19,0.29,0.4,0.5,0.61,0.71,0.82),y=0.8,size=18,color="black"))

HP.fig1

boxes <- data.frame(x=c(0.15,0.44,0.71),y=c(0.07,0.07,0.07))

HP.fig2<-ggdraw(HP.fig1) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +c(0.22,0.22,0.15), ymin = y, ymax = y +0.0005),colour = "black", fill = "black")

HP.fig3<-ggdraw(add_sub(HP.fig2,c(expression(italic("B. kinseyi")),expression(italic("B. fossoria")),expression(italic("B. treatae"))), 
                        x=c(0.27,0.55,0.78),y=1.8,size=22,color="black"))

HP.fig4<-ggdraw(add_sub(HP.fig3,c("Allopatry","Sympatry","Sympatry"), 
                        x=c(0.27,0.55,0.78),y=2,size=22,color="black"))
HP.fig4

#### save as 10x12 pdf #####
################# comapre habitat isolation ###################
hp.sum<-hp.sum[-9,]
hp.sum$pre.Qv<-hp.sum$pref
hp.sum$pre.Qv[hp.sum$species=="B. fossoria"]<-1-hp.sum$pref[hp.sum$species=="B. fossoria"]
hp.sum$species<-as.character(hp.sum$species)
hp.sum$Source<-as.character(hp.sum$Source)

###################################################################
host.pre.new$Qv.pre<-host.pre.new$preference
host.pre.new$Qv.pre[host.pre.new$HP=="Qg"]<-1-host.pre.new$preference[host.pre.new$HP=="Qg"]
########## bootstrap of HI combining all populations ##############
HI.bootsum<-data.frame(matrix(nrow=3,ncol=4))
colnames(HI.bootsum)<-c("pair","HI","low.HI","high.HI")
HI.bootsum$pair<-c("B. kinseyi x B. treatae","B. kinseyi x B. fossoria",
                   "B. treatae x B. fossoria")

HI.bootsum$pair<-factor(HI.bootsum$pair,level=c("B. kinseyi x B. treatae","B. kinseyi x B. fossoria",
                                                "B. treatae x B. fossoria"))

HI.bootsum$HI[1]<-abs(mean(host.pre.new$Qv.pre[host.pre.new$HP=="Qv-T"])-
                        mean(host.pre.new$Qv.pre[host.pre.new$HP=="Qv-F"]))
HI.bootsum$HI[2]<-abs(mean(host.pre.new$Qv.pre[host.pre.new$HP=="Qv-T"])-
                        mean(host.pre.new$Qv.pre[host.pre.new$HP=="Qg"]))
HI.bootsum$HI[3]<-abs(mean(host.pre.new$Qv.pre[host.pre.new$HP=="Qv-F"])-
                        mean(host.pre.new$Qv.pre[host.pre.new$HP=="Qg"]))

HI.TF<-c()
HI.TG<-c()
HI.FG<-c()
unique(host.pre.new$Source[host.pre.new$HP=="Qv-T"])
unique(host.pre.new$Source[host.pre.new$HP=="Qv-F"])
unique(host.pre.new$Source[host.pre.new$HP=="Qg"])
for (i in 1:10000){
  Qvt.new<-sample(host.pre.new$Qv.pre[host.pre.new$HP=="Qv-T"],replace=TRUE)
  Qvf.new<-sample(host.pre.new$Qv.pre[host.pre.new$HP=="Qv-F"],replace=TRUE)
  Qg.new<-sample(host.pre.new$Qv.pre[host.pre.new$HP=="Qg"],replace=TRUE)
  HI.TF[i]<-abs(mean(Qvt.new)-mean(Qvf.new))
  HI.TG[i]<-abs(mean(Qvt.new)-mean(Qg.new))
  HI.FG[i]<-abs(mean(Qvf.new)-mean(Qg.new))
}

HI.bootsum[1,c(3,4)]<-quantile(HI.TF,c(0.025,0.975))
HI.bootsum[2,c(3,4)]<-quantile(HI.TG,c(0.025,0.975))
HI.bootsum[3,c(3,4)]<-quantile(HI.FG,c(0.025,0.975))

##### generating Figure 2B ###############
HI_plot<-ggplot(data=HI.bootsum,aes(x=pair,y=HI,color=pair,fill=pair))+
  geom_errorbar(aes(ymin=low.HI,ymax=high.HI),width=0.08)+
  geom_point(size=6)+
  scale_color_manual(values=c("green4","blue3","burlywood3"))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.03,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=24,colour="black"),
        axis.title=element_text(size=28),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(c(1,0,3,0),"lines"),
        plot.title = element_text(size = 28, face = "bold",hjust = 0.5))+
  labs(x=" ",y="Habitat isolation")

HI_plot1<-ggdraw(add_sub(HI_plot,c(expression(paste(italic("B. kinseyi")," x ", italic("B. treatae"))),
                                   expression(paste(italic("B. kinseyi")," x ", italic("B. fossoria"))),
                                   expression(paste(italic("B. treatae")," x ", italic("B. fossoria")))), 
                         x=c(0.19,0.51,0.83),size=22,y=0.5,color="black"))
HI_plot1

HI_plot2<-ggdraw(add_sub(HI_plot1,c("Host plant: ","Same","Different","Different"), 
                         x=c(0.1,0.25,0.54,0.83),y=1.5,size=22,color="black"))

HI_plot2

HI_plot3<-ggdraw(add_sub(HI_plot2,c("Geography: ","Allopatry","Allopatry","Sympatry"), 
                         x=c(0.11,0.25,0.54,0.83),y=1.2,size=22,color="black"))
HI_plot3

plot_grid(HP.fig4,HI_plot3,nrow=2,labels=c("A","B"),label_size = 38)
### save as 15*18 #####
#######################################################################################
############################ sexual isolation  ########################################
#######################################################################################
#### read mating preference data ####
Mating.T<-read.csv("mate.preference.csv")
Mating.T<-Mating.T[,-1]
unique(Mating.T$Source..M.)
unique(Mating.T$Source..F.)
levels(Mating.T$Source..F.)
sf.data<-list()
Mating.T$pair5<-Mating.T$pair4
Mating.T$pair5<-as.character(Mating.T$pair5)
Mating.T$Source..F.<-as.character(Mating.T$Source..F.)
Mating.T$Source..M.<-as.character(Mating.T$Source..M.)
Mating.T$pair5[Mating.T$Source..F.!=Mating.T$Source..M. & Mating.T$pair4 =="same"]<-"different.pop"
Mating.T$pair6<-"same"
Mating.T$pair6[Mating.T$Source..M.!=Mating.T$Source..F.]<-"different"
Mating.T$pairs<-paste(Mating.T$Source..M.,Mating.T$Source..F.,sep="_")
Mating.Tnew<-Mating.T[,-c(5,6,7,12,17,36,37,38,39,40,42,43,44,45)]

colnames(Mating.Tnew)
colnames(Mating.Tnew)[25]<-"host.pair"
Mating.Tnew$sp.pair<-"same"
Mating.Tnew$sp.pair[Mating.Tnew$HP.M.!=Mating.Tnew$HP.F.]<-"different"
Mating.Tnew$geography<-"sympatry"
Mating.Tnew$geography[Mating.Tnew$HP.M. =="Qv-T"|Mating.Tnew$HP.F.=="Qv-T"]<-"allopatry"
Mating.Tnew$pair<-paste(Mating.Tnew$geography,Mating.Tnew$host.pair,sep="_")

################################ mate preference comparison #########################################
### sympatry,different host, Qv-F and Qg mate preference ####
head(Mating.T)

Mating.Tnew1<-Mating.Tnew[Mating.Tnew$pair %in% c("sympatry_same","sympatry_different"),]

M.m1<-glmer(Mt~sp.pair+(1|Source..M.)+(1|Source..F.),data=Mating.Tnew1,
            family="binomial")
summary(M.m1)

Mating.Tnew1.lsm<-data.frame(lsmeans(M.m1,~sp.pair))
Mating.Tnew1.lsm$lsmean<-exp(Mating.Tnew1.lsm$lsmean)/(exp(Mating.Tnew1.lsm$lsmean)+1)
Mating.Tnew1.lsm$asymp.LCL<-exp(Mating.Tnew1.lsm$asymp.LCL)/(exp(Mating.Tnew1.lsm$asymp.LCL)+1)
Mating.Tnew1.lsm$asymp.UCL<-exp(Mating.Tnew1.lsm$asymp.UCL)/(exp(Mating.Tnew1.lsm$asymp.UCL)+1)

Mating.Tnew1.sum<-ddply(Mating.Tnew1,c("HP.M.","HP.F."),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew1.sum$Mtcount<-Mating.Tnew1.sum$Mt*Mating.Tnew1.sum$count
Mating.Tnew1.sum

Mating.Tnew1.sum1<-ddply(Mating.Tnew1,c("HP.M.","HP.F.","Source..M.","Source..F."),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew1.sum1$pair<-"same"
Mating.Tnew1.sum1$pair[Mating.Tnew1.sum1$HP.M.=="Qv-F" & Mating.Tnew1.sum1$HP.F.=="Qg"]<-"different"
Mating.Tnew1.sum1$pair[Mating.Tnew1.sum1$HP.M.=="Qg" & Mating.Tnew1.sum1$HP.F.=="Qv-F"]<-"different"
Mating.Tnew1.sum1<-Mating.Tnew1.sum1[Mating.Tnew1.sum1$count>2,]

### sympatry,different host: Qv-F and Qg SI with bootstrap ###
head(Mating.Tnew1)

Qg<-c(rep(0,178),rep(1,38))
Qvf<-c(rep(0,38),rep(1,9))
Qvfg<-c(rep(0,64),rep(1,7))
Qvgf<-c(rep(0,88),rep(1,8))
home.fre<-(sum(Qvf)/length(Qvf)+sum(Qg)/length(Qg))/2
heto.fre<-(sum(Qvfg)/length(Qvfg)+sum(Qvgf)/length(Qvgf))/2
FGSI.mean<-1-2*heto.fre/(home.fre+heto.fre)
FG.SI<-c()
for (i in 1:10000){
  Qg.new<-sample(Qg,replace=TRUE)
  Qvf.new<-sample(Qvf,replace=TRUE)
  Qvfg.new<-sample(Qvfg,replace=TRUE)
  Qvgf.new<-sample(Qvgf,replace=TRUE)
  home.fre<-(sum(Qvf.new)/length(Qvf.new)+sum(Qg.new)/length(Qg.new))/2
  heto.fre<-(sum(Qvfg.new)/length(Qvfg.new)+sum(Qvgf.new)/length(Qvgf.new))/2
  FG.SI[i]<-1-2*heto.fre/(home.fre+heto.fre)
}
quantile(FG.SI, c(0.025,0.975)) 


##############################################################################
### allopatry, same host: Qv-T and Qv-F ###
Mating.Tnew2<-Mating.Tnew[Mating.Tnew$HPM=="Qv" & Mating.Tnew$HPF=="Qv",]
M.m2<-glmer(Mt~sp.pair+(1|Source..M.)+(1|Source..F.),data=Mating.Tnew[Mating.Tnew$HPM=="Qv" & Mating.Tnew$HPF=="Qv",],
            family="binomial")
summary(M.m2)

Mating.Tnew2.lsm<-data.frame(lsmeans(M.m2,~sp.pair))
Mating.Tnew2.lsm$lsmean<-exp(Mating.Tnew2.lsm$lsmean)/(exp(Mating.Tnew2.lsm$lsmean)+1)
Mating.Tnew2.lsm$asymp.LCL<-exp(Mating.Tnew2.lsm$asymp.LCL)/(exp(Mating.Tnew2.lsm$asymp.LCL)+1)
Mating.Tnew2.lsm$asymp.UCL<-exp(Mating.Tnew2.lsm$asymp.UCL)/(exp(Mating.Tnew2.lsm$asymp.UCL)+1)


ddply(Mating.Tnew2,c("sp.pair"),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew2.sum<-ddply(Mating.Tnew2,c("HP.M.","HP.F."),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew2.sum$Mtcount<-Mating.Tnew2.sum$Mt*Mating.Tnew2.sum$count
Mating.Tnew2.sum

Mating.Tnew2.sum1<-ddply(Mating.Tnew2,c("HP.M.","HP.F.","Source..M.","Source..F."),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew2.sum1$pair<-"same"
Mating.Tnew2.sum1$pair[Mating.Tnew2.sum1$HP.M.=="Qv-F" & Mating.Tnew2.sum1$HP.F.=="Qv-T"]<-"different"
Mating.Tnew2.sum1$pair[Mating.Tnew2.sum1$HP.M.=="Qv-T" & Mating.Tnew2.sum1$HP.F.=="Qv-F"]<-"different"
Mating.Tnew2.sum1<-Mating.Tnew2.sum1[Mating.Tnew2.sum1$count>2,]

### Qv-T and Qv-F SI with bootstrap ###
Qvt<-c(rep(0,104),rep(1,51))
Qvf<-c(rep(0,38),rep(1,9))
Qvtf<-c(rep(0,61),rep(1,21))
Qvft<-c(rep(0,30),rep(1,7))
home.fre<-(sum(Qvt)/length(Qvt)+sum(Qvf)/length(Qvf))/2
heto.fre<-(sum(Qvtf)/length(Qvtf)+sum(Qvft)/length(Qvft))/2
TFSI.mean<-1-2*heto.fre/(home.fre+heto.fre)
TF.SI<-c()
for (i in 1:10000){
  Qvt.new<-sample(Qvt,replace=TRUE)
  Qvf.new<-sample(Qvf,replace=TRUE)
  Qvtf.new<-sample(Qvtf,replace=TRUE)
  Qvft.new<-sample(Qvft,replace=TRUE)
  home.fre<-(sum(Qvt.new)/length(Qvt.new)+sum(Qvf.new)/length(Qvf.new))/2
  heto.fre<-(sum(Qvtf.new)/length(Qvtf.new)+sum(Qvft.new)/length(Qvft.new))/2
  TF.SI[i]<-1-2*heto.fre/(home.fre+heto.fre)
}

TF.SI
quantile(TF.SI, c(0.025,0.975)) 
mean(TF.SI)

unique(Mating.Tnew2$Source..M.[Mating.Tnew2$HP.M.=="Qv-T"]) ## "Rice" "GM"   "PY"  ###
unique(Mating.Tnew2$Source..M.[Mating.Tnew2$HP.M.=="Qv-F"]) ## "KRE"  "Okee" "Alva" "JI" ###

Mating.Tnew2.sum<-ddply(Mating.Tnew2,c("HP.M.","HP.F.","Source..M.","Source..F."),summarize,
                        Mt=mean(Mt),count=sum(count))

Mating.Tnew2.sum$Mt.count<-Mating.Tnew2.sum$Mt*Mating.Tnew2.sum$count

Mating.Tnew2.sum

#### allopatry, different host: Qv-T and Qg ####
Mating.Tnew3<-Mating.Tnew[Mating.Tnew$HP.M.!="Qv-F" & Mating.Tnew$HP.F.!="Qv-F",]
M.m3<-glmer(Mt~sp.pair+(1|Source..M.)+(1|Source..F.),data=Mating.Tnew[Mating.Tnew$HP.M.!="Qv-F" & Mating.Tnew$HP.F.!="Qv-F",],
            family="binomial")
summary(M.m3)
Mating.Tnew3.lsm<-data.frame(lsmeans(M.m3,~sp.pair))
Mating.Tnew3.lsm$lsmean<-exp(Mating.Tnew3.lsm$lsmean)/(exp(Mating.Tnew3.lsm$lsmean)+1)
Mating.Tnew3.lsm$asymp.LCL<-exp(Mating.Tnew3.lsm$asymp.LCL)/(exp(Mating.Tnew3.lsm$asymp.LCL)+1)
Mating.Tnew3.lsm$asymp.UCL<-exp(Mating.Tnew3.lsm$asymp.UCL)/(exp(Mating.Tnew3.lsm$asymp.UCL)+1)

Mating.Tnew3.sum1<-ddply(Mating.Tnew3,c("HP.M.","HP.F.","Source..M.","Source..F."),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew3.sum1$pair<-"same"
Mating.Tnew3.sum1$pair[Mating.Tnew3.sum1$HP.M.=="Qg" & Mating.Tnew3.sum1$HP.F.=="Qv-T"]<-"different"
Mating.Tnew3.sum1$pair[Mating.Tnew3.sum1$HP.M.=="Qv-T" & Mating.Tnew3.sum1$HP.F.=="Qg"]<-"different"
Mating.Tnew3.sum1<-Mating.Tnew3.sum1[Mating.Tnew3.sum1$count>2,]

### Qv-T and Qg SI with bootstrap ###
Mating.Tnew3.sum<-ddply(Mating.Tnew3,c("HP.M.","HP.F."),summarize,Mt=mean(Mt),count=sum(count))
Mating.Tnew3.sum$Mtcount<-Mating.Tnew3.sum$Mt*Mating.Tnew3.sum$count
Mating.Tnew3.sum

Qg<-c(rep(0,178),rep(1,38))
Qvt<-c(rep(0,104),rep(1,51))
Qvtg<-c(rep(0,191),rep(1,87))
Qvgt<-c(rep(0,133),rep(1,10))
home.fre<-(sum(Qvt)/length(Qvt)+sum(Qg)/length(Qg))/2
heto.fre<-(sum(Qvtg)/length(Qvtg)+sum(Qvgt)/length(Qvgt))/2
TGSI.mean<-1-2*heto.fre/(home.fre+heto.fre)
TG.SI<-c()
for (i in 1:10000){
  Qg.new<-sample(Qg,replace=TRUE)
  Qvt.new<-sample(Qvt,replace=TRUE)
  Qvtg.new<-sample(Qvtg,replace=TRUE)
  Qvgt.new<-sample(Qvgt,replace=TRUE)
  home.fre<-(sum(Qvt.new)/length(Qvt.new)+sum(Qg.new)/length(Qg.new))/2
  heto.fre<-(sum(Qvtg.new)/length(Qvtg.new)+sum(Qvgt.new)/length(Qvgt.new))/2
  TG.SI[i]<-1-2*heto.fre/(home.fre+heto.fre)
}

c(quantile(TG.SI, c(0.025,0.975)))
c(quantile(TF.SI, c(0.025,0.975)))
c(quantile(FG.SI, c(0.025,0.975)))
t.test(TF.SI,TG.SI)
t.test(TF.SI,FG.SI)
t.test(TG.SI,FG.SI)

#############################################################################################
######## making figure 3A #########
Mating.Tnew1.lsm$pair<-"B. treatae x B. fossoria"
Mating.Tnew2.lsm$pair<-"B. kinseyi x B. treatae"
Mating.Tnew3.lsm$pair<-"B. kinseyi x B. fossoria"

Mating.Tlsm<-rbind(Mating.Tnew1.lsm,Mating.Tnew2.lsm,Mating.Tnew3.lsm)
Mating.Tlsm$sp.pair<-factor(Mating.Tlsm$sp.pair,levels=c("same","different"))
Mating.Tlsm$pair<-factor(Mating.Tlsm$pair,levels=c("B. kinseyi x B. treatae","B. kinseyi x B. fossoria",
                                                   "B. treatae x B. fossoria"))
Mating.Tnew1.sum1$type<-"B. treatae x B. fossoria"
Mating.Tnew2.sum1$type<-"B. kinseyi x B. treatae"
Mating.Tnew3.sum1$type<-"B. kinseyi x B. fossoria"

Mating.Tnew.sum<-rbind(Mating.Tnew1.sum1,Mating.Tnew2.sum1,Mating.Tnew3.sum1)
Mating.Tnew.sum$type<-factor(Mating.Tnew.sum$type,levels=c("B. kinseyi x B. treatae","B. kinseyi x B. fossoria",
                                                           "B. treatae x B. fossoria"))
Mating.Tnew.sum$pair<-factor(Mating.Tnew.sum$pair,levels=c("same","different"))
Mating.Tlsm<-Mating.Tlsm[order(Mating.Tlsm$pair,Mating.Tlsm$sp.pair),]
Mating.Tnew.sum<-Mating.Tnew.sum[order(Mating.Tnew.sum$type,Mating.Tnew.sum$pair),]

Mate.sumcount<-ddply(Mating.Tnew,c("HP.M.","HP.F."),summarize,TN=sum(count))
Mate.sumcount$species.M[Mate.sumcount$HP.M.=="Qg"]<-"B. fossoria"
Mate.sumcount$species.M[Mate.sumcount$HP.M.=="Qv-F"]<-"B. treatae"
Mate.sumcount$species.M[Mate.sumcount$HP.M.=="Qv-T"]<-"B. kinseyi"

Mate.sumcount$species.F[Mate.sumcount$HP.F.=="Qg"]<-"B. fossoria"
Mate.sumcount$species.F[Mate.sumcount$HP.F.=="Qv-F"]<-"B. treatae"
Mate.sumcount$species.F[Mate.sumcount$HP.F.=="Qv-T"]<-"B. kinseyi"

Mating.Tlsm$TN<-NA
Mating.Tlsm$TN[1]<-202
Mating.Tlsm$TN[2]<-119
Mating.Tlsm$TN[3]<-216+115
Mating.Tlsm$TN[4]<-143+278
Mating.Tlsm$TN[5]<-216+47
Mating.Tlsm$TN[6]<-96+71

Mating.Tlsm$sd<-Mating.Tlsm$SE/sqrt(Mating.Tlsm$TN)
Mating.Tlsm$sp.pair<-factor(Mating.Tlsm$sp.pair,levels=c("same","different"))
Mating.Tlsm$pairs<-paste(Mating.Tlsm$pair,Mating.Tlsm$sp.pair,sep="_")
Mating.Tlsm$pairs<-factor(Mating.Tlsm$pairs,levels=c("B. kinseyi x B. treatae_same", "B. kinseyi x B. treatae_different",
                                                     "B. kinseyi x B. fossoria_same", "B. kinseyi x B. fossoria_different",
                                                     "B. treatae x B. fossoria_same", "B. treatae x B. fossoria_different"))

Mating.Tnew.sum$pairs<-paste(Mating.Tnew.sum$type,Mating.Tnew.sum$pair,sep="_")
colnames(Mating.Tnew.sum)[7]<-"sp.pair"
colnames(Mating.Tnew.sum)[8]<-"pair"

Mate.pre.fig<-ggplot()+
  geom_bar(data = Mating.Tlsm,aes(y=lsmean,x=pair,fill=pairs),
           stat="identity",position=position_dodge(width=0.8),width=0.7)+
  geom_errorbar(data = Mating.Tlsm,aes(ymin=lsmean-sd,ymax=lsmean+sd,x=pair,color=pairs),position=position_dodge(width=0.8), 
                stat="identity", width=.01)+
  scale_color_manual(values=c("black","black","black","black","black","black"))+
  scale_fill_manual(values=c("darkgreen","mediumspringgreen","blue3","lightblue","orange4","burlywood1"))+
  geom_point(data =Mating.Tnew.sum, aes(y=Mt,x=pair,color=pairs), position=position_dodge(width=0.8), 
             alpha =0.5,size=3)+
  scale_y_continuous(expand=c(0,0),limits=c(-0.002,0.8),breaks=seq(0,1,0.1))+
  ylab("Mating probability")+xlab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(colour="black"),
        axis.text.y = element_text(size=24),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=28,margin=margin(0,20,0,0)),
        plot.margin=unit(c(1,0,3,0),"lines"),
        axis.ticks.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")+
  geom_segment(aes(y=c(0.5,0.45,0.501,0.451,0.501,0.451),yend=c(0.5,0.45,0.48,0.43,0.48,0.43),
                   x=c(1.8,2.8,1.8,2.8,2.2,3.2),xend=c(2.2,3.2,1.8,2.8,2.2,3.2)), size=1)+
  annotate("text",x=c(2,3),y=c(0.51,0.46),label=c("*","*"),size=20)

Mate.pre.fig

Mate.pre.fig1<-ggdraw(add_sub(Mate.pre.fig,c(expression(paste(italic("Bk")," x ",italic("Bk"))),expression(paste(italic("Bk")," x ",italic("Bt"))),
                                             expression(paste(italic("Bk")," x ",italic("Bk"))),expression(paste(italic("Bk")," x ",italic("Bf"))),
                                             expression(paste(italic("Bt")," x ",italic("Bt"))),expression(paste(italic("Bt")," x ",italic("Bf")))), 
                              x=c(0.12,0.25,0.44,0.56,0.75,0.875),size=24,color="black"))

Mate.pre.fig1

Mate.pre.fig2<-ggdraw(add_sub(Mate.pre.fig1,c(expression(paste("& ",italic("Bt")," x ",italic("Bt"))),
                                              expression(paste("& ",italic("Bf")," x ",italic("Bf"))),
                                              expression(paste("& ",italic("Bf")," x ",italic("Bf")))),
                              x=c(0.17,0.46,0.76),y=1.8,size=24,color="black"))

Mate.pre.fig2

Mate.pre.fig3<-ggdraw(add_sub(Mate.pre.fig2,c("Host:","Same","Same","Same","Different","Same","Different"), 
                              x=c(0.1,0.18,0.3,0.47,0.59,0.77,0.88),y=1.4,size=26,color="black"))

Mate.pre.fig3

######## making figure 3B:  bootstrap sexual isolation comparison#########
SI<-data.frame(matrix(nrow=3,ncol=4))
colnames(SI)<-c("pair","mean.SI","lower.SI","upper.SI")
SI[,1]<-c("Qv-T_Qv-F","Qv-T_Qg","Qv-F_Qg")
SI[,2]<-c(TFSI.mean,TGSI.mean,FGSI.mean)
SI[1,3:4]<-c(quantile(TF.SI, c(0.025,0.975)))
SI[2,3:4]<-c(quantile(TG.SI, c(0.025,0.975)))
SI[3,3:4]<-c(quantile(FG.SI, c(0.025,0.975)))

SI$pair<-factor(SI$pair,levels=c("Qv-T_Qv-F","Qv-T_Qg","Qv-F_Qg"))

SI_plot<-ggplot(data=SI,aes(x=pair,y=mean.SI,color=pair,fill=pair))+
  geom_errorbar(aes(ymin=lower.SI,ymax=upper.SI),width=0.08)+
  geom_point(size=6)+
  geom_hline(yintercept=0,linetype="dashed")+
  scale_color_manual(values=c("green4","blue3","burlywood3"))+
  scale_y_continuous(expand=c(0,0),limits=c(-0.2,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),  
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=24,colour="black"),
        axis.title=element_text(size=28),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.margin=unit(c(1,0,3,0),"lines"),
        plot.title = element_text(size = 28, face = "bold",hjust = 0.5))+
  labs(x=" ",y="Sexual isolation")

SI_plot

SI_plot1<-ggdraw(add_sub(SI_plot,c(expression(paste(italic("B. kinseyi")," x ", italic("B. treatae"))),
                                   expression(paste(italic("B. kinseyi")," x ", italic("B. fossoria"))),
                                   expression(paste(italic("B. treatae")," x ", italic("B. fossoria")))), 
                         x=c(0.19,0.51,0.83),size=22,y=0.5,color="black"))
SI_plot1

SI_plot2<-ggdraw(add_sub(SI_plot1,c("Host plant: ","Same","Different","Different"), 
                         x=c(0.1,0.25,0.54,0.83),y=1.5,size=22,color="black"))

SI_plot2

SI_plot3<-ggdraw(add_sub(SI_plot2,c("Geography: ","Allopatry","Allopatry","Sympatry"), 
                         x=c(0.11,0.25,0.54,0.83),y=1.2,size=22,color="black"))
SI_plot3

plot_grid(Mate.pre.fig2,SI_plot3,nrow=2,labels=c("A","B"),label_size = 38) ### save as 12*18 ###


###########################################################################
######## male and female from Qv and Qg, test the mate preference for male and female#######
###########################################################################
Mating.T$Block.id<-paste(Mating.T$Year,Mating.T$Block.ID,sep="_")
Mating.T$Year<-as.factor(Mating.T$Year)
Mating.T.sym<-Mating.T[Mating.T$pair3 %in% c("Qg_Qg","Qv-F_Qg","Qv-F_Qv-F"),]
#########################################################################
### male wing buzz ######
### statistical analysis ###
Wp.m1<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
               data=Mating.T,
               family="binomial")
summary(Wp.m1)
Wp.m1<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
               data=Mating.T.sym,
               family="binomial")
summary(Wp.m1)
coef(summary(Wp.m1))$cond[6,4]
fixef(Wp.m1)
simulationOutput <- simulateResiduals(fittedModel = Wp.m1)
plot(simulationOutput)
plotResiduals(simulationOutput,Mating.T.sym$Year )
plotResiduals(simulationOutput,Mating.T.sym$HP.M. )
plotResiduals(simulationOutput,Mating.T.sym$pair4)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
lsmeans(Wp.m1,pairwise~HP.M.*pair4,adjust="tukey")

Wp.m1.anova<-data.frame(Anova(Wp.m1,type="III"))
### power analysis ####
Mating.T.sym.Qv<-Mating.T.sym[Mating.T.sym$HP.M.!="Qg",]

Qg.same<-Mating.T.sym[Mating.T.sym$HP.M.=="Qg" & Mating.T.sym$pair4=="same",]
Qg.diff<-Mating.T.sym[Mating.T.sym$HP.M.=="Qg" & Mating.T.sym$pair4=="different",]
same.no<-length(Qg.same[,1])
diff.no<-length(Qg.diff[,1])
z.list<-c()
p.list<-c()
p.list2<-c()
for (i in 1:1000){
  Qg.same.temp<-Qg.same[sample(1:same.no,47),]
  Qg.diff.temp<-Qg.diff[sample(1:diff.no,71),]
  Mating.T.sym.temp<-rbind(Mating.T.sym.Qv,Qg.same.temp,Qg.diff.temp)
  Wp.mnew<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
                   data=Mating.T.sym.temp,family="binomial")
  p.list2[i]<-data.frame(lsmeans(Wp.mnew,pairwise~HP.M.*pair4,adjust="tukey")$contrast)[5,6]
  p.list[i]<-coef(summary(Wp.mnew))$cond[6,4]
  z.list[i]<-coef(summary(Wp.mnew))$cond[6,3]
}
hist(z.list)
hist(p.list2)
length(p.list[p.list<0.05])
length(p.list[p.list<0.1])
length(p.list2[p.list2<0.05])
length(z.list[z.list>0])
##### Total Wing buzz period m2 ####
Wp.m2<-glmmTMB(T.W~Year+HP.M.*pair4+(1|Source..M.)+(1|Source..F.),
               data=Mating.T.sym,
               family="nbinom2")
summary(Wp.m2)
simulationOutput<- simulateResiduals(fittedModel = Wp.m2)
plot(simulationOutput)
plotResiduals(simulationOutput,Mating.T.sym$Year )
plotResiduals(simulationOutput,Mating.T.sym$HP.M. )
plotResiduals(simulationOutput,Mating.T.sym$pair4)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)
lsmeans(Wp.m2,pairwise~HP.M.*pair4,adjust="tukey")

Wp.m2.anova<-data.frame(Anova(Wp.m2,type="III"))
### Qv wing buzz plot1 ######
### least square means of wing buzz frequency ####
W.lsm<-lsmeans(Wp.m1,~HP.M.*pair4)
W.lsm<-data.frame(summary(W.lsm,type="response"))
W.lsm$HP.M.<-c("B. fossoria","B. treatae","B. fossoria","B. treatae")
W.lsm$HP.M.<-factor(W.lsm$HP.M.,levels=c("B. treatae","B. fossoria"))
W.lsm$pair4<-factor(W.lsm$pair4,levels=c("same","different"))
W.lsm<-W.lsm[order(W.lsm$HP.M.,W.lsm$pair4,decreasing=FALSE),]
### adding number of individuals per group ###
W.count<-ddply(Mating.T.sym,c("HP.M.","pair4"),summarize,count=sum(count))
W.count$HP.M.<-factor(W.count$HP.M.,levels=c("Qv-F","Qg"))
W.count$pair4<-factor(W.count$pair4,levels=c("same","different"))
W.count<-W.count[order(W.count$HP.M.,W.count$pair4,decreasing = FALSE),]
W.lsm$count<-W.count$count
#########################################
W.lsm$pair<-"conspecific"
W.lsm$pair[W.lsm$pair4=="different"]<-"heterospecific"
W.lsm$comp<-paste(W.lsm$HP.M.,W.lsm$pair4,sep="_")
W.lsm$HP.F.<-c("B. treatae","B. fossoria", "B. fossoria", "B. treatae")
W.lsm$HP.F.<-factor(W.lsm$HP.F.,levels=c("B. treatae","B. fossoria"))
W.lsm<-W.lsm[order(W.lsm$HP.M.,W.lsm$HP.F.),]
#################################################
Wing.plot<-ggplot(data=W.lsm,aes(x=HP.M.,group=HP.F.))+
  stat_summary(aes(pattern=pair,y=prob,fill=HP.M.,pattern_fill=HP.M.),
               fun = "mean", position=position_dodge(width=0.8),width=0.6,pattern_density=0.5,
               geom = "bar_pattern", colour="black")+
  scale_pattern_manual(values = c(conspecific = "none", heterospecific = "stripe"))+
  scale_fill_manual(values=c("blue3","orange2"))+
  scale_pattern_fill_manual(values=c("orange2","blue3"))+
  geom_errorbar(aes(ymin=prob-SE,ymax=prob+SE),position=position_dodge(width=0.8),width=0.2)+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=22,colour="black"),
    axis.title=element_text(size=26),
    axis.title.y=element_text(margin=margin(0,20,0,0)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,0,3,0),"lines"))+
  labs(x=" ",y="Male mate preference")+
  annotate("text",x=c(0.78,1.23,1.79,2.23),y=W.lsm$prob+W.lsm$SE+0.04,label=paste("n= ",W.lsm$count,sep=""),size=9)+
  geom_segment(aes(x=1.8,y=0.53,xend=2.2,yend=0.53),size=1)+annotate("text",x=2,y=0.56,label="**",size=11)

Wing.plot
### annotate("text",x=1.5,y=0.7,label="Male species x Species pair: ",size=9)+
### annotate("text",x=1.5,y=0.64,label=expression(paste(italic("z"),"= 2.664, ",italic("p")," =0.008")),size=9)+
### geom_segment(aes(x=1.8,y=0.56,xend=2.2,yend=0.56),size=1)+annotate("text",x=2,y=0.58,label="**",size=11)

Wing.plot1<-ggdraw(add_sub(Wing.plot,c("Female:",expression(italic("B. treatae")),expression(italic("B. fossoria")),
                                       expression(italic("B. treatae")),expression(italic("B. fossoria"))),x=c(0.03,0.19,0.37,0.63,0.81),y=0.4,size=20))
Wing.plot1
Wing.plot2<-ggdraw(add_sub(Wing.plot1,c("Male:",expression(italic("B. treatae")),expression(italic("B. fossoria"))),
                           x=c(0.17,0.38,0.77),y=1.6,size=20))
Wing.plot2
boxes <- data.frame( x = c(0.24,0.62),y = c(0.1,0.10))
Wing.plot3<-ggdraw(Wing.plot2) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +0.28, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
Wing.plot3 
Wing.plot4<-ggdraw(add_sub(Wing.plot3,c("Higher dispersal","Lower dispersal"),
                           x=c(0.38,0.77),y=1.6,size=20))
Wing.plot4
### save as 10*8 ###
###########################
###########################################################################
#### female copulation frequency ####
F.Mating.sym<-Mating.T.sym[Mating.T.sym$W==1,]
F.Mating.sym$Fa<-F.Mating.sym$Mo+F.Mating.sym$Mt
F.Mating.sym$Fa[F.Mating.sym$Fa!=0]<-1
######## male courtship behavior: including wing buzz and mounting ####
####### female mating behavior: copulation frequency when  the male courtship behavior presents ###
#### stastistical analysis ####
str(F.Mating.sym$Year)
Fsym.m1<-glmmTMB(Fa~Year+HP.F.*pair4+(1|Source..F.)+(1|Source..F.),data=F.Mating.sym,
                 family=binomial)
summary(Fsym.m1)

simulationOutput <- simulateResiduals(fittedModel = Fsym.m1)
plot(simulationOutput)
plotResiduals(simulationOutput,F.Mating.sym$Year )
plotResiduals(simulationOutput,F.Mating.sym$HP.F. )
plotResiduals(simulationOutput,F.Mating.sym$pair4)
testDispersion(simulationOutput)
testZeroInflation(simulationOutput)

lsmeans(Fsym.m1,pairwise~HP.F.*pair4,adjust="tukey")

Fsym.m1.anova<-data.frame(Anova(Fsym.m1,type="III"))
##### least square means of copulation frequency ####
F.lsm<-lsmeans(Fsym.m1,~HP.F.*pair4)
F.lsm<-data.frame(summary(F.lsm,type="response"))
F.lsm$HP.F.<-factor(F.lsm$HP.F.,levels=c("Qv-F","Qg"))
F.lsm$pair4<-factor(F.lsm$pair4,levels=c("same","different"))
F.lsm<-F.lsm[order(F.lsm$HP.F.,F.lsm$pair4,decreasing = FALSE),]
F.lsm$HP.F.<-c("B. treatae","B. treatae","B. fossoria","B. fossoria")
F.lsm$HP.F.<-factor(F.lsm$HP.F.,levels=c("B. treatae","B. fossoria"))
F.lsm$pair<-c("conspecific","heterospecific","conspecific","heterospecific")
##### add the number of individuals #####
F.count<-ddply(F.Mating.sym,c("HP.F.","pair4"),summarize,count=sum(count))

ddply(Mating.T,c("HP.F.","pair4"),summarize,count=sum(count))

F.count$HP.F.<-factor(F.count$HP.F.,levels=c("Qv-F","Qg"))
F.count$pair4<-factor(F.count$pair4,levels=c("same","different"))
F.count<-F.count[order(F.count$HP.F.,F.count$pair4,decreasing = FALSE),]
F.lsm$count<-F.count$count
F.lsm$HP.M.<-c("B. treatae","B. fossoria","B. fossoria","B. treatae")
F.lsm$HP.M.<-factor(F.lsm$HP.M.,levels=c("B. treatae","B. fossoria"))
F.lsm<-F.lsm[order(F.lsm$HP.F.,F.lsm$HP.M.),]

F.plot<-ggplot(data=F.lsm,aes(y=prob,x=HP.F.,group=HP.M.))+
  stat_summary(aes(pattern=pair,y=prob,fill=HP.F.,pattern_fill=HP.F.),
               fun = "mean", position=position_dodge(width=0.8),width=0.6,pattern_density=0.5,
               geom = "bar_pattern", colour="black")+
  scale_pattern_manual(values = c(conspecific = "none", heterospecific = "stripe"))+
  scale_fill_manual(values=c("blue3","orange2"))+
  scale_pattern_fill_manual(values=c("orange2","blue3"))+
  geom_errorbar(aes(ymin=prob-SE,ymax=prob+SE),position=position_dodge(width=0.8),width=0.2)+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=22,colour="black"),
    axis.title=element_text(size=26),
    axis.title.y=element_text(margin=margin(0,20,0,0)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,0,3,0),"lines"))+
  labs(x=" ",y="Female mate preference")+
  annotate("text",x=c(0.79,1.24,1.78,2.24),y=F.lsm$prob+F.lsm$SE+0.04,label=paste("n= ",F.lsm$count,sep=""),size=9)
F.plot

F.plot1<-ggdraw(add_sub(F.plot,c("Male:",expression(italic("B. treatae")),expression(italic("B. fossoria")),
                                 expression(italic("B. treatae")),expression(italic("B. fossoria"))),x=c(0.03,0.19,0.37,0.63,0.81),y=0.4,size=20))
F.plot1
F.plot2<-ggdraw(add_sub(F.plot1,c("Female:",expression(italic("B. treatae")),expression(italic("B. fossoria"))),
                        x=c(0.16,0.4,0.76),y=1.6,size=20))
F.plot2
boxes <- data.frame( x = c(0.24,0.62),y = c(0.1,0.1))
F.plot3<-ggdraw(F.plot2) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +0.28, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
F.plot3 
F.plot4<-ggdraw(add_sub(F.plot3,c("Higher dispersal","Lower dispersal"),
                        x=c(0.42,0.77),y=1.6,size=20))
F.plot4

###############################
### combine male and female plot ####
plot_grid(Wing.plot4,F.plot4,ncol=2,labels=c("A","B"),label_size = 38)
#######################################################################
########## test of reinforcement #######################
Mating.T$Block.id<-paste(Mating.T$Year,Mating.T$Block.ID,sep="_")
Mating.T$Year<-as.factor(Mating.T$Year)
#########################################################################
### reinforcement for male mate preference
Mating.Tnew<-Mating.T[Mating.T$pair2 %in% c("Qv-T_Qv-T","Qg_Qg",
                                            "Qv-T_Qg","Qg_Qv-T"),]

Wp.Qv.1<-glmmTMB(W~Year+HP.M.*pair4+(1|Source..M.),
                 data=Mating.Tnew,
                 family="binomial")
summary(Wp.Qv.1)
simulationOutput <- simulateResiduals(fittedModel = Wp.Qv.1)
plot(simulationOutput)
lsmeans(Wp.Qv.1,pairwise~HP.M.*pair4,adjust="tukey")
Wp.Qv.1.anova<-data.frame(Anova(Wp.Qv.1,type="III"))
WQv<-lsmeans(Wp.Qv.1,~HP.M.*pair4)
WQv<-data.frame(summary(WQv,type="response"))
WQv$pair4<-factor(WQv$pair4,levels=c("same","different"))
WQv<-WQv[order(WQv$HP.M.,WQv$pair4,decreasing=FALSE),]
WQv$species<-c("B. fossoria","B. fossoria","B. kinyesi","B. kinyesi")
WQv$species<-factor(WQv$species,levels=c("B. fossoria","B. kinyesi"))
WQv$HP.F.<-c("B. fossoria","B. kinyesi","B. kinyesi","B. fossoria")
WQv$HP.F.<-factor(WQv$HP.F.,levels=c("B. fossoria","B. kinyesi"))


### adding number of individuals per group ###
WQv.count<-ddply(Mating.Tnew,c("HP.M.","pair4"),summarize,count=sum(count))
WQv.count$HP.M.<-factor(WQv.count$HP.M.,levels=c("Qv-F","Qg","Qv-T"))
WQv.count$pair4<-factor(WQv.count$pair4,levels=c("same","different"))
WQv.count<-WQv.count[order(WQv.count$HP.M.,WQv.count$pair4,decreasing = FALSE),]
WQv$count<-WQv.count$count

WQv<-WQv[order(WQv$species,WQv$HP.F.,decreasing=FALSE),]
WQv
#########################################
Wing.Qv.plot<-ggplot(data=WQv,aes(y=prob,x=HP.M.,group=HP.F.))+
  stat_summary(aes(pattern=pair4,y=prob,fill=HP.M.,pattern_fill=HP.M.),
               fun = "mean", position=position_dodge(width=0.8),width=0.6,pattern_density=0.5,
               geom = "bar_pattern", colour="black")+
  scale_pattern_manual(values = c(same = "none", different = "stripe"))+
  scale_pattern_fill_manual(values=c("mediumspringgreen","orange2"))+
  geom_errorbar(aes(ymin=prob-SE,ymax=prob+SE),position=position_dodge(width=0.8),width=0.2)+
  scale_fill_manual(values=c("orange2","mediumspringgreen"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,1),breaks=seq(0,1,0.2))+
  theme_bw()+
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),  
    panel.background=element_blank(),
    axis.line=element_line(colour="black"),
    axis.text=element_text(size=22,colour="black"),
    axis.title=element_text(size=26),
    axis.title.y=element_text(margin=margin(0,20,0,0)),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    plot.margin=unit(c(1,0,3,0),"lines"))+
  labs(x=" ",y="Male mate preference")+
  annotate("text",x=c(0.76,1.24,1.76,2.24),y=WQv$prob+WQv$SE+0.06,label=paste("n= ",WQv$count,sep=""),size=7)+
  geom_segment(aes(x=0.8,y=0.56,xend=1.2,yend=0.56),size=1)+annotate("text",x=1,y=0.58,label="***",size=12)
Wing.Qv.plot

Wing.Qv.plot1<-ggdraw(add_sub(Wing.Qv.plot,c("Female:",expression(italic("B. fossoria")),expression(italic("B. kinseyi")),
                                             expression(italic("B. fossoria")),expression(italic("B. kinseyi"))),x=c(0.03,0.19,0.37,0.63,0.81),y=0.4,size=20))
Wing.Qv.plot1
Wing.Qv.plot2<-ggdraw(add_sub(Wing.Qv.plot1,c("Male:",expression(paste(italic("B. fossoria")," (sympatry)")),expression(paste(italic("B. kinseyi")," (allopatry)"))),
                              x=c(0.14,0.37,0.76),y=1.6,size=20))
Wing.Qv.plot2
boxes <- data.frame( x = c(0.22,0.61),y = c(0.12,0.12))
Wing.Qv.plot3<-ggdraw(Wing.Qv.plot2) + 
  geom_rect(data = boxes, aes(xmin = x, xmax = x +0.29, ymin = y, ymax = y +0.0018),colour = "black", fill = "black")
Wing.Qv.plot3


### save image 12.89*8.7 #######
###########################################################################
