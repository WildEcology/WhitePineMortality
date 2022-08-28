
library(tidyverse)
library(magrittr)
library(ggeffects)
library(sjPlot)


mort=read_csv("data_mortality_July16_19.csv") ##newer version ##Lots of -99!!
dat=read.csv("model_data_July3_2019.csv")
elevation=dat%>%
  dplyr::select(plot, elevation)%>%
  group_by(plot)%>%
  dplyr::summarize(elevation=paste(unique(elevation)))%>%
  mutate(elevation=as.numeric(elevation))

vpd=read_csv("meanvpd_dewpt.csv")
drought=read_csv("dewpointvpd.csv")

##drought years
drought=drought%>%
  filter(year>2011&year<2016)%>%
  group_by(plot)%>%
  dplyr::summarize(droughtvpd=mean(vpdmax), droughtdew=mean(dewpoint))

##mort dat correct dim
mortlate=mort%>%
  filter(year>2000)%>%
  dplyr::select(plot, species, tree_id, year, status, dbh)%>%
  set_colnames(c("plot", "species", "tree_id", "year00", "status00", "dbh00"))

mort_sm=mort%>%
  filter(year<2000)%>%
  dplyr::select(plot, species, tree_id, year, status, dbh)%>%
  full_join(mortlate)

dat2=dat%>%
  dplyr::select(plot, tree_id, fire_only, slope, ppt_post95,  strm_dist,
                species, beetles, south, baplot95,  era95_inc_tot,tmin_post95, 
                inc_tot,  tmax_post95,tmin_pre95,lake_dist)%>%
  mutate(species=as.character(species))


combdat=mort_sm%>%
  left_join(dat2, by=c("plot", "tree_id", "species"))%>%
  mutate(replace(status, status==-99, 0))%>%
  mutate(mort=ifelse(status00==1, 0, 1))%>%
  left_join(drought)%>%
  left_join(elevation, by="plot")%>%
  left_join(vpd)

summary(combdat$elevation)


##adding elevational bands every ~2500ft
combdat2=combdat%>%
  #filter(species!="PIBA")%>%
  mutate(elevband=ifelse(elevation<7000, 7000, 
                         ifelse(elevation<9500, 9500, 11000)))



##================================= all mortality===========================================
allmort=combdat2%>%
  filter(year00>2014) %>% 
  dplyr::select(fire_only ,dbh, elevation,era95_inc_tot , inc_tot, strm_dist, lake_dist , 
                ppt_post95,droughtvpd, droughtdew, beetles , south , slope ,baplot95 ,
                tmin_post95,tmin_pre95,plot, mort, tree_id, year, species, status, meantvpd16 )%>%
  mutate_at(scale, .vars = vars(-plot, -mort,-fire_only,-era95_inc_tot,-beetles, -tree_id, -year, -species, -status, -meantvpd16))%>%
  as.data.frame(.)

#write.csv(allmort, "robertsdat.csv")




allmortmod=glmer(mort ~ fire_only * meantvpd16 + era95_inc_tot*meantvpd16 + beetles* meantvpd16+
                   fire_only*beetles + fire_only*era95_inc_tot + beetles*era95_inc_tot + 
                   (1|plot/species), control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                            optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                 family=binomial, data=allmort)

allmortmod2=glmer(mort ~ dbh+fire_only*beetles + fire_only*era95_inc_tot + beetles*era95_inc_tot + 
                    fire_only*beetles*meantvpd16 + fire_only*era95_inc_tot*meantvpd16 + beetles*era95_inc_tot*meantvpd16 + 
                    (1|plot/species), control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                             optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                  family=binomial, data=allmort)

summary(allmortmod2)

tab_model(allmortmod2, show.est=T, show.ci = F, show.se = T)


allmortmod_long=glmer(mort ~ fire_only * meantvpd16 + era95_inc_tot*meantvpd16 + beetles* meantvpd16+
                        fire_only*beetles + fire_only*era95_inc_tot + beetles*era95_inc_tot + beetles*era95_inc_tot*fire_only+
                        fire_only*beetles*meantvpd16+fire_only*era95_inc_tot*meantvpd16+
                        (1|plot/species), control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                                 optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                      family=binomial, data=allmort)


allmortmod_long=glmer(mort ~ fire_only * era95_inc_tot * beetles + meantvpd16+
                        (1|plot/species), control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                                 optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                      family=binomial, data=allmort)



summary(allmortmod_long)
AIC(allmortmod, allmortmod_long)


predplot=ggpredict(allmortmod2, terms = c("fire_only[all]","era95_inc_tot[all]"))
plot(predplot)
predplot=ggpredict(allmortmod2, terms = c("meantvpd16", "beetles[all]","era95_inc_tot[all]"))
plot(predplot)



##transforming to get positive odds
library(sjstats)
odds=odds_to_rr(allmortmod) 

odds2=as.data.frame(odds)
names(odds2)[1:3]=c( "estimate", "lower", "upper")
odds2$var=c("intercept", "Fire-only","Rust'95", "Beetles", "Beetles*Fire",
            "Rust*fire","Rust*Beetles", "Rust*fire*beetles")

odds3=odds2%>%
  mutate(pvalue=c(.001, .001,.001, .001,.001, .001, .001,.001,.001, 1))%>%
  mutate(upper2=upper-estimate, lower2=estimate-lower)%>%
  filter(var!="intercept")


oddsgraph<-ggplot(odds3, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  scale_fill_manual(values = c("#E69F00","palevioletred","#56B4E9", "orchid4",
                               "cadetblue4", "lightsalmon2", "darkred", "blue", "red")) +
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.01, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  ggtitle(label="Drivers of mortality and their interactions")+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 3.5)+
  ylab("Transformed logistic regression estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=15), plot.title = element_text(hjust = .5, size=18))

oddsgraph


library(effects)
e=allEffects(allmortmod)
print(e)
plot(e)

cols <- c("Fire"="#E69F00","Rust"="palevioletred", "Beetles"="#56B4E9", 
          "Beetles*fire"="cadetblue4", "Rust*fire"="lightsalmon2",
          "Rust*Beetles"="orchid4", "Rust*fire*beetles"="darkred")
#"Mintemps"="slateblue", "DBH"= "mediumturquoise")

intdat=data.frame(summary(allmortmod)$coefficients)
intdat=round(intdat, digits=3)
names(intdat)[1:4]=c("estimate", "sterr", "zscore", "pvalue")
intdat$var=c("intercept", "Fire", "Rust", "Beetles", 
             "Beetles*fire","Rust*fire", "Rust*Beetles", "Rust*fire*beetles")
intdat=intdat%>%
  filter(var!="intercept")


intgraph<-ggplot(intdat, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  scale_fill_manual(values = cols)+
  scale_x_discrete(limits=c( "Fire", "Rust", "Beetles", 
                             "Beetles*fire","Rust*fire", "Rust*Beetles", 
                             "Rust*fire*beetles"),
                   labels = c("Fire", "Rust", "Beetles", 
                              "Beetles*fire","Rust*fire", 
                              "Rust*Beetles", "Rust*fire*beetles"))+
  geom_errorbar(aes(ymin=estimate-sterr, ymax=estimate+sterr), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.01, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  #ggtitle(label="Drivers of mortality and their interactions")+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(-.3, 2.5)+
  ylab("Logistic regression estimates of mortality")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=15), plot.title = element_text(hjust = .5, size=15))

intgraph


beetint=glmer(mort~beetles+(beetles:era95_inc_tot)+(beetles:fire_only)+(1|plot/species), 
              control=glmerControl(optimizer="bobyqa"),family=binomial, data=allmort)

summary(beetint)

allmortsqrt=glmer(mort ~  I(sqrt(era95_inc_tot+1*fire_only+1))+ I(sqrt(beetles+1*era95_inc_tot+1))+
                    I(sqrt(beetles+1*fire_only+1))+(1|plot), 
                  control=glmerControl(optimizer="bobyqa"),family=binomial, data=allmort)
summary(allmortsqrt)


##====================================ELEVATION BANDS==========================================
library(optimx) ##this optimizer seems to work well on these complex models

low=subset(combdat2, elevband==7000)
med=subset(combdat2, elevband==9500)
high=subset(combdat2, elevband==11000)
high=high%>%
  mutate(fire_only=replace_na(fire_only, 0))


length(unique(high$tree_id))
summary(unique(high$tmin_post95))
summary(unique(med$tmin_post95))
summary(unique(low$tmin_post95))


##low elevation mortality
lowmort=low%>% 
  dplyr::select(fire_only ,elevation,dbh, era95_inc_tot , inc_tot, strm_dist, lake_dist , ppt_post95,droughtvpd, droughtdew,
                beetles , south , slope ,baplot95 ,tmin_post95,plot, mort, tree_id, year, species, status )%>%
  mutate_at(scale, .vars = vars(-plot, -mort, -tree_id, -year, -species, -status,-fire_only,-beetles, -era95_inc_tot, -inc_tot))%>%
  as.data.frame(.)

lowmortmodall=glmer(mort ~ fire_only + era95_inc_tot + beetles + tmin_post95 + baplot95+
                      fire_only*beetles + era95_inc_tot*fire_only + beetles*era95_inc_tot+
                      tmin_post95*beetles + tmin_post95*era95_inc_tot + tmin_post95*fire_only+
                      (1|plot), control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                                       optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                    family=binomial, data=lowmort)
summary(lowmortmod)
lowmortint=glmer(mort ~ fire_only + era95_inc_tot + beetles + 
                   fire_only*beetles + era95_inc_tot*fire_only + beetles*era95_inc_tot+
                   fire_only*beetles*era95_inc_tot+(1|plot), 
                 control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                        optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                 family=binomial, data=lowmort)


summary(lowmortint)

oddslow=odds_to_rr(allmortmod) 

oddslow2=as.data.frame(oddslow)
names(oddslow2)[1:3]=c( "estimate", "lower", "upper")
oddslow2$var=c("intercept", "Fire-only","Rust'95", "Beetles",  "Beetles*fire",
               "Rust*fire","Rust*Beetles", "Rust*fire*beetles")

oddslow3=oddslow2%>%
  mutate(pvalue=c(.001, .001,.001, .001, .001,.002,.002, 1))%>%
  mutate(upper2=upper-estimate, lower2=estimate-lower)%>%
  filter(var!="intercept")%>%
  mutate(type=c("direct", "direct", "direct", "ind", "ind", "ind","ind"))%>%
  mutate(type2=c(1,1,1,2,2,2,2))

colsdir <- c("Fire-only"="#E69F00","Rust'95"="palevioletred", "Beetles"="#56B4E9")
colsind=c("Beetles*fire"="orchid4","Rust*fire"="cadetblue4","Rust*Beetles"="lightsalmon2", "Rust*fire*beetles"="darkred")

oddslowgraph<-ggplot(oddslow3, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  scale_fill_manual(values = cols)+
  scale_x_discrete(limits=c("Fire-only","Rust'95", "Beetles",  "Beetles*fire",
                            "Rust*fire","Rust*Beetles", "Rust*fire*beetles"),
                   labels = c("Fire-only","Rust'95", "Beetles",  "Beetles*fire",
                              "Rust*fire","Rust*Beetles", "Rust*fire*beetles"))+
  #ggtitle(label="High minimum temperatures")+
  
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 5)+
  ylab("Transformed logistic regression estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=15), plot.title = element_text(hjust = .5, size=18))

oddslowgraphdirect<-ggplot(subset(oddslow3, subset=type=="direct"), aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  scale_fill_manual(values = cols)+
  scale_x_discrete(limits=c("Fire-only","Rust'95", "Beetles"),
                   labels = c("Fire-only","Rust'95", "Beetles"))+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 5)+
  ylab("Transformed estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18), plot.title = element_text(hjust = .5, size=18))

oddslowgraphdirect

oddslowgraphind<-ggplot(subset(oddslow3, subset=type2==2), aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  scale_fill_manual(values = colsind)+
  scale_x_discrete(limits=c( "Beetles*fire","Rust*fire","Rust*Beetles", "Rust*fire*beetles"),
                   labels = c("Beetles*fire", "Rust*fire","Rust*Beetles", "Rust*fire*beetles"))+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 5)+
  ylab("Transformed estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18), plot.title = element_text(hjust = .5, size=18))

oddslowgraphind

##medium mortality==================================================================================
##==================================================================================================
medmort=med%>% 
  dplyr::select(fire_only ,elevation,dbh, era95_inc_tot , inc_tot, strm_dist, lake_dist , ppt_post95,droughtvpd, droughtdew,
                beetles , south , slope ,baplot95 ,tmin_post95,plot, mort, tree_id, year, species, status )%>%
  mutate_at(scale, .vars = vars(-plot, -mort, -tree_id, -year, -species, -status,-fire_only,-beetles, -era95_inc_tot, -inc_tot))%>%
  as.data.frame(.)

medmortmodall=glmer(mort ~ fire_only + era95_inc_tot + beetles + tmin_post95 +
                      tmin_post95*beetles + tmin_post95*era95_inc_tot + tmin_post95*fire_only+
                      fire_only*beetles + era95_inc_tot*fire_only + beetles*era95_inc_tot+ 
                      (1|plot), family=binomial(link=logit), data=medmort,
                    control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                           optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(medmortmodall)

##interactions only
medmortint=glmer(mort ~ fire_only + era95_inc_tot + beetles + 
                   fire_only*beetles + era95_inc_tot*fire_only + beetles*era95_inc_tot+
                   fire_only*beetles*era95_inc_tot+(1|plot), 
                 control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                        optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                 family=binomial, data=medmort)

summary(medmortint)

oddsmed=odds_to_rr(medmortint) 

oddsmed2=as.data.frame(oddsmed)
oddsmed2=oddsmed2%>%
  mutate(upper.ci=replace_na(upper.ci, 0))

names(oddsmed2)[1:3]=c( "estimate", "lower", "upper")
oddsmed2$var=c("intercept", "Fire-only","Rust'95", "Beetles",  "Beetles*fire",
               "Rust*fire","Rust*Beetles", "Rust*fire*beetles")

oddsmed3=oddsmed2%>%
  mutate(pvalue=c(.001, .001,.001, .001, .08,.004,.02, 1))%>%
  mutate(upper2=upper-estimate, lower2=estimate-lower)%>%
  filter(var!="intercept")%>%
  mutate(type=c("direct", "direct", "direct", "ind", "ind", "ind","ind"))


oddsmedgraph<-ggplot(oddsmed3, aes(x=var, y=estimate,  fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-4, 
            color="black", position=position_dodge(0.9), size=5)+
  #ggtitle(label="Medium minimum temperatures")+
  scale_fill_manual(values=cols)+
  scale_x_discrete(limits=c("Fire-only","Rust'95", "Beetles",  "Beetles*fire",
                            "Rust*fire","Rust*Beetles", "Rust*fire*beetles"),
                   labels = c("Fire-only","Rust'95", "Beetles",  "Beetles*fire",
                              "Rust*fire","Rust*Beetles", "Rust*fire*beetles")) +
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 5)+
  ylab("Transformed logistic regression estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=13),
        axis.title=element_text(size=15), plot.title = element_text(hjust = .5, size=18))

oddsmedgraph

oddsmedgraphdirect<-ggplot(subset(oddsmed3, subset=type=="direct"), aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  scale_fill_manual(values = cols)+
  scale_x_discrete(limits=c("Fire-only","Rust'95", "Beetles"),
                   labels = c("Fire-only","Rust'95", "Beetles"))+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 5)+
  ylab("Transformed estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18), plot.title = element_text(hjust = .5, size=18))

oddsmedgraphdirect

oddsmedgraphind<-ggplot(subset(oddsmed3, subset=type=="ind"), aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-3, 
            color="black", position=position_dodge(0.9), size=5)+
  scale_fill_manual(values = colsind)+
  scale_x_discrete(limits=c( "Beetles*fire","Rust*fire","Rust*Beetles", "Rust*fire*beetles"),
                   labels = c("Beetles*fire", "Rust*fire","Rust*Beetles", "Rust*fire*beetles"))+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 5)+
  ylab("Transformed estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18), plot.title = element_text(hjust = .5, size=18))

oddsmedgraphind





##==================================================##high mortality==========================================
##============================================================================================================
highmort=high%>% 
  dplyr::select(fire_only ,elevation,dbh, era95_inc_tot , inc_tot, strm_dist, lake_dist , ppt_post95,droughtvpd, droughtdew,
                beetles , south , slope ,baplot95 ,tmin_post95,plot, mort, tree_id, year, species, status )%>%
  mutate_at(scale, .vars = vars(-plot, -mort, -tree_id, -year, -species, -status,-fire_only,-beetles, -era95_inc_tot, -inc_tot))%>%
  as.data.frame(.)


highmortmod=glmer(mort ~ fire_only + era95_inc_tot + beetles + tmin_post95 + baplot95+
                    tmin_post95*beetles + tmin_post95*era95_inc_tot + tmin_post95*fire_only+
                    fire_only*beetles + era95_inc_tot*fire_only + beetles*era95_inc_tot+ 
                    (1|plot), family=binomial(link=logit), data=highmort,
                  control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                         optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

summary(highmortmod)

highmortint=glmer(mort ~ fire_only + era95_inc_tot + beetles + 
                    fire_only*beetles + era95_inc_tot*fire_only + beetles*era95_inc_tot+
                    fire_only*beetles*era95_inc_tot+(1|plot), 
                  control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                         optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                  family=binomial, data=highmort)

highmortint=glmer(mort ~ era95_inc_tot + beetles + 
                    (1|plot), 
                  control = glmerControl(optimizer = "optimx", calc.derivs = FALSE,
                                         optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)),
                  family=binomial, data=highmort)
summary(highmortint)

oddshigh=odds_to_rr(highmortint) 

oddshigh2=as.data.frame(oddshigh)
names(oddshigh2)[1:3]=c( "estimate", "lower", "upper")
oddshigh2$var=c("intercept", "Rust'95", "Beetles")

oddshigh3=oddshigh2%>%
  mutate(pvalue=c(.9, .9, .001))%>%
  mutate(upper=replace_na(upper, 0))%>%
  mutate(upper2=upper-estimate, lower2=estimate-lower)%>%
  filter(var!="intercept")

"Rust'95"="palevioletred", "Beetles"="#56B4E9"
oddshighgraph<-ggplot(oddshigh3, aes(x=var, y=estimate, fill=var))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.6)+
  scale_fill_manual(values = c("#56B4E9","palevioletred"), labels = c("Beetles", "Rust'95")) +
  geom_errorbar(aes(ymin=estimate+upper2, ymax=estimate-lower2), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=ifelse(pvalue<.002, "**",
                             ifelse(pvalue<.05, "*",
                                    ifelse(pvalue<.099,"+", "")))), vjust=-4, 
            color="black", position=position_dodge(0.9), size=5)+
  #ggtitle(label="Low minimum temperatures")+
  theme_tufte(ticks=F)+
  guides(fill=FALSE)+
  ylim(0, 40)+
  ylab("Transformed estimates")+
  theme(axis.title.x = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=18), plot.title = element_text(hjust = .5, size=18))

oddshighgraph
oddsmedgraph
