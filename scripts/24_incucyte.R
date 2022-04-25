library(here)


library(reshape)
library(ggplot2)

source(here("scripts","00_pretty_plots.R"))


library(pzfx)
library(plyr)

pzfx_tables("data/incucyte/ARF_Passage_paper_COMBINED.pzfx")

incucyte_df_SC <- read_pzfx("data/incucyte/ARF_Passage_paper_COMBINED.pzfx","T528 & T527 SC P1-P8")
incucyte_df_TI <- read_pzfx("data/incucyte/ARF_Passage_paper_COMBINED.pzfx","T528 & T527 TI P1-P8")

incucyte_df_SC<-melt(incucyte_df_SC, id="Time (Days)")
incucyte_df_SC$variable<-as.character(incucyte_df_SC$variable)
incucyte_df_SC$passage<-sapply(1:nrow(incucyte_df_SC), function(x) strsplit(incucyte_df_SC$variable[x]," P|_")[[1]][2])
incucyte_df_SC$sample<-sapply(1:nrow(incucyte_df_SC), function(x) strsplit(incucyte_df_SC$variable[x]," P|_")[[1]][3])
colnames(incucyte_df_SC)[1]<-"time"


incucyte_df_TI<-melt(incucyte_df_TI, id="Time (Days)")
incucyte_df_TI$variable<-as.character(incucyte_df_TI$variable)
incucyte_df_TI$passage<-sapply(1:nrow(incucyte_df_TI), function(x) strsplit(incucyte_df_TI$variable[x]," P|_")[[1]][2])
incucyte_df_TI$sample<-sapply(1:nrow(incucyte_df_TI), function(x) strsplit(incucyte_df_TI$variable[x]," P|_")[[1]][3])
colnames(incucyte_df_TI)[1]<-"time"


ggplot(incucyte_df_SC, aes(time,value, group=passage ))+geom_point()+
  stat_summary(fun.y=mean, colour="red", geom="line", aes(group = passage))

std <- function(x) sd(x, na.rm=T)/sqrt(length(x))

                  # agg = ddply(incucyte_df_SC, .(time, passage), function(x){
                  #   c(mean=mean(x$value, na.rm=T), sterr = std(x$value))
                  # })
                  # 
                  # ggplot(agg, aes(x=time, y=mean, group=passage)) + 
                  #   geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
                  #   geom_line(aes(colour=passage)) +
                  #   geom_point(aes(colour=passage)) +
                  #   labs(x = "Time (Days)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+ylim(0,600)+th_present+theme_bw()+
                  #   scale_color_manual(values=pass_col,name="Passage\nNumber")
                  # 
                  # 
                  # agg = ddply(incucyte_df_TI, .(time, passage), function(x){
                  #   c(mean=mean(x$value, na.rm=T), sterr = std(x$value))
                  # })
                  # 
                  # ggplot(agg, aes(x=time, y=mean, group=passage)) + 
                  #   geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
                  #   geom_line(aes(colour=passage)) +
                  #   geom_point(aes(colour=passage)) +
                  #   labs(x = "Time (Days)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+ylim(0,600)+th_present+theme_bw()+
                  #   scale_color_manual(values=pass_col,name="Passage\nNumber")


incucyte_df_TI$segment<-"TI"
incucyte_df_SC$segment<-"SC"

incucyte_df<-rbind(incucyte_df_SC, incucyte_df_TI)

agg = ddply(incucyte_df, .(time, passage, segment), function(x){
  c(mean=mean(x$value, na.rm=T), sterr = std(x$value))
})

ggplot(agg, aes(x=time, y=mean, group=passage)) + 
  geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
  geom_line(aes(colour=passage)) +
  geom_point(aes(colour=passage)) +
  labs(x = "Time (Days)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+ylim(0,500)+th_present+theme_bw()+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+facet_wrap(~segment)

ggsave("figs/incucyte_growth_allpassage.pdf", width = 8, height=3)
ggsave("figs/jpeg/incucyte_growth_allpassage.jpeg", width = 8, height=3)



incucyte_df$hilo<-sapply(1:nrow(incucyte_df), function(x){
  if(incucyte_df$passage[x]<=4){"Low\n(1-4)"}else{"High\n(5-8)"}
})

incucyte_df$hilo<-factor(incucyte_df$hilo, levels=c("Low\n(1-4)","High\n(5-8)"))

agg = ddply(incucyte_df, .(time, hilo, segment), function(x){
  c(mean=mean(x$value, na.rm=T), sterr = std(x$value))
})

ggplot(agg, aes(x=time, y=mean, group=hilo)) + 
  geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
  geom_line(aes(colour=hilo)) +
  geom_point(aes(colour=hilo), size=0.75) +
  labs(x = "Time (Days)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+th_present+theme_bw()+
  scale_color_manual(values=c("#0055ff","#ff2b00"),name="Passage")+facet_wrap(~segment)

ggsave("figs/incucyte_growth_hilo.pdf", width = 5, height=2)
ggsave("figs/jpeg/incucyte_growth_hilo.jpeg", width = 6, height=2)


aov(agg$mean ~ agg$hilo + agg$segment)
summary(aov(agg$mean ~ agg$hilo + agg$segment))



############################### 
## forskolin
############################### 
pzfx_tables("data/incucyte/Forskolin assay for the passage paper 25.03.22.pzfx")

incucyte_df_early <- read_pzfx("data/incucyte/Forskolin assay for the passage paper 25.03.22.pzfx","early")
incucyte_df_late <- read_pzfx("data/incucyte/Forskolin assay for the passage paper 25.03.22.pzfx","late")


incucyte_df_early<-melt(incucyte_df_early, id="ROWTITLE")
incucyte_df_early$variable<-as.character(incucyte_df_early$variable)
incucyte_df_early$treatment<-sapply(1:nrow(incucyte_df_early), function(x) strsplit(incucyte_df_early$variable[x],"_")[[1]][1])
incucyte_df_early$passage<-"early"

incucyte_df_late<-melt(incucyte_df_late, id="ROWTITLE")
incucyte_df_late$variable<-as.character(incucyte_df_late$variable)
incucyte_df_late$treatment<-sapply(1:nrow(incucyte_df_late), function(x) strsplit(incucyte_df_late$variable[x],"_")[[1]][1])
incucyte_df_late$passage<-"late"

incucyte_df<-rbind(incucyte_df_late, incucyte_df_early)


ggplot(incucyte_df, aes(ROWTITLE,value, group=passage ))+geom_point()+
  stat_summary(fun.y=mean, colour="red", geom="line", aes(group = passage))

std <- function(x) sd(x, na.rm=T)/sqrt(length(x))


agg = ddply(incucyte_df, .(ROWTITLE, passage, treatment), function(x){
  c(mean=mean(x$value, na.rm=T), sterr = std(x$value))
})


ggplot(agg, aes(x=ROWTITLE, y=mean, group=passage)) +
  geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
  geom_line(aes(colour=passage)) +
  geom_point(aes(colour=passage),size=0.75) +
  labs(x = "Time (Hours)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+ylim(0,500)+th_present+theme_bw()+
  scale_color_manual(values=c("#0055ff","#ff2b00"),name="Passage")+facet_wrap(~treatment)

ggsave("figs/incucyte_forskolin_hilo.pdf", width = 5, height=3)
ggsave("figs/jpeg/incucyte_forskolin_hilo.jpeg", width = 5, height=3)


agg_main<-agg[which(agg$treatment%in%c("IFNg+Forskolin","TNFa+Forskolin")),]
agg_supp<-agg[which(!(agg$treatment%in%c("IFNg+Forskolin","TNFa+Forskolin"))),]

ggplot(agg_main, aes(x=ROWTITLE, y=mean, group=passage)) +
  geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
  geom_line(aes(colour=passage)) +
  geom_point(aes(colour=passage), size=0.75) +
  labs(x = "Time (Hours)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+th_present+theme_bw()+
  scale_color_manual(values=c("#0055ff","#ff2b00"),name="Passage")+facet_wrap(~treatment)

ggsave("figs/incucyte_forskolin_hilo_main.pdf", width = 6, height=2)
ggsave("figs/jpeg/incucyte_forskolin_hilo_main.jpeg", width = 6, height=2)


ggplot(agg_supp, aes(x=ROWTITLE, y=mean, group=passage)) +
  geom_errorbar(aes(ymin=mean+sterr, ymax=mean-sterr), color="lightgrey", width=.2) +
  geom_line(aes(colour=passage)) +
  geom_point(aes(colour=passage),size=0.75) +
  labs(x = "Time (Hours)", y = "Organoid Area\n(normalized to 0hr (%))", colour = "Passage")+th_present+theme_bw()+
  scale_color_manual(values=c("#0055ff","#ff2b00"),name="Passage")+facet_wrap(~treatment)

ggsave("figs/incucyte_forskolin_hilo_supp.pdf", width = 9, height=2)
ggsave("figs/jpeg/incucyte_forskolin_hilo_supp.jpeg", width = 9, height=2)

aov(agg$mean ~ agg$passage + agg$treatment)
summary(aov(agg$mean ~ agg$passage + agg$treatment))
