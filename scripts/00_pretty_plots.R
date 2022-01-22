suppressMessages(library(scales))


## color by sample site
myColors_sampsite <- c("#1a9850","#a6d96a","cornflowerblue")
color_possibilities_sampsite<-c( "AC","SC","TI")
names(myColors_sampsite) <- color_possibilities_sampsite
fillscale_sampsite <- scale_fill_manual(name="Sample Site",
                                         values = myColors_sampsite, drop = T)
colscale_sampsite <- scale_color_manual(name="Sample Site",
                                        values = myColors_sampsite, drop = T)


## color by passage (values for all 3 datasets used)
pass_col<-c("grey30","#9E0142", "#D53E4F", "#F46D43", "#FDAE61","#FEC776","#FEE08B", "#E6F598","#ABDDA4","#94D4A4","#7DCBA5","#66C2A5","#4CA5B1","#3288BD", "#5E4FA2","#762a83","#3f007d")
names(pass_col)<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,21,23)


# Theme setting for all plots 
th <-   theme(axis.text=element_text(size=10),
              axis.title=element_text(size=12),
              strip.text = element_text(size = 12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))

th_present <- theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    strip.text = element_text(size = 12),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14))
