#Plot distribution of AS event types for genes associated with survival in MES or PN

library(ggplot2)

data = data.frame(direction=character(),event=character(),value=numeric())
data = rbind(data,data.frame(direction="MES GSCs",event="SE",value=(2/12*100)))
data = rbind(data,data.frame(direction="PN GSCs",event="SE",value=(7/12*100)))
data = rbind(data,data.frame(direction="MES GSCs",event="RI",value=(0/12*100)))
data = rbind(data,data.frame(direction="PN GSCs",event="MXE",value=(2/12*100)))
data = rbind(data,data.frame(direction="MES GSCs",event="ASS",value=(2/12*100)))

pdf("plot_AS_events_survival.pdf",height=2,width=5)
ggplot(data, aes(x = event, y = value, fill = direction)) + geom_bar(stat = 'identity', position = 'stack') + scale_fill_manual(values = c("#FF6F65E0","#00989CCC")) + theme(text = element_text(size=20),legend.title = element_blank()) + scale_y_continuous(labels = function(x) paste0(x, "%"), limits=c(0,100), breaks=seq(0,100,20), expand=c(0,0)) + xlab("") + ylab("percentage of AS events") + coord_flip() + theme_classic()
dev.off()
