library(ggplot2)
library(MASS)
library(forecast)


dataFrame <- read.csv("../data/traces_reports_csv/control_time_full.csv")
dataFrame$executionTimeMS = dataFrame$executionTime*1e-6
ggplot(dataFrame, aes(x=executionTimeMS)) +
  geom_histogram(binwidth = 0.002, aes(y=after_stat(density)), fill="black") +
  xlab("Execution time (ms)") +
  ylab("Density") + theme_bw() +  theme(text = element_text(size = 14)) + 
  theme(axis.text.x= element_text(size = 14)) + 
  theme(axis.text.y= element_text(size = 14)) + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  scale_x_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5), limits=c(0, 0.55), expand=c(0,0)) +
  scale_y_continuous(breaks=c(25, 50, 75), limits=c(0,100),expand=c(0,0))
filename = "trace_density.png"
ggsave(filename, width=10, height=5)
filename = "trace_density.eps"
ggsave(filename, width=10, height=5)

dataFrame$index = 1:nrow(dataFrame)
ggplot(dataFrame, aes(x=index, y=executionTimeMS)) +
  geom_point(size=0.3, shape=3) +
  #ylim(0.125, 0.525) + 
  scale_x_continuous(limits=c(0, 48000), expand=c(0,0)) +
  scale_y_continuous(limits=c(0, 0.55), expand=c(0,0)) +
  ylab("Execution time (ms)") +
  xlab("Job index" )+ 
  theme_bw() +  theme(text = element_text(size = 14)) +
  theme(axis.text.x= element_text(size = 14)) + 
  theme(axis.text.y= element_text(size = 14)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
filename = "trace_sequence.png"
ggsave(filename, height = 5, width = 10)
filename = "trace_sequence.eps"
ggsave(filename, height = 5, width = 10)


dataFrame <- read.csv("../data/traces_reports_csv/control_time_full.csv")

ggAcf(
  dataFrame$executionTime,
  lag.max = 20,
  plot = TRUE,
) + scale_x_continuous(breaks=c(5, 10, 15), limits=c(0, 20), expand=c(0,0)) +
  scale_y_continuous(breaks=c(0.1, 0.2, 0.3, 0.4, 0.5), limits=c(0, 0.6), expand=c(0,0)) +
  theme_bw() +  theme(text = element_text(size = 14)) +
  theme(axis.text.x= element_text(size = 14)) + 
  theme(axis.text.y= element_text(size = 14)) + 
  theme(plot.title=element_blank()) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

filename = "ggplot_acorr.png"
ggsave(filename, width=10, height=5)
filename = "ggplot_acorr.eps"
ggsave(filename, width=10, height=5)



      
