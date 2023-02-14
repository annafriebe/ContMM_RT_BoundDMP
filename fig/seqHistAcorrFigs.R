library(ggplot2)
library(MASS)
library(forecast)


dataFrame <- read.csv("../data/traces_reports_csv/control_time_full.csv")
dataFrame$executionTimeMS = dataFrame$executionTime*1e-6
ggplot(dataFrame, aes(x=executionTimeMS)) +
  geom_histogram(binwidth = 0.001, aes(y=after_stat(density)), fill="grey") +
  xlab("Execution time (ms)") +
  ylab("Density") + theme_bw() +  theme(text = element_text(size = 28)) + 
  theme(axis.text.x= element_text(size = 28)) + 
  theme(axis.text.y= element_text(size = 28)) + 
  scale_x_continuous(breaks=c(0.15, 0.2, 0.25), limits=c(0.1, 0.3))
filename = "trace_density.png"
ggsave(filename, width=6, height=6)
filename = "trace_density.eps"
ggsave(filename, width=6, height=6)

dataFrame$index = 1:nrow(dataFrame)
ggplot(dataFrame, aes(x=index, y=executionTimeMS)) +
  geom_point(size=0.1) +
  ylim(0.125, 0.525) + 
  ylab("Execution time (ms)") +
  xlab("Job index" )+ 
  theme_bw() +  theme(text = element_text(size = 14)) +
  theme(axis.text.x= element_text(size = 14)) + 
  theme(axis.text.y= element_text(size = 14)) 
filename = "trace_sequence.png"
ggsave(filename, height = 3, width = 6)
filename = "trace_sequence.eps"
ggsave(filename, height = 3, width = 6)


dataFrame <- read.csv("../data/traces_reports_csv/control_time_full.csv")

ggAcf(
  dataFrame$executionTime,
  lag.max = 10,
  plot = TRUE,
) + theme_bw() +  theme(text = element_text(size = 28)) +
  theme(axis.text.x= element_text(size = 28)) + 
  theme(axis.text.y= element_text(size = 28)) + 
  theme(plot.title=element_blank())

filename = "ggplot_acorr.png"
ggsave(filename, width=6, height=6)
filename = "ggplot_acorr.eps"
ggsave(filename, width=6, height=6)



      
