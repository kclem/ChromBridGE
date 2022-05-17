library(ggplot2)
pdf('03_plot.r.pdf',width=12,height=8)
d = read.table('02_summarize_sims.py.results.txt',head=T,check.names=F)
#colnames(d) = c('cut_pos',1:(ncol(d)-1))

ggplot(d, aes(fill=level, y=count, x=idx)) + 
           geom_bar(position="stack", stat="identity") +
	   geom_vline(xintercept=53, linetype='dotted') +
	   geom_segment(aes(x=50,xend=50+23,y=-1,yend=-1)) +
	   guides(fill=guide_legend(title="CRISPECTOR classification")) +
	   xlab('Translocation position') +
	   ylab('Read classification (count)')

dev.off()


