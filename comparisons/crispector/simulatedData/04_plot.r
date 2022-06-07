library(ggplot2)
pdf('04_plot.r.pdf',width=12,height=8)

d = read.table('02_summarize_sims.py.results.txt',head=T,check.names=F,sep="\t")
p1 = ggplot(d, aes(fill=level, y=count, x=idx)) + 
           geom_bar(position="stack", stat="identity") +
	   geom_vline(xintercept=53, linetype='dotted') +
	   geom_segment(aes(x=50,xend=50+23,y=-1,yend=-1)) +
	   guides(fill=guide_legend(title="CRISPECTOR classification",title.position="top", title.hjust = 0.5)) +
	   xlab('Translocation position') +
	   ylab('Read classification (count)') + 
	   ggtitle('ChromBridGE analysis of simulated reads') +
	   theme(legend.position="bottom")


d = read.table('03_run_chrombridge.py.results.txt',head=T,check.names=F,sep="\t")
#colnames(d) = c('cut_pos',1:(ncol(d)-1))

p2 = ggplot(d, aes(fill=level, y=count, x=idx)) + 
           geom_bar(position="stack", stat="identity") +
	   geom_vline(xintercept=53, linetype='dotted') +
	   geom_segment(aes(x=50,xend=50+23,y=-1,yend=-1)) +
	   guides(fill=guide_legend(title="ChromBridGE classification",title.position="top", title.hjust = 0.5)) +
	   xlab('Translocation position') +
	   ylab('Read classification (count)') +
	   ggtitle('ChromBridGE analysis of simulated reads') +
	   theme(legend.position="bottom")

print(p1)
print(p2)

library(cowplot)
plot_grid(p1, p2, labels = c('A', 'B'))

dev.off()


