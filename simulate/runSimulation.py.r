library(ggplot2)
pdf('runSimulation.py.r.pdf')
m = data.frame()
for (i in 1:21)
{
	input_file = paste0('simulations/center_guideMuts.',i,'.ChromBridGE.txt')
	dd = read.table(input_file,head=1,sep="\t")
	dd$num_mutations = i
	m = rbind(m,dd)
}

ggplot(m,aes(x=breakpoints_seqA,group=num_mutations,fill=num_mutations)) +
	geom_density(alpha=0.4)+
	ylab("Density") +
	xlab("Position of predicted cut")

