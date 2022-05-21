from ChromBridGE import Simulation
import argparse
import random
import os
import re
from collections import defaultdict
from ChromBridGE import ChromBridGE_aln

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))

cut_a_pos = 146
ref_a = "GATGTGAGGACTAAGATTCTGCTTTGAAAGCAGAAAGGGAAGAGTTGCCACCATTTTAACATTTAGTGACTAGCTTATTTATATCCCCCTTGGCCACAGTGACAAAAACACAGCTTTATTTCTAAAGACAACCACACGGGCCATTATTCAGCCAGCCTTCTCAACTGATGATGACGTGAACTGTGAATCCACATTTTTTCCTGTTAATAGACAACTGAAATGGAACAAATATCCTGAACTCAC"

cut_b_pos = 104
ref_b = "GATCTCATTTTGGATGTAAGCAGGGAAGAGCCTTGAAGAAAGTATTTATGCCTTGTGCTATTTGGTGGGGGGAGGGGTAGAAGGAGCAGGGAAGCCTGGCTGAATAAAGGGTAAGTGAGAGAAATGAGACAAGGCTCCTCCATCTCCTGTGGAGGGTGAGTTAGAGCACAGACTCAGGAAACTATGTCAAGAAATGAGCG"

summ_file = os.path.basename(__file__)+'.summary.txt'
details_file = os.path.basename(__file__)+'.details.txt'
plot_file = os.path.basename(__file__)+'.plot.r'
max_max_pos = 0

all_res_seen_pos_a_counts = defaultdict(int)
all_res_tx_pos_a_counts = defaultdict(int)
all_res_seen_pos_b_counts = defaultdict(int)
all_res_tx_pos_b_counts = defaultdict(int)

likely_AB_count = 0
likely_BA_count = 0
unlikely_tx_count = 0

plot_entries = []
with open(summ_file,'w') as summ_out, open(details_file,'w') as details_out:
    crispector_tx_output = "CRISPECTOR/crispector_output/tx_reads_with_primer_inconsistency.csv"

    tx_line_count = 0
    tx_read_count = 0
    with open(crispector_tx_output,'r') as fin:
        head_line = fin.readline().strip()
        head_line_els = head_line.split(",")
        summ_out.write("\t".join(head_line_els) + "\tcuts\tChromBridGE_status\n")
        for line in fin:
            tx_line_count += 1
            line_els = line.strip().split(",")
            read_count = int(line_els[1])
            tx_read_count += read_count
            read = line_els[5]
            crispector_tx_status = line_els[13] == 'True'

            read_aln,refA_aln,refB_aln,breakpoints,aln_score,read_path = ChromBridGE_aln.nw_breakpoint(read,ref_a,ref_b,
                            ref1_cut_pos=cut_a_pos,
                            ref2_cut_pos=cut_b_pos)


#            print(f"{aln_score=}")
#            print(f"{breakpoints=}")
#            print(f"{read_aln=}")
#            print(f"{refA_aln=}")
#            print(f"{refB_aln=}")
#
            tx_status = 'Not Tx'
            #discard alignments with multiple breakpoints
            if len(breakpoints[0]) == 1:
                for idx,b in enumerate(breakpoints[0]):
                    all_res_seen_pos_a_counts[breakpoints[1][idx]] += read_count
                    all_res_seen_pos_b_counts[breakpoints[2][idx]] += read_count
                    if breakpoints[1][idx] <= cut_a_pos and breakpoints[2][idx] >= cut_b_pos and read_path[idx] == 1 and read_path[idx+1] == 2:
                        likely_AB_count += 1
                        tx_status = 'Tx A>B'
                        all_res_tx_pos_a_counts[breakpoints[1][idx]] += read_count
                        all_res_tx_pos_b_counts[breakpoints[2][idx]] += read_count
                    elif breakpoints[1][idx] >= cut_a_pos and breakpoints[2][idx] <= cut_b_pos and read_path[idx] == 2 and read_path[idx+1] == 1:
                        likely_BA_count += 1
                        tx_status = 'Tx B>A'
                        all_res_tx_pos_a_counts[breakpoints[1][idx]] += read_count
                        all_res_tx_pos_b_counts[breakpoints[2][idx]] += read_count
                    else:
                        unlikely_tx_count += 1

                    if read_path[idx] == 1:
                        plot_entry = (1,breakpoints[1][idx],breakpoints[2][idx],tx_status)
                    else:
                        plot_entry = (2,breakpoints[2][idx],breakpoints[1][idx],tx_status)
                    plot_entries.append(plot_entry)
#            print('TX')
#            print(f"{ref_a=}")
#            print(f"{ref_b=}")
#            print(f"{read_aln=}")
#            print(f"{refA_aln=}")
#            print(f"{refB_aln=}")
#            print(f"{breakpoints=}")
#            print(f"{aln_score=}")
#
#            asdf()
            print('CRISPECTOR: ' + str(crispector_tx_status) + ' ' + line_els[0] + ' ChromBridGE: ' + tx_status)

            summ_out.write("\t".join(line_els) + "\t"+";".join([str(x) for x in breakpoints]) + "\t"+tx_status+"\n")
            details_out.write(line)
            details_out.write(f"{refA_aln=}\n")
            details_out.write(f"{read_aln=}\n")
            details_out.write(f"{refB_aln=}\n")
            details_out.write("ChromBridGE: "+ tx_status +"\n")
            details_out.write(f"{breakpoints=}\n")
            details_out.write(f"{read_path=}\n")
            details_out.write(f"{aln_score=}\n")
            details_out.write("===================\n\n")

all_res_a_file = os.path.basename(__file__)+'.results_a.txt'
with open(all_res_a_file,'w') as res_out:
    res_out.write('idx\tlevel\tcount\n')
    for idx in range(len(ref_a)):
        tx_pos_count = all_res_tx_pos_a_counts[idx]
        seen_pos_count = all_res_seen_pos_a_counts[idx]
        unlikely_tx_pos_count = seen_pos_count - tx_pos_count
        res_out.write("%d\t%s\t%d\n"%(idx,'translocation',tx_pos_count))
        res_out.write("%d\t%s\t%d\n"%(idx,'unlikely translocation',unlikely_tx_pos_count))

all_res_b_file = os.path.basename(__file__)+'.results_b.txt'
with open(all_res_b_file,'w') as res_out:
    res_out.write('idx\tlevel\tcount\n')
    for idx in range(len(ref_b)):
        tx_pos_count = all_res_tx_pos_b_counts[idx]
        seen_pos_count = all_res_seen_pos_b_counts[idx]
        unlikely_tx_pos_count = seen_pos_count - tx_pos_count
        res_out.write("%d\t%s\t%d\n"%(idx,'translocation',tx_pos_count))
        res_out.write("%d\t%s\t%d\n"%(idx,'unlikely translocation',unlikely_tx_pos_count))

b_bp = 20 #buffer
likely_tx_count = likely_AB_count + likely_BA_count
with open(plot_file,'w') as plot_out:
    xmax = max(len(ref_a),len(ref_b))
    ymax = len(plot_entries)*2 + 2
    plot_out.write('pdf("'+plot_file+'.pdf")\n')
    plot_out.write('layout(matrix(c(1,2),1,2,byrow=T),widths=c(1,3))\n')
    plot_out.write(f"barplot(c({likely_AB_count},{likely_BA_count},{unlikely_tx_count}),names=c('Tx A>B','Tx B>A','Unlikely Tx'),ylab='Count',xlab='Category')\n")
    plot_out.write(f"plot(0,type='n',ylim=c(0-0.5,{ymax}+0.5),xlim=c(0,{xmax}),yaxt='n',ylab='Translocations between sequences',xlab='bp')\n")
    plot_out.write(f"axis(side=2,at=c(0,{ymax}),labels=c('RAG2_5','RAG2_1'),las=2)\n")
    plot_out.write("#reference 1\n")
    plot_out.write(f"segments(0,{ymax},{len(ref_a)},{ymax},lwd=2)\n")
    plot_out.write(f"segments({cut_a_pos},{ymax}+0.5,{cut_a_pos},{ymax}-0.5,lwd=2,lty=2)\n")
    plot_out.write("#reference 2\n")
    plot_out.write(f"segments(0,0,{len(ref_b)},0,lwd=2)\n")
    plot_out.write(f"segments({cut_b_pos},0+0.5,{cut_b_pos},0-0.5,lwd=2,lty=2)\n")
    plot_out.write(f"cols = c('black','red3')\n")
    plot_out.write(f"s = seq(3)\n")
    plot_out.write(f"legend('right',bty='n',col=cols,legend=c('Unlikely Tx N={unlikely_tx_count}','Likely Tx N={likely_tx_count}'))\n")

    seq1_y = ymax-1
    seq2_y = seq1_y - len(plot_entries)
    for plot_entry in plot_entries:
        #starts on seq 1
        xs = [plot_entry[1]-b_bp,plot_entry[1],plot_entry[2],plot_entry[2]+b_bp]
        if plot_entry[0] == 1:
            ys = [seq1_y,seq1_y,seq2_y,seq2_y]
        else:
            ys = [seq2_y,seq2_y,seq1_y,seq1_y]
        plot_out.write(f"xs = c("+",".join([str(x) for x in xs])+")\n")
        plot_out.write(f"ys = c("+",".join([str(y) for y in ys])+")\n")
        this_col_ind = 1
        if plot_entry[3] != 'Not Tx':
            this_col_ind = 2
        plot_out.write(f"segments(xs[s],ys[s],xs[s+1],ys[s+1],col=cols[{this_col_ind}])\n")
        seq1_y -= 1
        seq2_y -= 1

