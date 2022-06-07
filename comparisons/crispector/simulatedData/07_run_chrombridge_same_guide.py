from ChromBridGE import Simulation
import argparse
import random
import os
import re
from collections import defaultdict
from ChromBridGE import ChromBridGE_aln


summ_file = os.path.basename(__file__)+'.summary.txt'
max_max_pos = 0
all_res_seen_pos_counts = defaultdict(int)
all_res_tx_pos_counts = defaultdict(int)
with open(summ_file,'w') as summ_out:
    this_run_id = 0
    while(True):
        print('running ' + str(this_run_id))
        summ_out.write('running ' + str(this_run_id))
        tx_reads_file = 'sims_same_guide/05_simulate_reads_same_guide.py.sim_%d.tx.fq'%this_run_id
        sim_info_file = 'sims_same_guide/05_simulate_reads_same_guide.py.sim_%d.info'%this_run_id
        if not os.path.exists(tx_reads_file):
            break

        with open(sim_info_file,'r') as fin:
            cut_a_line = fin.readline().strip().split("\t")
            cut_a_pos = int(cut_a_line[1])
            seq_a_line = fin.readline().strip().split("\t")
            seq_a = seq_a_line[1]

            cut_b_line = fin.readline().strip().split("\t")
            cut_b_pos = int(cut_b_line[1])
            seq_b_line = fin.readline().strip().split("\t")
            seq_b = seq_b_line[1]

        read_cache = {}
        seen_pos_counts = defaultdict(int)
        tx_pos_counts = defaultdict(int)
        seqs_read = 0
        matches_read = 0
        with open(tx_reads_file,'r') as fin:
            id_line = fin.readline()
            seq_line = fin.readline().strip()
            plus_line = fin.readline()
            qual_line = fin.readline()
            while id_line:
                seqs_read += 1
                match = re.search('@cut_(\d+)_rep_\d',id_line)
                if match:
                    cut_pos = int(match.group(1))
                    if cut_pos > max_max_pos:
                        max_max_pos = cut_pos
                    matches_read += 1
                    if seq_line not in read_cache:

                        read_aln,refA_aln,refB_aln,breakpoints,aln_score,read_path = ChromBridGE_aln.nw_breakpoint(seq_line,seq_a,seq_b,
                            ref1_cut_pos=cut_a_pos,
                            ref2_cut_pos=cut_b_pos)
                        read_cache[seq_line] = breakpoints
                    else:
                        breakpoints = read_cache[seq_line]
#                    if cut_a_pos in breakpoints[0]:
#                        print(f"{read_aln=}")
#                        print(f"{refA_aln=}")
#                        print(f"{refB_aln=}")
#                        print(f"{breakpoints=}")
#
#                        asdf()

                    if cut_pos not in seen_pos_counts:
                        print('cutpos: ' + str(cut_pos) + ' txpos: ' + str(breakpoints[0]))
                        summ_out.write('\tcutpos: ' + str(cut_pos) + ' txpos: ' + str(breakpoints[0]))
                        seen_pos_counts[cut_pos] += 1
                        if cut_a_pos in breakpoints[0]:
                            tx_pos_counts[cut_pos] += 1

                id_line = fin.readline()
                seq_line = fin.readline().strip()
                plus_line = fin.readline()
                qual_line = fin.readline()

        for idx in range(max_max_pos):
            if seen_pos_counts[idx] > 0:
                all_res_seen_pos_counts[idx] += 1
            if tx_pos_counts[idx] > 0:
                all_res_tx_pos_counts[idx] += 1

        this_run_id += 1
        print('Finished ' + str(this_run_id))

all_res_file = os.path.basename(__file__)+'.results.txt'
with open(all_res_file,'w') as res_out:
    res_out.write('idx\tlevel\tcount\n')
    for idx in range(max_max_pos):
        tx_pos_count = all_res_tx_pos_counts[idx]
        seen_pos_count = all_res_seen_pos_counts[idx]
        unlikely_tx_pos_count = seen_pos_count - tx_pos_count
        res_out.write("%d\t%s\t%d\n"%(idx,'translocation',tx_pos_count))
        res_out.write("%d\t%s\t%d\n"%(idx,'unlikely translocation',unlikely_tx_pos_count))

print('Finished')

