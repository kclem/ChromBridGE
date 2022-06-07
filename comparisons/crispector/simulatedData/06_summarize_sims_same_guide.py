import os
import re
import gzip
from collections import defaultdict

os.makedirs('results_same_guide', exist_ok=True)

summ_file = os.path.basename(__file__)+'.summary.txt'
all_res_file = os.path.basename(__file__)+'.results.txt'
max_max_pos = 0
all_res_counts = defaultdict(lambda:defaultdict(int))
all_res_status_levels = defaultdict(int) #keep track of unique seen levels

outroot = "sims_same_guide/05_simulate_reads_same_guide.py"
with open(summ_file,'w') as summ_out:
    this_run_id = 0
    while(True):
        tx_reads_file = outroot + '.sim_%d.tx.fq'%this_run_id
        if not os.path.exists(tx_reads_file):
            break
        crispector_tx_file = 'CRISPECTOR_OUT/'+outroot+'.sim_%d/crispector_output/tx_reads_with_primer_inconsistency.csv'%this_run_id
        crispector_ref1_file = 'CRISPECTOR_OUT/'+outroot+'.sim_%d/crispector_output/ref1/treatment_aligned_reads.csv.gz'%this_run_id
        crispector_ref2_file = 'CRISPECTOR_OUT/'+outroot+'.sim_%d/crispector_output/ref2/treatment_aligned_reads.csv.gz'%this_run_id
        crispector_treatment_unmatched_file = 'CRISPECTOR_OUT/'+outroot+'.sim_%d/crispector_output/treatment_unmatched_reads.fa.gz'%this_run_id

        print('Reading simulated reads')
        sim_reads_dict = {}
        sim_cut_idx_dict = {}
        matches_read = 0
        seqs_read = 0
        max_pos = 0
        sim_count_dict = defaultdict(int)
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
                    # some sequences will be the same for different cuts - e.g. if the bases match
                    # e.g. mappings of seq_line > cut_pos will not have unique cut_pos values
#                    if seq_line in sim_reads_dict and sim_reads_dict[seq_line] != cut_pos:
#                        print('mismatch for %s:%d vs %d\n'%(seq_line,sim_reads_dict[seq_line],cut_pos))

                    sim_reads_dict[seq_line] = cut_pos
                    # but all mappings of cut_pos>seq_line should be the same
                    if cut_pos in sim_cut_idx_dict and sim_cut_idx_dict[cut_pos] != seq_line:
                        print('mismatch for %s:\n%s\n%s\n'%(cut_pos,sim_cut_idx_dict[cut_pos],seq_line))
                    sim_cut_idx_dict[cut_pos] = seq_line
                    matches_read += 1
                    sim_count_dict[seq_line] += 1
                    if cut_pos > max_pos:
                        max_pos = cut_pos
                    if cut_pos > max_max_pos:
                        max_max_pos = cut_pos

                id_line = fin.readline()
                seq_line = fin.readline().strip()
                plus_line = fin.readline()
                qual_line = fin.readline()

        print('got ' + str(matches_read) + ' matches (' + str(len(sim_reads_dict)) + ' unique) from ' + str(seqs_read) + ' simulated tx sequences')

        print('Reading unaligned reads')
        crispector_status_dict = {}
        unaligned_count = 0
        unaligned_line_count = 0
        crispector_count_dict = defaultdict(int)
        with gzip.open(crispector_treatment_unmatched_file,'rt') as fin:
            info_line = fin.readline()
            seq_line = fin.readline().strip()
            while info_line:
                unaligned_line_count += 1
                match = re.search('unaligned read with (\d+) copies in the original',info_line)
                if match:
                    read_count = int(match.group(1))
                    unaligned_count += read_count
                    crispector_count_dict[seq_line] += read_count
                crispector_status_dict[seq_line] = 'Unaligned'
                info_line = fin.readline()
                seq_line = fin.readline().strip()

        print('Read ' + str(unaligned_count) + ' unaligned reads (' + str(unaligned_line_count) + ' lines)')

        print('Reading tx reads')
        tx_line_count = 0
        tx_read_count = 0
        with open(crispector_tx_file,'r') as fin:
            head_line = fin.readline()
            for line in fin:
                tx_line_count += 1
                line_els = line.strip().split(",")
                read_count = int(line_els[1])
                tx_read_count += read_count
                read = line_els[5]
                tx_status = line_els[13] == 'True'
                crispector_count_dict[read] += read_count

                if read in crispector_status_dict:
                    print('error: read is already in dict')
                if tx_status:
                    crispector_status_dict[read] = 'Translocated'
                else:
                    crispector_status_dict[read] = 'Not translocated'
        print('Read ' + str(tx_read_count) + ' tx reads (' + str(tx_line_count) + ' lines)')

        print('Reading aln ref1 reads')
        r1_line_count = 0
        r1_read_count = 0
        with gzip.open(crispector_ref1_file,'rt') as fin:
            head_line = fin.readline()
            for line in fin:
                r1_line_count += 1
                line_els = line.strip().split(",")
                read_count = int(line_els[1])
                r1_read_count += read_count
                read = line_els[0]
                edited_status = line_els[18] == 'True'
                crispector_count_dict[read] += read_count

                if read in crispector_status_dict:
                    print('error: read is already in dict')
                if edited_status:
                    crispector_status_dict[read] = 'Ref1_edited'
                else:
                    crispector_status_dict[read] = 'Ref1_unedited'
        print('Read ' + str(r1_read_count) + ' ref1 reads (' + str(r1_line_count) + ' lines)')

        print('Reading aln ref2 reads')
        r2_line_count = 0
        r2_read_count = 0
        with gzip.open(crispector_ref2_file,'rt') as fin:
            head_line = fin.readline()
            for line in fin:
                r2_line_count += 1
                line_els = line.strip().split(",")
                read_count = int(line_els[1])
                r2_read_count += r2_read_count
                read = line_els[0]
                edited_status = line_els[18] == 'True'
                crispector_count_dict[read] += read_count

                if read in crispector_status_dict:
                    print('error: read is already in dict')
                if edited_status:
                    crispector_status_dict[read] = 'Ref2_edited'
                else:
                    crispector_status_dict[read] = 'Ref2_unedited'
        print('Read ' + str(r2_read_count) + ' ref2 reads (' + str(r2_line_count) + ' lines)')


        res_array = ['Not found']*(max_pos+1)
        res_file = os.path.join('results',str(this_run_id)+'.txt')
        with open(res_file,'w') as fout:
            for idx in range(max_pos):
                read_seq = sim_cut_idx_dict[idx]
                crispector_status = 'Unknown'
                if read_seq in crispector_status_dict:
                    crispector_status = crispector_status_dict[read_seq]
                fout.write("%s\t%s\t%s\t%s\t%s\n"%(read_seq,idx,crispector_status,sim_count_dict[read_seq],crispector_count_dict[read_seq]))
                res_array[idx] = crispector_status
                all_res_counts[idx][crispector_status] += 1
                all_res_status_levels[crispector_status] += 1
        print('Wrote ' + res_file)

        summ_out.write(str(this_run_id)+"\t"+"\t".join(res_array)+"\n")

        this_run_id += 1

levels = sorted(all_res_status_levels.keys())
with open(all_res_file,'w') as res_out:
    res_out.write('idx\tlevel\tcount\n')
    for idx in range(max_max_pos):
        for level in levels:
            count = all_res_counts[idx][level]
            res_out.write("%d\t%s\t%d\n"%(idx,level,count))

print('Finished')
