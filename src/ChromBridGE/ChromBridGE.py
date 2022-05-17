import argparse
import os
import re
import gzip
from ChromBridGE import ChromBridGE_aln



def main():

    parser = argparse.ArgumentParser(description='ChromBridGE: Translocation detection in genome-edited reads.')
    parser.add_argument('-f','--fastq', help='Input fastq file', required=True)
    parser.add_argument('-a','--sequence_a', help='Input sequence a', required=True)
    parser.add_argument('-b','--sequence_b', help='Input sequence b', required=True)
    parser.add_argument('--seqA_cut_pos', type=int, help='Index in sequence a of predicted cut site',default=None)
    parser.add_argument('--seqB_cut_pos', type=int, help='Index in sequence b of predicted cut site',default=None)
    parser.add_argument('--match_score', help='Match score for alignment',default=3)
    parser.add_argument('--mismatch_score', help='Mismatch score for alignment',default=-1)
    parser.add_argument('--gap_score', help='Gap score for alignment',default=-2)
    parser.add_argument('--jump_score', help='Jump score for alignment',default=-3)
    parser.add_argument('--cut_pos_incentive_score', help='Incentive for jumping at a predicted cut site',default=1)
    parser.add_argument('-o','--output_file', help='Output file to write results',default=None)
    args = parser.parse_args()

    if not os.path.isfile(args.fastq):
        raise Exception('File ' + args.fastq + ' does not exist')

    output_file = args.output_file
    if output_file is None:
        root = args.fastq
        root = re.sub(".fq.gz$","",root)
        root = re.sub(".fq$","",root)
        root = re.sub(".fastq.gz$","",root)
        root = re.sub(".fastq$","",root)
        output_file = root+".ChromBridGE.fa"

    if args.fastq.endswith('.gz'):
        f_in = io.BufferedReader(gzip.open(args.fastq,'rt'))
    else:
        f_in = open(args.fastq,'rt')

    if output_file.endswith('.gz'):
        f_out = gzip.open(output_file, 'wt')
    else:
        f_out = open(output_file, 'wt')

    f_out.write("read_id\tbreakpoints\tbreakpoint_count\tbreakpoint_cumulative_distance_from_cut\tread_aln\trefA_aln\trefB_aln\n")

    total_read_count = 0
    while (1):
        id_line   = f_in.readline().strip()
        seq_line  = f_in.readline().strip()
        plus_line = f_in.readline()
        qual_line = f_in.readline().strip()

        if not qual_line : break

        total_read_count+=1
        print('total read count: ' + str(total_read_count))
        if total_read_count % 1000 == 0:
            print('total read count: ' + str(total_read_count))

        read_aln,refA_aln,refB_aln,breakpoints = ChromBridGE_aln.nw_breakpoint(seq_line,args.sequence_a,args.sequence_b,
                match_score=args.match_score,
                mismatch_score=args.mismatch_score,
                gap_score=args.gap_score,
                jump_score=args.jump_score,
                cut_pos_incentive_score=args.cut_pos_incentive_score,
                ref1_cut_pos=args.seqA_cut_pos,
                ref2_cut_pos=args.seqB_cut_pos)

        breakpoints_count = len(breakpoints[0])
        breakpoints_out = ",".join([str(x) for x in breakpoints[0]])
        breakpoint_cumulative_distance = "NA"
        if args.seqA_cut_pos and args.seqB_cut_pos:
            breakpoint_cumulative_distance = 0
            for bp1 in breakpoints[1]:
                breakpoint_cumulative_distance += abs(bp1-args.seqA_cut_pos)
            for bp2 in breakpoints[2]:
                breakpoint_cumulative_distance += abs(bp2-args.seqB_cut_pos)
            breakpoint_cumulative_distance = str(breakpoint_cumulative_distance)
        f_out.write(id_line+"\t"+breakpoints_out+"\t"+str(breakpoints_count)+"\t"+breakpoint_cumulative_distance+"\t"+read_aln+"\t"+refA_aln+"\t"+refB_aln+"\n")


if __name__ == "__main__":
    main()
