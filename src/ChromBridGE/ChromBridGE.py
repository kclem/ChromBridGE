import argparse
import os
import re
import gzip
from ChromBridGE import ChromBridGE_aln
from ChromBridGE import ChromBridGE_tx


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

        aln_info = ChromBridGE_aln.nw_breakpoint(seq_line,args.sequence_a,args.sequence_b,
                match_score=args.match_score,
                mismatch_score=args.mismatch_score,
                gap_score=args.gap_score,
                jump_score=args.jump_score,
                cut_pos_incentive_score=args.cut_pos_incentive_score,
                ref1_cut_pos=args.seqA_cut_pos,
                ref2_cut_pos=args.seqB_cut_pos)

        f_out.write(id_line+"\t"+breakpoints_out+"\t"+str(breakpoints_count)+"\t"+breakpoint_cumulative_distance+"\t"+read_aln+"\t"+refA_aln+"\t"+refB_aln+"\n")


def analyze_read(read_seq, ref1_seq, ref2_seq,
                    ref1_cut_pos=None,
                    ref2_cut_pos=None,
                    match_score=3,
                    mismatch_score=-1,
                    gap_score=-2,
                    perimeter_gap_extension_score=0,
                    jump_score=-12, # four matches
                    cut_pos_jump_incentive_score=1,
                    min_num_bases_beyond_cut=4,
                    min_num_bases_before_cut=4,
                    mismatch_tolerance=0,
                    gap_tolerance=0,
		            debug=False):
    """
    Computes the optimal alignment of a read to two sequences, locating the optimal break between the two reads, refines the site of translocation by identifying sites with at least num_bases_to_check exact matches, and determines whether the sequences are compatible with a translocation event.

    params:
        read_seq: read to align to the other two sequences
        ref1_seq: first sequence to align to
        ref2_seq: second sequence to align to

        ref1_cut_pos: position of predicted cut site in ref1
        ref2_cut_pos: position of predicted cut site in ref2

        min_num_bases_beyond_cut: Min number of matching bases that must be seen beyond the cut site for a breakpoint to be reported beyond the cut site
        min_num_bases_before_cut: Min number of matching bases that must be seen before the cut site - the position of this first match will be reported
        mismatch_tolerance: int How many mismatches to tolerate before returning false
        gap_tolerance: int How many gaps to tolerate before returning false

        match_score: score for adding a match in alignment (positive)
        mismatch_score: score for adding a mismatch in alignment
        gap_score: score for adding a gap in alignment
        perimeter_gap_extension_score: score for adding a gap in the first/last column/row, corresponding to gaps at the beginning or ends of sequences
        jump_score: score for jumping between ref1 and ref2
        cut_pos_jump_incentive_score: score incentive for jumping at a predicted cut position
        debug: print intermediate debug information

    returns:
        tuple of:
        aln_info: dict with keys:
            read_aln: aligned sequence of read
            ref1_aln: sequence of ref1 aligned to read
            ref2_aln: sequence of ref2 aligned to read
            breakpoints_read: indices in read of the breakpoints discovered
            breakpoints_ref1: indices in ref1 of breakpoints in the optimal alignment
            breakpoints_ref2: indices in ref2 of breakpoints in the optimal alignment
            aln_score: score of alignment
            read_path: index of ref that the read is aligned to, corresponding to the break points (there will be len(breakpoints)+1 items in read_path)

        tx_info: dict with keys:
            final_read_str: string of read alignment
            final_ref1_str: string of bases for which read aligns to ref1 - including '~' for trimmed insertions at translocation sites
            final_ref2_str: string of bases for which read aligns to ref2 - including '~' for trimmed insertions at translocation sites
            final_path: indices of the ref that the read is aligned to
            final_breakpoints_ref1: locations in the ref1 corresponding to translocations in ref1 (same length as final_breakpoints_read_ref1)
            final_breakpoints_ref2: locations in the ref2 corresponding to translocations in ref2 (same length as final_breakpoints_read_ref2)
            bp_match_ref1: number of bp in final alignment that match read and ref1 exactly
            bp_match_ref2: number of bp in final alignment that match read and ref2 exactly
            bp_insertion: number of bp that are inserted (match neither reference)

            is_tx: boolean for whether the read looks like a translocation
            tx_status: string with details for tx result
            left_distance: int, number of bp the read extends beyond the cut to the right (from the left-identified reference)
            right_distance: int, number of bp the read extends beyond the cut to the left (from the right-identified reference)
                Read: AAAA|TTTT
                Ref:  AAAATT  > left-distance is 2 because it extended beyond the cut by 2bp
                Read: AAAA|TTTT
                Ref:  AA      > left-distance is -2 because it was within the cut by 2bp
            tx_lucky_insertions: sum of the left- and right- distances if they extend beyond the cut. If the cut actually happened, these would be lucky insertions that happened to match the uncut reference sequence
    """
    aln_info = ChromBridGE_aln.nw_breakpoint(
                read_seq,
                ref1_seq,
                ref2_seq,
                match_score=match_score,
                mismatch_score=mismatch_score,
                gap_score=gap_score,
                perimeter_gap_extension_score=perimeter_gap_extension_score,
                jump_score=jump_score,
                cut_pos_jump_incentive_score=cut_pos_jump_incentive_score,
                ref1_cut_pos = ref1_cut_pos,
                ref2_cut_pos = ref2_cut_pos,
                debug=debug
    )


    tx_info = ChromBridGE_tx.analyze_tx_alignment(
            read_aln_str=aln_info['read_aln'],
            ref1_aln_str=aln_info['ref1_aln'],
            ref2_aln_str=aln_info['ref2_aln'],
            breakpoints_read=aln_info['breakpoints_read'],
            breakpoints_ref1=aln_info['breakpoints_ref1'],
            breakpoints_ref2=aln_info['breakpoints_ref2'],
            read_path=aln_info['read_path'],
            ref1_cut_pos=ref1_cut_pos,
            ref2_cut_pos=ref2_cut_pos,
            min_num_bases_beyond_cut=min_num_bases_beyond_cut,
            min_num_bases_before_cut=min_num_bases_before_cut,
            mismatch_tolerance=mismatch_tolerance,
            gap_tolerance=gap_tolerance
    )

    return aln_info, tx_info


if __name__ == "__main__":
    main()
