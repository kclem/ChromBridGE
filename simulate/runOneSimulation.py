import Simulation
import argparse
import random
import ChromBridGE_aln

parser = argparse.ArgumentParser(description='Simulation for ChromBridGE',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--left_length_read',help='Number of bp left of translocation position in read',type=int,default = 30)
parser.add_argument('--right_length_read',help='Number of bp right of translocation position in read',type=int,default = 30)
parser.add_argument('--left_length_A',help='Number of bp left of cut site in sequence A',type=int,default = 70)
parser.add_argument('--right_length_A',help='Number of bp right of cut site in sequence A',type=int,default = 70)
parser.add_argument('--left_length_B',help='Number of bp left of cut site in sequence B',type=int,default = 70)
parser.add_argument('--right_length_B',help='Number of bp right of cut site in sequence B',type=int,default = 70)
parser.add_argument('--num_mutations_in_guide',help='Number of mutations in guide sequence',type=int,default = 1)
parser.add_argument('--match_score', help='Match score for alignment',default=1)
parser.add_argument('--mismatch_score', help='Mismatch score for alignment',default=-1)
parser.add_argument('--gap_score', help='Gap score for alignment',default=-2)
parser.add_argument('--jump_score', help='Jump score for alignment',default=-3)
parser.add_argument('--cut_pos_incentive_score', help='Incentive for jumping at a predicted cut site',default=0)
parser.add_argument('--debug',help='Whether to produce simplified simulated reads for debugging',action='store_true')

args = parser.parse_args()

(guide_left_A,guide_right_A,guide_left_B,guide_right_B) = Simulation.makeRandomGuides(args.num_mutations_in_guide,debug=args.debug)
print ("guideA: " + guide_left_A + "|" + guide_right_A)
print ("guideB: " + guide_left_B + "|" + guide_right_B)

(left_A,right_A,left_B,right_B) = Simulation.makeRandomAmplicons(guide_left_A,guide_right_A,guide_left_B,guide_right_B,args.left_length_A,args.right_length_A,args.left_length_B,args.right_length_B,debug=args.debug)

(read_A_left,read_A_right,read_B_left,read_B_right) = Simulation.makeSimulatedRead(left_A,right_A,left_B,right_B,args.left_length_read,args.right_length_read,0,debug=args.debug)

read_seq = read_A_left+read_A_right+read_B_left+read_B_right
seq_A = left_A + right_A
cut_A = len(left_A)

seq_B = left_B + right_B
cut_B = len(left_B)

read_aln,refA_aln,refB_aln,breakpoints = ChromBridGE_aln.nw_breakpoint(read_seq,seq_A,seq_B,
        match_score=args.match_score,
        mismatch_score=args.mismatch_score,
        gap_score=args.gap_score,
        jump_score=args.jump_score,
        cut_pos_incentive_score=args.cut_pos_incentive_score,
        ref1_cut_pos=cut_A,
        ref2_cut_pos=cut_B)



print("read_seq:"+read_seq)
nucStr = "  "
for i in range(15):
    nucStr += "        "+str(10*i)
print(nucStr)
print("seq_A:   "+seq_A)
print("seq_A:   "+refA_aln)
print("read_aln:"+read_aln)
print("seq_B:   "+refB_aln)
print("seq_B:   "+seq_B)
print('breakpoints:'+str(breakpoints))


print("Finished")

