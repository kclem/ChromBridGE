import Simulation
import argparse
import random

parser = argparse.ArgumentParser(description='Simulation for ChromBridGE',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--output_root', help='Output file to write results',required=True)
parser.add_argument('-n', help='Number of sequences to simulate',type=int,default=100)
parser.add_argument('--left_length_read',help='Number of bp left of translocation position in read',type=int,default = 30)
parser.add_argument('--right_length_read',help='Number of bp right of translocation position in read',type=int,default = 30)
parser.add_argument('--left_length_A',help='Number of bp left of cut site in sequence A',type=int,default = 70)
parser.add_argument('--right_length_A',help='Number of bp right of cut site in sequence A',type=int,default = 70)
parser.add_argument('--left_length_B',help='Number of bp left of cut site in sequence B',type=int,default = 70)
parser.add_argument('--right_length_B',help='Number of bp right of cut site in sequence B',type=int,default = 70)
parser.add_argument('--num_mutations_in_guide',help='Number of mutations in guide sequence',type=int,default = 1)
parser.add_argument('--debug',help='Whether to produce simplified simulated reads for debugging',action='store_true')

args = parser.parse_args()

(guide_left_A,guide_right_A,guide_left_B,guide_right_B) = Simulation.makeRandomGuides(args.num_mutations_in_guide,debug=args.debug)

(left_A,right_A,left_B,right_B) = Simulation.makeRandomAmplicons(guide_left_A,guide_right_A,guide_left_B,guide_right_B,args.left_length_A,args.right_length_A,args.left_length_B,args.right_length_B,debug=args.debug)

with open(args.output_root+".amps",'w') as fout:
    fout.write(guide_left_A+guide_right_A+"\n")
    fout.write(guide_left_B+guide_right_B+"\n")
    fout.write(left_A+right_A+"\n")
    fout.write(left_B+right_B+"\n")
    fout.write('left_length_read:'+str(args.left_length_read)+"\n")
    fout.write('right_length_read:'+str(args.right_length_read)+"\n")
    fout.write('left_length_A:'+str(args.left_length_A)+"\n")
    fout.write('right_length_A:'+str(args.right_length_A)+"\n")
    fout.write('left_length_B:'+str(args.left_length_B)+"\n")
    fout.write('right_length_B:'+str(args.right_length_B)+"\n")
    fout.write('num_mutations_in_guide:'+str(args.num_mutations_in_guide)+"\n")

with open (args.output_root+".fa",'w') as fout:
    for i in range(args.n):
        breakpoint = random.choice(range(-30,10))
        fout.write('@seq'+str(i)+'_bp'+str(breakpoint)+"\n")
        (read_A_left,read_A_right,read_B_left,read_B_right) = Simulation.makeSimulatedRead(left_A,right_A,left_B,right_B,args.left_length_read,args.right_length_read,breakpoint,debug=args.debug)
        read_seq = read_A_left+read_A_right+read_B_left+read_B_right
        qual = "H"*len(read_seq)
        fout.write(read_seq+"\n+\n"+qual+"\n")


print("Finished")

