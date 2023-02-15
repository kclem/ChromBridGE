from ChromBridGE import Simulation
from ChromBridGE import ChromBridGE_aln
import os
import random


total_count = 100
correct_count_pos_not_given = 0
correct_count_pos_given = 0
num_seq_reps = 100 #number of times each seq is repeated in the fastq file
num_ctl_reads = 5000

cmd_file = os.path.basename(__file__)+'.sh'
with open(cmd_file,'w') as cmd_out:
    os.makedirs('sims_same_guide', exist_ok=True)
    for i in range(total_count):
        this_root = "sims_same_guide/"+os.path.basename(__file__) + ".sim_"+str(i)

        (guide_left_A,guide_right_A,guide_left_B,guide_right_B) = Simulation.makeRandomGuides(num_mutations_in_guide=0,guide_length_bp=20,predicted_cut_position=17)

        (left_A,right_A,left_B,right_B) = Simulation.makeRandomAmplicons(
                guide_left_A = guide_left_A,
                guide_right_A = guide_right_A,
                guide_left_B = guide_left_B,
                guide_right_B = guide_right_B,
                left_length_A = 50,
                right_length_A = 50,
                left_length_B = 50,
                right_length_B = 50
                )

        seq_A = left_A + right_A
        cut_A = len(left_A)

        seq_B = left_B + right_B
        cut_B = len(left_B)

        with open(this_root+".info",'w') as fout:
            fout.write("cutA\t"+str(cut_A)+"\n")
            fout.write("seqA\t"+seq_A+"\n")
            fout.write("cutB\t"+str(cut_B)+"\n")
            fout.write("seqB\t"+seq_B+"\n")

        with open(this_root+".tx.fq",'w') as fout:
            for cut_ind in range(min(len(seq_A),len(seq_B))):
                read = seq_A[0:cut_ind] + seq_B[cut_ind:]
                qual = "H"*len(read)
                name = 'cut_%d'%cut_ind
                for rep in range(num_seq_reps):
                    this_read = read
    #                if random.random() < 0.01:
    #                    this_nuc = random.choice(['A','T','C','G'])
    #                    this_read = this_read[0:29] + this_nuc + this_read[31:]
                    fout.write("@%s_rep_%d\n"%(name,rep))
                    fout.write("%s\n+\n%s\n"%(this_read,qual))

            read = seq_A
            qual = "H"*len(read)
            name = 'ref1'
            for rep in range(num_ctl_reads):
                this_read = read
    #            if random.random() < 0.01:
    #                this_nuc = random.choice(['A','T','C','G'])
    #                this_read = this_read[0:29] + this_nuc + this_read[31:]
                fout.write("@%s_rep_%d\n"%(name,rep))
                fout.write("%s\n+\n%s\n"%(this_read,qual))

            #write some with indels
            read = seq_A
            name = 'ref1_mod'
            for rep in range(num_ctl_reads):
                this_read = read
                this_del_start = int(20+random.random()*(len(this_read)-40))
                this_del_end = int(20+random.random()*(len(this_read)-40))
                this_del_start,this_del_end = sorted([this_del_start,this_del_end])
                this_read = this_read[0:this_del_start] + this_read[this_del_end:]
                this_qual = "H"*len(this_read)
                fout.write("@%s_rep_%d\n"%(name,rep))
                fout.write("%s\n+\n%s\n"%(this_read,this_qual))

            read = seq_B
            qual = "H"*len(read)
            name = 'ref2'
            for rep in range(num_ctl_reads):
                this_read = read
    #            if random.random() < 0.01:
    #                this_nuc = random.choice(['A','T','C','G'])
    #                this_read = this_read[0:29] + this_nuc + this_read[31:]
                fout.write("@%s_rep_%d\n"%(name,rep))
                fout.write("%s\n+\n%s\n"%(this_read,qual))

            #write some with indels
            read = seq_B
            name = 'ref2_mod'
            for rep in range(num_ctl_reads):
                this_read = read
                this_del_start = int(20+random.random()*(len(this_read)-40))
                this_del_end = int(20+random.random()*(len(this_read)-40))
                this_del_start,this_del_end = sorted([this_del_start,this_del_end])
                this_read = this_read[0:this_del_start] + this_read[this_del_end:]
                this_qual = "H"*len(this_read)
                fout.write("@%s_rep_%d\n"%(name,rep))
                fout.write("%s\n+\n%s\n"%(this_read,this_qual))

        with open(this_root+".mock.fq",'w') as fout:
            read = seq_A
            qual = "H"*len(read)
            name = 'ref1'
            for rep in range(num_ctl_reads):
                this_read = read
    #            if random.random() < 0.01:
    #                this_nuc = random.choice(['A','T','C','G'])
    #                this_read = this_read[0:29] + this_nuc + this_read[31:]
                fout.write("@%s_rep_%d\n"%(name,rep))
                fout.write("%s\n+\n%s\n"%(this_read,qual))

            read = seq_B
            qual = "H"*len(read)
            name = 'ref2'
            for rep in range(num_ctl_reads):
                this_read = read
    #            if random.random() < 0.01:
    #                this_nuc = random.choice(['A','T','C','G'])
    #                this_read = this_read[0:29] + this_nuc + this_read[31:]
                fout.write("@%s_rep_%d\n"%(name,rep))
                fout.write("%s\n+\n%s\n"%(this_read,qual))

        with open(this_root+".crispectorConfig.csv",'w') as fout:
            fout.write('SiteName,AmpliconReference,gRNA,OnTarget\n')
            fout.write(",".join([
                'ref1',
                seq_A,
                guide_left_A+guide_right_A,
                'True'
                ])+"\n")
            fout.write(",".join([
                'ref2',
                seq_B,
                guide_left_B+guide_right_B,
                'False'
                ])+"\n")

        with open(this_root+".crispector.sh",'w') as fout:
            fout.write('crispector -t_r1 '+this_root+'.tx.fq -m_r1 '+this_root+'.mock.fq -c '+this_root+'.crispectorConfig.csv -o CRISPECTOR_OUT/'+this_root+' > ' + this_root + '.crispector.sh.log')
        cmd_out.write('bash ' + this_root+".crispector.sh\n")


print("Finished")
print('Run bash ' + cmd_file)

