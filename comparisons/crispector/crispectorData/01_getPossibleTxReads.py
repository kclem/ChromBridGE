import gzip

with open ('sites.txt','r') as fin:
    r1_loc = fin.readline().strip()
    r1_seq = fin.readline().strip()
    r2_loc = fin.readline().strip()
    r2_seq = fin.readline().strip()


primer1L = r1_seq[1:20].upper()
primer1R = r1_seq[len(r1_seq)-20:].upper()
primer2L = r2_seq[1:20].upper()
primer2R = r2_seq[len(r2_seq)-20:].upper()


total_count = 0
L1R1_count = 0
L2R2_count = 0
L1R2_count = 0
L2R1_count = 0

print('Mock')
with gzip.open('CRISPRessoPooled_on_EMX1_Mock/out.extendedFrags.fastq.gz','rt') as fin:
    while True:
        read_id = fin.readline()
        if not read_id:
            break
        read_seq = fin.readline().strip()
        read_plus = fin.readline()
        read_qual = fin.readline()

        total_count += 1


        left_20 = read_seq[1:20].upper()
        right_20 = read_seq[len(read_seq)-20:].upper()

#        print('comparing ' + left_20 + ' to ' + primer1L + ' and ' + right_20 + ' to ' + primer1R)
#        asdf()
        if left_20 == primer1L and right_20 == primer1R:
            L1R1_count += 1
        elif left_20 == primer2L and right_20 == primer2R:
            L2R2_count += 1
        elif left_20 == primer1L and right_20 == primer2R:
            L1R2_count += 1
        elif left_20 == primer2L and right_20 == primer1R:
            L2R1_count += 1

print('total: '+str(total_count))
print('L1R1: '+str(L1R1_count))
print('L2R2: '+str(L2R2_count))
print('L1R2: '+str(L1R2_count))
print('L2R1: '+str(L2R1_count))

total_count = 0
L1R1_count = 0
L2R2_count = 0
L1R2_count = 0
L2R1_count = 0

print('XT2P')
with gzip.open('CRISPRessoPooled_on_EMX1_XT2p/out.extendedFrags.fastq.gz','rt') as fin, open('01_getPossibleTxReads.py.reads','w') as fout:
    while True:
        read_id = fin.readline()
        if not read_id:
            break
        read_seq = fin.readline().strip()
        read_plus = fin.readline()
        read_qual = fin.readline()

        total_count += 1


        left_20 = read_seq[1:20].upper()
        right_20 = read_seq[len(read_seq)-20:].upper()

#        print('comparing ' + left_20 + ' to ' + primer1L + ' and ' + right_20 + ' to ' + primer1R)
#        asdf()
        if left_20 == primer1L and right_20 == primer1R:
            L1R1_count += 1
        elif left_20 == primer2L and right_20 == primer2R:
            L2R2_count += 1
        elif left_20 == primer1L and right_20 == primer2R:
            L1R2_count += 1
        elif left_20 == primer2L and right_20 == primer1R:
            L2R1_count += 1
        elif left_20 == primer1L and read_seq[len(read_seq)-4:].upper() != 'AACA' and 'ATATCTCTTTTACTGTG' not in read_seq:
            fout.write(read_seq+"\n")

print('total: '+str(total_count))
print('L1R1: '+str(L1R1_count))
print('L2R2: '+str(L2R2_count))
print('L1R2: '+str(L1R2_count))
print('L2R1: '+str(L2R1_count))
