import numpy as np

def nw_breakpoint(read_seq,ref1_seq,ref2_seq,match_score=3,mismatch_score=-1,gap_score=-2,perimeter_gap_extension_score=0,jump_score=-3,cut_pos_incentive_score=0,ref1_cut_pos=None,ref2_cut_pos=None):
    """
    Computes the optimal alignment of a read to two seqences, locating the optimal break between the two reads.
    The alignment score will be the sum of the match, mismatch, gap, and jump scores.

    params:
        read_seq: read to align to the other two sequences
        ref1_seq: first sequence to align to
        ref2_seq: second sequence to align to
        match_score: score for adding a match in alignment (positive)
        mismatch_score: score for adding a mismatch in alignment
        gap_score: score for adding a gap in alignment
        perimeter_gap_extension_score: score for adding a gap in the first/last column/row, corresponding to gaps at the beginning or ends of sequences
        jump_score: score for jumping between ref1 and ref2
        cut_pos_incentive_score: score for jumping at a predicted cut position
        ref1_cut_pos: position of predicted cut site in ref1
        ref2_cut_pos: position of predicted cut site in ref2

    returns:
        read: aligned sequence of read
        ref1: sequence of ref1 aligned to read
        ref2: sequence of ref2 aligned to read
        breakpoints: indices in read of the breakpoints discovered
            tuple of (read, ref1, ref2) indices
            read indices show the positions of all breakpoints
            ref1 and ref2 indices show the positions at which the optimal alignment switches to that reference.

    """

    pointer_match = 1
    pointer_gap_ref = 2
    pointer_gap_read = 3
    pointer_jump = 4

    #read is columns, refs are rows
    len_read = len(read_seq)
    len_ref1 = len(ref1_seq)
    len_ref2 = len(ref2_seq)

    #set jump incentive arrays (where jumping is less penalized)
    jump_incentive_ref1 = np.zeros(len_ref1 + 1)
    if ref1_cut_pos is not None:
        jump_incentive_ref1[ref1_cut_pos+1] = cut_pos_incentive_score
    jump_incentive_ref2 = np.zeros(len_ref2 + 1)
    if ref2_cut_pos is not None:
        jump_incentive_ref2[ref2_cut_pos+1] = cut_pos_incentive_score

    # Optimal score at each possible pair of characters.
    score1 = np.zeros((len_ref1 + 1, len_read + 1))
    score1[:,0] = np.linspace(0, len_ref1*perimeter_gap_extension_score, len_ref1 + 1)
    score1[0,:] = np.linspace(0, len_read*perimeter_gap_extension_score, len_read + 1)

    score2 = np.zeros((len_ref2 + 1, len_read + 1))
    score2[:,0] = np.linspace(0, len_ref2*perimeter_gap_extension_score, len_ref2 + 1)
    score2[0,:] = np.linspace(0, len_read*perimeter_gap_extension_score, len_read + 1)

    #this hack keeps references from sliding all the way to the end
    score1[1,0] = gap_score
    score1[0,1] = gap_score

    score2[1,0] = gap_score
    score2[0,1] = gap_score

    # Pointers to trace through an optimal aligment.
    pointer1 = np.zeros((len_ref1 + 1, len_read + 1))
    pointer1[:,0] = pointer_gap_read
    pointer1[0,:] = pointer_gap_ref

    pointer2 = np.zeros((len_ref2 + 1, len_read + 1))
    pointer2[:,0] = pointer_gap_read
    pointer2[0,:] = pointer_gap_ref


    #keep track of where the maximum is for jumping
    colmaxes1 = np.copy(score1[0,:])
    colmaxesInd1 = [0]*(len_read + 1)
    colmaxes2 = np.copy(score2[0,:])
    colmaxesInd2 = [0]*(len_read + 1)


    for idx_read in range(1,len_read+1):


        #do table 1 for this column
        for idx_ref1 in range(1,len_ref1+1):
            this_match_score = 0
            if read_seq[idx_read-1] == ref1_seq[idx_ref1-1]:
                this_match_score = score1[idx_ref1-1,idx_read-1] + match_score
            else:
                this_match_score = score1[idx_ref1-1,idx_read-1] + mismatch_score

            this_gap_up_score = gap_score
            if idx_read == len_read: #if the last column, no gap penalty
                this_gap_up_score = 0
                if idx_ref1 ==len_ref1:#except for the bottom right cell with full penalty to shift alignments to the middle
                    this_gap_up_score = gap_score

            this_gap_left_score = gap_score
            if idx_ref1 == len_ref1: #if the last row, no gap penalty
                this_gap_left_score = 0
                if idx_read == len_read:
                    this_gap_left_score = gap_score

            this_ref_gap_score = score1[idx_ref1,idx_read-1] + this_gap_left_score
            this_read_gap_score = score1[idx_ref1-1,idx_read] + this_gap_up_score
            this_jump_score = colmaxes2[idx_read-1] + jump_score + jump_incentive_ref1[idx_ref1]
            tmax = np.max([this_match_score,this_ref_gap_score,this_read_gap_score,this_jump_score])
#            print('idx_read: ' + str(idx_read))
#            print('idx_ref1: ' + str(idx_ref1))
#            print('testing ' + read_seq[idx_read-1]+ ' and ' + ref1_seq[idx_ref1-1])
#            print('this_match_score: ' + str(this_match_score))
#            print('this_ref_gap_score: ' + str(this_ref_gap_score))
#            print('this_read_gap_score: ' + str(this_read_gap_score))
#            print('this_jump_score: ' + str(this_jump_score))
#            print('tmax: ' + str(tmax))


            score1[idx_ref1,idx_read] = tmax
            if this_match_score == tmax:
                pointer1[idx_ref1,idx_read] = pointer_match
            if this_ref_gap_score == tmax:
                pointer1[idx_ref1,idx_read] = pointer_gap_ref
            if this_read_gap_score == tmax:
                pointer1[idx_ref1,idx_read] = pointer_gap_read
            if this_jump_score == tmax:
                pointer1[idx_ref1,idx_read] = pointer_jump

#            print('comparing tmax ' + str(tmax) + ' to val: ' + str(colmaxes1[idx_read]))
            if tmax > colmaxes1[idx_read]:
                colmaxes1[idx_read] = tmax
                colmaxesInd1[idx_read] = idx_ref1

        for idx_ref2 in range(1,len_ref2+1):
            this_match_score = 0
            if read_seq[idx_read-1] == ref2_seq[idx_ref2-1]:
                this_match_score = score2[idx_ref2-1,idx_read-1] + match_score
            else:
                this_match_score = score2[idx_ref2-1,idx_read-1] + mismatch_score
            this_gap_up_score = gap_score
            if idx_read == len_read: #if the last column, no gap penalty
                this_gap_up_score = 0
                if idx_ref2 ==len_ref2:#except for the bottom right cell with full penalty to shift alignments to the middle
                    this_gap_up_score = gap_score

            this_gap_left_score = gap_score
            if idx_ref2 == len_ref2: #if the last row, no gap penalty
                this_gap_left_score = 0
                if idx_read == len_read:
                    this_gap_left_score = gap_score

            this_ref_gap_score = score2[idx_ref2,idx_read-1] + this_gap_left_score
            this_read_gap_score = score2[idx_ref2-1,idx_read] + this_gap_up_score
            this_jump_score = colmaxes1[idx_read-1] + jump_score + jump_incentive_ref2[idx_ref2]
            tmax = np.max([this_match_score,this_ref_gap_score,this_read_gap_score,this_jump_score])
#            print('idx_read: ' + str(idx_read))
#            print('idx_ref2: ' + str(idx_ref2))
#            print('testing ' + read_seq[idx_read-1]+ ' and ' + ref2_seq[idx_ref2-1])
#            print('this_match_score: ' + str(this_match_score))
#            print('this_ref_gap_score: ' + str(this_ref_gap_score))
#            print('this_read_gap_score: ' + str(this_read_gap_score))
#            print('this_jump_score: ' + str(this_jump_score))
#            print('tmax: ' + str(tmax))

            score2[idx_ref2,idx_read] = tmax
            if this_match_score == tmax:
                pointer2[idx_ref2,idx_read] = pointer_match
            if this_ref_gap_score == tmax:
                pointer2[idx_ref2,idx_read] = pointer_gap_ref
            if this_read_gap_score == tmax:
                pointer2[idx_ref2,idx_read] = pointer_gap_read
            if this_jump_score == tmax:
                pointer2[idx_ref2,idx_read] = pointer_jump

#            print('comparing tmax ' + str(tmax) + ' to val: ' + str(colmaxes2[idx_read]))
            if tmax > colmaxes2[idx_read]:
                colmaxes2[idx_read] = tmax
                colmaxesInd2[idx_read] = idx_ref2



#    print('score1:')
#    print(score1)
#    print('pointer1:')
#    print(pointer1)
#    print('colmaxes1:')
#    print(colmaxes1)
#
#    print('score2:')
#    print(score2)
#    print('pointer2:')
#    print(pointer2)
#    print('colmaxes2:')
#    print(colmaxes2)

    # Trace through an optimal alignment.
    idx_read = len_read
    idx_ref = len_ref1
    curr_matrix = 1
    if score1[len_ref1,len_read] < score2[len_ref2,len_read]:
        idx_ref = len_ref2
        curr_matrix = 2

    final_read_aln = []
    final_ref1_aln = []
    final_ref2_aln = []
    breakpoints_read = []
    breakpoints_ref1 = []
    breakpoints_ref2 = []
    while idx_read > 0 or idx_ref > 0:
#        print('curr matrix: ' + str(curr_matrix))
#        print('idx_read: ' + str(idx_read))
#        print('idx_ref: ' + str(idx_ref))
        if curr_matrix == 1:
#            print('pointer1: ' + str(pointer1[idx_ref,idx_read]))
            if pointer1[idx_ref,idx_read] == pointer_match:
                final_read_aln.append(read_seq[idx_read-1])
                final_ref1_aln.append(ref1_seq[idx_ref-1])
                final_ref2_aln.append(" ")
                idx_read -= 1
                idx_ref -= 1
            elif pointer1[idx_ref,idx_read] == pointer_gap_ref:
                final_read_aln.append(read_seq[idx_read-1])
                final_ref1_aln.append("-")
                final_ref2_aln.append(" ")
                idx_read -= 1
            elif pointer1[idx_ref,idx_read] == pointer_gap_read:
                final_read_aln.append("-")
                final_ref1_aln.append(ref1_seq[idx_ref-1])
                final_ref2_aln.append(" ")
                idx_ref -= 1
            elif pointer1[idx_ref,idx_read] == pointer_jump:
                final_read_aln.append(read_seq[idx_read-1])
                final_ref1_aln.append(ref1_seq[idx_ref-1])
                final_ref2_aln.append(" ")
                idx_read -= 1
                curr_matrix = 2
                breakpoints_ref2.append(idx_ref-1)
                idx_ref = colmaxesInd2[idx_read]
                breakpoints_read.append(idx_read)
                breakpoints_ref1.append(idx_ref)

        elif curr_matrix == 2:
#            print('pointer2: ' + str(pointer2[idx_ref,idx_read]))
            if pointer2[idx_ref,idx_read] == pointer_match:
                final_read_aln.append(read_seq[idx_read-1])
                final_ref1_aln.append(" ")
                final_ref2_aln.append(ref2_seq[idx_ref-1])
                idx_read -= 1
                idx_ref -= 1
            elif pointer2[idx_ref,idx_read] == pointer_gap_ref:
                final_read_aln.append(read_seq[idx_read-1])
                final_ref1_aln.append(" ")
                final_ref2_aln.append("-")
                idx_read -= 1
            elif pointer2[idx_ref,idx_read] == pointer_gap_read:
                final_read_aln.append("-")
                final_ref1_aln.append(" ")
                final_ref2_aln.append(ref2_seq[idx_ref-1])
                idx_ref -= 1
            elif pointer2[idx_ref,idx_read] == pointer_jump:
                final_read_aln.append(read_seq[idx_read-1])
                final_ref1_aln.append(" ")
                final_ref2_aln.append(ref2_seq[idx_ref-1])
                idx_read -= 1
                curr_matrix = 1
                breakpoints_ref1.append(idx_ref-1)
                idx_ref = colmaxesInd1[idx_read]
                breakpoints_read.append(idx_read)
                breakpoints_ref2.append(idx_ref)

    # Reverse the strings.
    final_read_aln = ''.join(final_read_aln)[::-1]
    final_ref1_aln = ''.join(final_ref1_aln)[::-1]
    final_ref2_aln = ''.join(final_ref2_aln)[::-1]
#    print('read aln : ' + final_read_aln)
#    print('read ref1: ' + final_ref1_aln)
#    print('read ref2: ' + final_ref2_aln)
    return(final_read_aln,final_ref1_aln,final_ref2_aln,
            (breakpoints_read,breakpoints_ref1,breakpoints_ref2))

if __name__ == "__main__":
    print('Performing tests..')
    read,ref1,ref2,breakpoints = nw_breakpoint(
                'AAATTT',
                'AAAAAA',
                'TTTTTT')
#    print('read: ' + read)
#    print('ref1: ' + ref1)
#    print('ref2: ' + ref2)
    if read !=  'AAATTT' or \
        ref1 != 'AAA   ' or \
        ref2 != '   TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'TTTAAA',
                'AAAAAA',
                'TTTTTT')
    if read !=  'TTTAAA' or \
        ref1 != '   AAA' or \
        ref2 != 'TTT   ' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'AAAATTT',
                'AAAAAA',
                'TTTTTTT')
    if read !=  'AAAATTT' or \
        ref1 != 'AAAA   ' or \
        ref2 != '    TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'AAAATTT',
                'AAAAA',
                'TTTTTTT')
    if read !=  'AAAATTT' or \
        ref1 != 'AAAA   ' or \
        ref2 != '    TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'AAAATTT',
                'AAA',
                'TTTTTTT')
    if read !=  'AAAATTT' or \
        ref1 != 'AAA    ' or \
        ref2 != '   TTTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'AACCTTGG',
                'AAA',
                'CCCTGG')
    if read !=  'AACCTTGG' or \
        ref1 != 'AA      ' or \
        ref2 != '  CCCTGG':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'ACTGACTGACTG',
                'ACTGACTGACTG',
                'CCCTGG')
    if read !=  'ACTGACTGACTG' or \
        ref1 != 'ACTGACTGACTG' or \
        ref2 != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'ACTGACTGACTG',
                'ACTGCTGACTG',
                'CCCTGG')
    if read !=  'ACTGACTGACTG' or \
        ref1 != 'ACTG-CTGACTG' or \
        ref2 != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'ACTGACTGACTG',
                'CTGCTGACT',
                'CCCTGG')
    if read !=  'ACTGACTGACTG' or \
        ref1 != '-CTG-CTGACT-' or \
        ref2 != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints = nw_breakpoint(
                'GTCGACTGACTG',
                'GTCAGTCAGTCA',
                'ACTGACTGACTG')
    if read !=  'GTCGACTGACTG' or \
        ref1 != 'GTC         ' or \
        ref2 != '   GACTGACTG':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    print("Tests passed")
