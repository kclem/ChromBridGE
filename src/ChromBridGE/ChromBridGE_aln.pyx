import numpy as np
import cython

cdef int mymax4(int s1, int s2, int s3, int s4):
    cdef int mymax = s1
    if s2 > mymax:
        mymax = s2
    if s3 > mymax:
        mymax = s3
    if s4 > mymax:
        mymax = s4
    return mymax


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.initializedcheck(False)
cpdef nw_breakpoint(str read_seq_py,
                    str ref1_seq_py,
                    str ref2_seq_py,
                    int match_score=3,
                    int mismatch_score=-1,
                    int gap_score=-2,
                    int perimeter_gap_extension_score=0,
                    int jump_score=-12, # four matches
                    int cut_pos_jump_incentive_score=1,
                    ref1_cut_pos=None,
                    ref2_cut_pos=None,
		    bint debug=False):
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
        cut_pos_jump_incentive_score: score incentive for jumping at a predicted cut position
        ref1_cut_pos: position of predicted cut site in ref1
        ref2_cut_pos: position of predicted cut site in ref2

    returns:
        dict containing:
        read_aln: aligned sequence of read
        ref1_aln: sequence of ref1 aligned to read
        ref2_aln: sequence of ref2 aligned to read
        breakpoints_read: indices in read of the breakpoints discovered
        breakpoints_ref1: indices in ref1 of breakpoints in the optimal alignment
        breakpoints_ref2: indices in ref2 of breakpoints in the optimal alignment
        aln_score: score of alignment
        read_path: index of ref that the read is aligned to, corresponding to the break points (there will be len(breakpoints)+1 items in read_path)
    """

    score_type = np.intc

    cdef const unsigned char[:] read_seq = read_seq_py.encode()
    cdef const unsigned char[:] ref1_seq = ref1_seq_py.encode()
    cdef const unsigned char[:] ref2_seq = ref2_seq_py.encode()

    cdef char pointer_match = 1
    cdef char pointer_gap_read = 2 #move from left
    cdef char pointer_gap_ref = 3 #move from up
    cdef char pointer_jump = 4

    #read is columns, refs are rows
    cdef int len_read = len(read_seq)
    cdef int len_ref1 = len(ref1_seq)
    cdef int len_ref2 = len(ref2_seq)

    cdef int prefer_cut_ref1_idx = len_ref1 #if above/lower than this idx prefer match/mismatch over jump. If below/greater than this idx prefer jump over match/mismatch
    cdef int prefer_cut_ref2_idx = len_ref2

    #set jump incentive arrays (where jumping is less penalized at cut sites)
    cdef int[:] jump_incentive_ref1 = np.zeros(len_ref1 + 1, dtype=score_type)
    if ref1_cut_pos is not None:
        jump_incentive_ref1[ref1_cut_pos] = cut_pos_jump_incentive_score
        prefer_cut_ref1_idx = ref1_cut_pos + 1
    cdef int[:] jump_incentive_ref2 = np.zeros(len_ref2 + 1, dtype=score_type)
    if ref2_cut_pos is not None:
        jump_incentive_ref2[ref2_cut_pos] = cut_pos_jump_incentive_score
        prefer_cut_ref2_idx = ref2_cut_pos + 1

    # Optimal score at each possible pair of characters.
    score1_py = np.zeros((len_ref1 + 1, len_read + 1), dtype=score_type)
    score1_py[:,0] = np.linspace(0, len_ref1*perimeter_gap_extension_score, len_ref1 + 1)
    score1_py[0,:] = np.linspace(0, len_read*perimeter_gap_extension_score, len_read + 1)

    score2_py = np.zeros((len_ref2 + 1, len_read + 1), dtype=score_type)
    score2_py[:,0] = np.linspace(0, len_ref2*perimeter_gap_extension_score, len_ref2 + 1)
    score2_py[0,:] = np.linspace(0, len_read*perimeter_gap_extension_score, len_read + 1)

    #this hack keeps references from sliding all the way to the end
    score1_py[1,0] = gap_score
    score1_py[0,1] = gap_score

    score2_py[1,0] = gap_score
    score2_py[0,1] = gap_score

    # Pointers to trace through an optimal aligment.
    pointer1_py = np.zeros((len_ref1 + 1, len_read + 1),dtype=np.byte)
    pointer1_py[:,0] = pointer_gap_ref
    pointer1_py[0,:] = pointer_gap_read

    pointer2_py = np.zeros((len_ref2 + 1, len_read + 1),dtype=np.byte)
    pointer2_py[:,0] = pointer_gap_ref
    pointer2_py[0,:] = pointer_gap_read


    #keep track of where the maximum is for jumping
    #colmaxesInd keep track of the index (row) which had the max value. It's an array because there could be multiple values
    cdef int[:] colmaxes1 = np.copy(score1_py[0,:])
    cdef int[:] colmaxesInd1 = np.zeros(len_read + 1,dtype=score_type)
    cdef int[:] colmaxes2 = np.copy(score2_py[0,:])
    cdef int[:] colmaxesInd2 = np.zeros(len_read + 1,dtype=score_type)

    #crazy python bug caught here:
    # incorrect:
    # colmaxesInd2 = [[0]]*(len_read + 1)
    # because this produces referencs to a single array...
# >>> d = [5]*5
# >>> d
# [5, 5, 5, 5, 5]
# >>> d[1] += 5
# >>> d
# [5, 10, 5, 5, 5]
# >>> d = [[5]]*5
# >>> d
# [[5], [5], [5], [5], [5]]
# >>> d[1][0] += 5
# >>> d
# [[10], [10], [10], [10], [10]]



    cdef int idx_read, idx_ref1, idx_ref2
    cdef int this_match_or_mismatch_score, this_match_score, this_gap_up_score, this_gap_left_score, this_read_gap_score, this_ref_gap_score, this_jump_score, tmax, tmax_plus_jump
    cdef int[:,:] score1 = score1_py #cython memory view
    cdef int[:,:] score2 = score2_py
    cdef char[:,:] pointer1 =  pointer1_py
    cdef char[:,:] pointer2 =  pointer2_py

    for idx_read in range(1,len_read+1):


        #do table 1 for this column
        for idx_ref1 in range(1,len_ref1+1):
            this_match_or_mismatch_score = mismatch_score #keep this separate for the jump score below
            if read_seq[idx_read-1] == ref1_seq[idx_ref1-1]:
                this_match_or_mismatch_score = match_score
            this_match_score = score1[idx_ref1-1,idx_read-1] + this_match_or_mismatch_score

            this_gap_up_score = gap_score
            if idx_read == len_read: #if the last column, no gap penalty
                this_gap_up_score = perimeter_gap_extension_score
                if idx_ref1 ==len_ref1:#except for the bottom right cell with full penalty to shift alignments to the middle
                    this_gap_up_score = gap_score

            this_gap_left_score = gap_score
            if idx_ref1 == len_ref1: #if the last row, no gap penalty
                this_gap_left_score = perimeter_gap_extension_score
                if idx_read == len_read:
                    this_gap_left_score = gap_score

            this_read_gap_score = score1[idx_ref1,idx_read-1] + this_gap_left_score
            this_ref_gap_score = score1[idx_ref1-1,idx_read] + this_gap_up_score
            #technically, a 'jump' is a 'jump and consume' so it's two steps, but because you would never have two jumps in a row, we can consume a base from the read sequence and do two steps (jump and consume) in one step based on the max values from the last column
            this_jump_score = colmaxes2[idx_read-1] + jump_score + jump_incentive_ref1[idx_ref1 -1] + this_match_or_mismatch_score

            tmax = mymax4(this_match_score,this_ref_gap_score,this_read_gap_score,this_jump_score)
#            print('SCORE1')
#            print('idx_read: ' + str(idx_read))
#            print('idx_ref1: ' + str(idx_ref1))
#            print('testing ' + read_seq[idx_read-1]+ ' and ' + ref1_seq[idx_ref1-1])
#            print('this_match_score (diagonal): ' + str(this_match_score))
#            print('this_ref_gap_score (left): ' + str(this_ref_gap_score))
#            print('this_read_gap_score (up): ' + str(this_read_gap_score))
#            print('this_jump_score: ' + str(this_jump_score))
#            print('this_jump_score is colmax=' + str(colmaxes2[idx_read-1]) + ' + jumpscore='+ str(jump_score) + ' + incentive=' + str(jump_incentive_ref1[idx_ref1-1]) + ' + match_or_mismatch=' + str(this_match_or_mismatch_score))
#            print('tmax: ' + str(tmax))
#            print('score1:')
#            score1[idx_ref1,idx_read] = '-100'
#            print(score1)
#            print('pointer1:')
#            print(pointer1)


            score1[idx_ref1,idx_read] = tmax
            if idx_ref1 < prefer_cut_ref1_idx:
                # prefer jump over match/mismatch if above
                if this_jump_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_jump
                elif this_match_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_match
                elif this_ref_gap_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_gap_ref
                elif this_read_gap_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_gap_read
            else:
                # prefer match/mismatch/gap if below
                if this_match_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_match
                elif this_ref_gap_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_gap_ref
                elif this_read_gap_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_gap_read
                elif this_jump_score == tmax:
                    pointer1[idx_ref1,idx_read] = pointer_jump

#            print('comparing tmax ' + str(tmax) + ' to val: ' + str(colmaxes1[idx_read]))
            tmax_plus_jump = tmax + jump_incentive_ref1[idx_ref1]
            if tmax_plus_jump > colmaxes1[idx_read]:
                colmaxes1[idx_read] = tmax_plus_jump
                colmaxesInd1[idx_read] = idx_ref1


        for idx_ref2 in range(1,len_ref2+1):
            this_match_or_mismatch_score = mismatch_score
            if read_seq[idx_read-1] == ref2_seq[idx_ref2-1]:
                this_match_or_mismatch_score = match_score
            this_match_score = score2[idx_ref2-1,idx_read-1] + this_match_or_mismatch_score

            this_gap_up_score = gap_score
            if idx_read == len_read: #if the last column, no gap penalty
                this_gap_up_score = perimeter_gap_extension_score
                if idx_ref2 ==len_ref2:#except for the bottom right cell with full penalty to shift alignments to the middle
                    this_gap_up_score = gap_score

            this_gap_left_score = gap_score
            if idx_ref2 == len_ref2: #if the last row, no gap penalty
                this_gap_left_score = perimeter_gap_extension_score
                if idx_read == len_read:
                    this_gap_left_score = gap_score

            this_read_gap_score = score2[idx_ref2,idx_read-1] + this_gap_left_score
            this_ref_gap_score = score2[idx_ref2-1,idx_read] + this_gap_up_score
            this_jump_score = colmaxes1[idx_read-1] + jump_score + jump_incentive_ref2[idx_ref2 - 1] + this_match_or_mismatch_score


            tmax = mymax4(this_match_score,this_ref_gap_score,this_read_gap_score,this_jump_score)
#            print('SCORE2')
#            print('idx_read: ' + str(idx_read))
#            print('idx_ref2: ' + str(idx_ref2))
#            print('testing ' + read_seq[idx_read-1]+ ' and ' + ref2_seq[idx_ref2-1])
#            print('this_match_score: ' + str(this_match_score))
#            print('this_ref_gap_score (left): ' + str(this_ref_gap_score))
#            print('this_read_gap_score (up): ' + str(this_read_gap_score))
#            print('this_jump_score: ' + str(this_jump_score))
#            print('this_jump_score is colmax=' + str(colmaxes1[idx_read-1]) + ' + jumpscore='+ str(jump_score) + ' + incentive=' + str(jump_incentive_ref2[idx_ref2 - 1]) + ' + match_or_mismatch=' + str(this_match_or_mismatch_score))
#            print('tmax: ' + str(tmax))
#            print('score2:')
#            score2[idx_ref2,idx_read] = '-100'
#            print(score2)
#            print('pointer2:')
#            print(pointer2)

            score2[idx_ref2,idx_read] = tmax
            if idx_ref2 < prefer_cut_ref2_idx:
                if this_jump_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_jump
                elif this_match_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_match
                elif this_ref_gap_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_gap_ref
                elif this_read_gap_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_gap_read
            else:
                if this_match_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_match
                elif this_ref_gap_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_gap_ref
                elif this_read_gap_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_gap_read
                elif this_jump_score == tmax:
                    pointer2[idx_ref2,idx_read] = pointer_jump

#            print('comparing tmax ' + str(tmax) + ' to val: ' + str(colmaxes2[idx_read]))
            tmax_plus_jump = tmax + jump_incentive_ref2[idx_ref2]
            if tmax_plus_jump > colmaxes2[idx_read]:
                colmaxes2[idx_read] = tmax_plus_jump
                colmaxesInd2[idx_read] = idx_ref2


    if debug:
        np.set_printoptions(threshold=np.inf)
        print('jump_incentive_ref1')
        print(np.array(jump_incentive_ref1))
        print('score1:')
        print(np.array(score1))
        print('pointer1:')
        print(np.array(pointer1))
        print('colmaxes1:')
        print(np.array(colmaxes1))
        print('colmaxesInd1:')
        print(np.array(colmaxesInd1))

        print('jump_incentive_ref2')
        print(np.array(jump_incentive_ref2))
        print('score2:')
        print(np.array(score2))
        print('pointer2:')
        print(np.array(pointer2))
        print('colmaxes2:')
        print(np.array(colmaxes2))
        print('colmaxesInd2:')
        print(np.array(colmaxesInd2))

    # Trace through an optimal alignment.
    idx_read = len_read
    idx_ref = len_ref1
    curr_matrix = 1
    max_score = score1[len_ref1,len_read]
    if score1[len_ref1,len_read] < score2[len_ref2,len_read]:
        idx_ref = len_ref2
        curr_matrix = 2
        max_score = score2[len_ref2,len_read]

    final_read_aln = []
    final_ref1_aln = []
    final_ref2_aln = []
    breakpoints_read = []
    breakpoints_ref1 = []
    breakpoints_ref2 = []
    read_path = [curr_matrix]

    while idx_read > 0 or idx_ref > 0:
#        print('curr matrix: ' + str(curr_matrix))
#        print('idx_read: ' + str(idx_read))
#        print('idx_ref: ' + str(idx_ref))
#        print('     breakpoints_read: ' + str(breakpoints_read))
#        print('     final_read_aln: ' + ''.join(final_read_aln[::-1]))
#        print('     final_ref1_aln: ' + ''.join(final_ref1_aln[::-1]))
#        print('     final_ref2_aln: ' + ''.join(final_ref2_aln[::-1]))
        if curr_matrix == 1:
#            print('     pointer1: ' + str(pointer1[idx_ref,idx_read]))
            if pointer1[idx_ref,idx_read] == pointer_match:
                final_read_aln.append(chr(read_seq[idx_read-1]))
                final_ref1_aln.append(chr(ref1_seq[idx_ref-1]))
                final_ref2_aln.append(" ")
                idx_read -= 1
                idx_ref -= 1
            elif pointer1[idx_ref,idx_read] == pointer_gap_read:
                final_read_aln.append(chr(read_seq[idx_read-1]))
                final_ref1_aln.append("-")
                final_ref2_aln.append(" ")
                idx_read -= 1
            elif pointer1[idx_ref,idx_read] == pointer_gap_ref:
                final_read_aln.append("-")
                final_ref1_aln.append(chr(ref1_seq[idx_ref-1]))
                final_ref2_aln.append(" ")
                idx_ref -= 1
            elif pointer1[idx_ref,idx_read] == pointer_jump:
                final_read_aln.append(chr(read_seq[idx_read-1]))
                final_ref1_aln.append(chr(ref1_seq[idx_ref-1]))
                final_ref2_aln.append(" ")
                idx_read -= 1
                curr_matrix = 2
                breakpoints_ref1.append(idx_ref-1)
                idx_ref = colmaxesInd2[idx_read]
                breakpoints_read.append(idx_read)
                breakpoints_ref2.append(idx_ref)
                read_path.append(curr_matrix)

        elif curr_matrix == 2:
#            print('     pointer2: ' + str(pointer2[idx_ref,idx_read]))
            if pointer2[idx_ref,idx_read] == pointer_match:
                final_read_aln.append(chr(read_seq[idx_read-1]))
                final_ref1_aln.append(" ")
                final_ref2_aln.append(chr(ref2_seq[idx_ref-1]))
                idx_read -= 1
                idx_ref -= 1
            elif pointer2[idx_ref,idx_read] == pointer_gap_read:
                final_read_aln.append(chr(read_seq[idx_read-1]))
                final_ref1_aln.append(" ")
                final_ref2_aln.append("-")
                idx_read -= 1
            elif pointer2[idx_ref,idx_read] == pointer_gap_ref:
                final_read_aln.append("-")
                final_ref1_aln.append(" ")
                final_ref2_aln.append(chr(ref2_seq[idx_ref-1]))
                idx_ref -= 1
            elif pointer2[idx_ref,idx_read] == pointer_jump:
                final_read_aln.append(chr(read_seq[idx_read-1]))
                final_ref2_aln.append(chr(ref2_seq[idx_ref-1]))
                final_ref1_aln.append(" ")
                idx_read -= 1
                curr_matrix = 1
                breakpoints_ref2.append(idx_ref-1)
                idx_ref = colmaxesInd1[idx_read]
                breakpoints_read.append(idx_read)
                breakpoints_ref1.append(idx_ref)
                read_path.append(curr_matrix)

#                final_read_aln.append(read_seq[idx_read-1])
#                final_ref1_aln.append(" ")
#                final_ref2_aln.append(ref2_seq[idx_ref-1])
#                idx_read -= 1
#                curr_matrix = 1
#                idx_ref = colmaxesInd1[idx_read]
#                breakpoints_ref2.append(idx_ref-1)
#                breakpoints_read.append(idx_read)
#                breakpoints_ref2.append(idx_ref)

    # Reverse the strings.
    final_read_aln_str = ''.join(final_read_aln)[::-1]
    final_ref1_aln_str = ''.join(final_ref1_aln)[::-1]
    final_ref2_aln_str = ''.join(final_ref2_aln)[::-1]
    final_read_path = read_path[::-1]
    final_breakpoints_read = breakpoints_read[::-1]
    final_breakpoints_ref1 = breakpoints_ref1[::-1]
    final_breakpoints_ref2 = breakpoints_ref2[::-1]

    return({
        "read_aln":final_read_aln_str,
        "ref1_aln":final_ref1_aln_str,
        "ref2_aln":final_ref2_aln_str,
        "breakpoints_read":final_breakpoints_read,
        "breakpoints_ref1":final_breakpoints_ref1,
        "breakpoints_ref2":final_breakpoints_ref2,
        "aln_score":max_score, 
        "read_path":final_read_path
        })

