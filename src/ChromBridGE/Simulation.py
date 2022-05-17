import random

nucs = ['A','T','C','G']

def makeRandomGuides(num_mutations_in_guide = 1,predicted_cut_position=3,guide_length_bp=23,debug=False):
    """
    Produces two guide sequences with specified number of mutations between them

    params:
        num_mutations_in_guide: number of bases that differ in the guide beteween sequence A and B
        predicted_cut_position: position for predicted cut (number of bp from left end of guide)
        guide_length_bp: guide length in number of bp
        debug: if true, dummy sequences for debugging are produced

    returns:
        guide_left_A: left part of guide A (before predicted cut site
        guide_right_A: right part of guide A (after predicted cut site
        guide_left_B: left part of guide B (before predicted cut site with mutaitons compared to guide A
        guide_right_B: right part of guide B (after predicted cut site with mutations compared to guide B
    """
    guide_seq = random.choices(nucs,k=guide_length_bp) #array of nucleotides
    guide_seq_A = guide_seq[:]
    guide_seq_B = guide_seq[:]

    # introduce mutations in guide region for sequence B
    locs_to_mutate = random.sample(range(guide_length_bp),guide_length_bp)
    for i in range(num_mutations_in_guide):
        old_nuc = guide_seq_B[locs_to_mutate[i]]
        new_nuc = random.choice(nucs)
        while old_nuc == new_nuc:
            new_nuc = random.choice(nucs)
        guide_seq_B[locs_to_mutate[i]] = new_nuc

    if debug:
        guide_seq_A = ['G']*guide_length_bp
        guide_seq_B = ['H']*guide_length_bp
#    guide_left_A = "".join(guide_seq_A[0:17])
#    guide_right_A = "".join(guide_seq_A[17:guide_length_bp])

    guide_left_A = "".join(guide_seq_A[0:predicted_cut_position])
    guide_right_A = "".join(guide_seq_A[predicted_cut_position:guide_length_bp])

    guide_left_B = "".join(guide_seq_B[0:predicted_cut_position])
    guide_right_B = "".join(guide_seq_B[predicted_cut_position:guide_length_bp])

    return (guide_left_A,guide_right_A,guide_left_B,guide_right_B)

def makeRandomAmplicons(
        guide_left_A,
        guide_right_A,
        guide_left_B,
        guide_right_B,
        left_length_A = 50,
        right_length_A = 50,
        left_length_B = 50,
        right_length_B = 50,
        debug=False):
    """
    Produces two simulated reference sequences surrounding given guide sequences
    sequence A:
    {left_length_A} {guide} {right_length_A}
    sequence B:
    {left_length_B} {guide (with num_mutations_in_guide mutations} {right_length_B}
    read:
    {left from A} {guide A/B} {right from B}

    params:
        guide_left_A: left part of guide A (before predicted cut site)
        guide_right_A: right part of guide A (after predicted cut site)
        guide_left_B: left part of guide B (before predicted cut site with mutaitons compared to guide A)
        guide_right_B: right part of guide B (after predicted cut site with mutations compared to guide B)
        left_length_A: number of bases on the left side of the guide for sequence A
        right_length_A: number of bases on the right side of the guide for sequence A
        left_length_B: number of bases on the left side of the guide for sequence B
        right_length_B: number of bases on the right side of the guide for sequence B
        debug: if true, dummy sequences for debugging are produced

    returns:
        left_A: left part of sequence A including guide (up to cut site)
        right_A: right part of sequence A including guide (starting from cut site)
        left_B: left part of sequence B including guide (up to cut site)
        right_B: right part of sequence B including guide (starting from cut site)
    """

    left_A = "".join(random.choices(nucs,k=left_length_A)) + guide_left_A
    right_A = guide_right_A + "".join(random.choices(nucs,k=right_length_A))

    left_B = "".join(random.choices(nucs,k=left_length_B)) + guide_left_B
    right_B = guide_right_B + "".join(random.choices(nucs,k=right_length_B))

    if debug:
        left_A = "".join(random.choices(['A'],k=left_length_A)) + guide_left_A
        right_A = guide_right_A + "".join(random.choices(['B'],k=right_length_A))

        left_B = "".join(random.choices(['C'],k=left_length_B)) + guide_left_B
        right_B = guide_right_B + "".join(random.choices(['D'],k=right_length_B))

    return (left_A,right_A,left_B,right_B)

def makeSimulatedRead(
        left_A,
        right_A,
        left_B,
        right_B,
        left_length_read = 20,
        right_length_read = 20,
        translocation_position = 0,
        debug=False):
    """
    Creates a read from the given sequence and a translocation position

    params:
        left_A: left part of sequence A including guide (up to cut site)
        right_A: right part of sequence A including guide (starting from cut site)
        left_B: left part of sequence B including guide (up to cut site)
        right_B: right part of sequence B including guide (starting from cut site)
        left_length_read: number of bases to the left of the translocation position for the read
        right_length_read: number of bases to the right of the translocation position for the read
        translocation_position: position relative to the known cut site (0 is the cut site, -1 would 1bp left, 5 would be 5bp right)
        debug: if true, dummy sequences for debugging are produced

    returns:
        read_A_left: the part of the read coming from the left part of seq A
        read_A_right: the part of the read coming from the left part of seq A
        read_B_left: the part of the read coming from the left part of seq B
        read_B_right: the part of the read coming from the left part of seq B

    """

    a_left_start = min(len(left_A),len(left_A)+translocation_position - left_length_read)
    a_left_end = min(len(left_A),len(left_A)+translocation_position)
    read_A_left = left_A[a_left_start:a_left_end]

    a_right_start = max(0,translocation_position - left_length_read)
    a_right_end = max(0,translocation_position)
    read_A_right = right_A[a_right_start:a_right_end]


    b_left_start = min(len(left_B),len(left_B) + translocation_position)
    b_left_end = min(len(left_B),len(left_B) + translocation_position + right_length_read)
    read_B_left = left_B[b_left_start:b_left_end]

    b_right_start = max(0,translocation_position)
    b_right_end = max(0,translocation_position + right_length_read)
    read_B_right = right_B[b_right_start:b_right_end]

    if debug:
        print('makeSimulatedRead:')
        print('\ta_left_start: ' + str(a_left_start))
        print('\ta_left_end: ' + str(a_left_end))
        print('\ta_right_start: ' + str(a_right_start))
        print('\ta_right_end: ' + str(a_right_end))
        print('\tb_left_start: ' + str(b_left_start))
        print('\tb_left_end: ' + str(b_left_end))
        print('\tb_right_start: ' + str(b_right_start))
        print('\tb_right_end: ' + str(b_right_end))


    return(read_A_left,read_A_right,read_B_left,read_B_right)


def simulateOne(num_mutations_in_guide = 1,
        predicted_cut_position=17,
        guide_length_bp=23,
        left_length_read = 20,
        right_length_read = 20,
        left_length_A = 50,
        right_length_A = 50,
        left_length_B = 50,
        right_length_B = 50,
        translocation_position = 0,
        debug = False):
    """
    Produces simulated sequences of the following form:
    sequence A:
    {left_length_A} {guide} {right_length_A}
    sequence B:
    {left_length_B} {guide (with num_mutations_in_guide mutations)} {right_length_B}
    read:
    {left from A} {guide A/B} {right from B}

    params:
        num_mutations_in_guide: number of bases that differ in the guide beteween sequence A and B
        predicted_cut_position: position for predicted cut (number of bp from left end of guide)
        guide_length_bp: guide length in number of bp
        translocation_position: position relative to the known cut site (0 is the cut site, -1 would 1bp left, 5 would be 5bp right)
        left_length_read: number of bases to the left of the translocation position for the read
        right_length_read: number of bases to the right of the translocation position for the read
        left_length_A: number of bases on the left side of the guide for sequence A
        right_length_A: number of bases on the right side of the guide for sequence A
        left_length_B: number of bases on the left side of the guide for sequence B
        right_length_B: number of bases on the right side of the guide for sequence B
        debug: if true, dummy sequences for debugging are produced
    """
    (guide_left_A,guide_right_A,guide_left_B,guide_right_B) = makeRandomGuides(num_mutations_in_guide=num_mutations_in_guide, predicted_cut_position=predicted_cut_position,guide_length_bp=guide_length_bp,debug=debug)

    (left_A,right_A,left_B,right_B) = makeRandomAmplicons(guide_left_A,guide_right_A,guide_left_B,guide_right_B,left_length_A,right_length_A,left_length_B,right_length_B,debug)

    (read_A_left,read_A_right,read_B_left,read_B_right) = makeSimulatedRead(left_A,right_A,left_B,right_B,left_length_read,right_length_read,translocation_position,debug)

    if debug:
        print('A:'+left_A+"|"+right_A)
        print('B:'+left_B+"|"+right_B)
        print('left read: ' + read_A_left+"|"+read_A_right)
        print('right read: ' + read_B_left+"|"+read_B_right)

    return(read_A_left+read_A_right+read_B_left+read_B_right, left_A+right_A,left_B+right_B)


if __name__ == "__main__":
    print('Performing tests..')

    read,refA,refB = simulateOne(num_mutations_in_guide = 1,
            left_length_read = 25,
            right_length_read = 15,
            left_length_A = 40,
            right_length_A = 60,
            left_length_B = 20,
            right_length_B = 30,
            translocation_position = -8,
            debug = True)
    print('refA: ' + str(refA))
    print('refB: ' + str(refB))
    print('read: ' + str(read))
    if read != "AAAAAAAAAAAAAAAAGGGGGGGGGHHHHHHHHHHHHHHD":
        raise Exception('TEST DID NOT PASS')

    read,refA,refB = simulateOne(num_mutations_in_guide = 1,
            left_length_read = 20,
            right_length_read = 20,
            left_length_A = 50,
            right_length_A = 50,
            left_length_B = 50,
            right_length_B = 50,
            translocation_position = -8,
            debug = True)
    if read != "AAAAAAAAAAAGGGGGGGGGHHHHHHHHHHHHHHDDDDDD":
        raise Exception('TEST DID NOT PASS')

    read,refA,refB = simulateOne(num_mutations_in_guide = 1,
            left_length_read = 20,
            right_length_read = 20,
            left_length_A = 50,
            right_length_A = 50,
            left_length_B = 50,
            right_length_B = 50,
            translocation_position = -3,
            debug = True)
    if read != "AAAAAAGGGGGGGGGGGGGGHHHHHHHHHDDDDDDDDDDD":
        raise Exception('TEST DID NOT PASS')

    read,refA,refB = simulateOne(num_mutations_in_guide = 1,
            left_length_read = 20,
            right_length_read = 20,
            left_length_A = 50,
            right_length_A = 50,
            left_length_B = 50,
            right_length_B = 50,
            translocation_position = 0,
            debug = True)
    if read != "AAAGGGGGGGGGGGGGGGGGHHHHHHDDDDDDDDDDDDDD":
        raise Exception('TEST DID NOT PASS')

    read,refA,refB = simulateOne(num_mutations_in_guide = 1,
            left_length_read = 20,
            right_length_read = 20,
            left_length_A = 50,
            right_length_A = 50,
            left_length_B = 50,
            right_length_B = 50,
            translocation_position = 3,
            debug = True)
    if read != "GGGGGGGGGGGGGGGGGGGGHHHDDDDDDDDDDDDDDDDD":
        raise Exception('TEST DID NOT PASS')

    read,refA,refB = simulateOne(num_mutations_in_guide = 1,
            left_length_read = 20,
            right_length_read = 20,
            left_length_A = 50,
            right_length_A = 50,
            left_length_B = 50,
            right_length_B = 50,
            translocation_position = 8,
            debug = True)
    if read != "GGGGGGGGGGGGGGGGGGBBDDDDDDDDDDDDDDDDDDDD":
        raise Exception('TEST DID NOT PASS')

    print("Tests passed")
    


