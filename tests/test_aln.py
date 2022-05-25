from ChromBridGE.ChromBridGE_aln import nw_breakpoint

if __name__ == "__main__":
    print('Performing tests..')

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=None,
                ref2_cut_pos=1)
    print('read: ' + read)
    print('ref1: ' + ref1)
    print('ref2: ' + ref2)
    if read !=  'AGGA' or \
        ref1 != 'A   ' or \
        ref2 != ' GGA' or \
        tx_info[0] != False:
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=1,
                ref2_cut_pos=None)
    if read !=  'AGGA' or \
        ref1 != 'A   ' or \
        ref2 != ' GGA':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=2,
                ref2_cut_pos=None)
    if read !=  'AGGA' or \
        ref1 != 'AG  ' or \
        ref2 != '  GA':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=None,
                ref2_cut_pos=2)
    if read !=  'AGGA' or \
        ref1 != 'AG  ' or \
        ref2 != '  GA':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=None)
    if read !=  'AGGA' or \
        ref1 != 'AGG ' or \
        ref2 != '   A':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=None,
                ref2_cut_pos=3)
    if read !=  'AGGA' or \
        ref1 != 'AGG ' or \
        ref2 != '   A':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    #here, the jump-match-mismatch (-6 +3 -3 = -6) is better than the jump at the breakpoint gap-gap-jump (-3+1 -3+1 -6+2 = -8)
    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGTGGA',
                'AGCGG'.replace(" ",""),
                ' GTCGA'.replace(" ",""),
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-6,
                cut_pos_jump_incentive_score=2,
                ref1_cut_pos=2,
                ref2_cut_pos=3,
                cut_pos_insertion_incentive_score=1)
    if read !=  'AGTGGA' or \
        ref1 != 'AG    ' or \
        ref2 != '  TCGA':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    #but here, the insertion incentive makes the jump-match-mismatch (-6 +3 -3 = -6) worse than the jump at the breakpoint gap-gap-jump (-3+3 -3+3 -6+2 = -4)
    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGTGGA',
                'AGCGG'.replace(" ",""),
                ' GTCGA'.replace(" ",""),
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-6,
                cut_pos_jump_incentive_score=2,
                ref1_cut_pos=2,
                ref2_cut_pos=3,
                cut_pos_insertion_incentive_score=3)
    if read !=  'AGTGGA' or \
        ref1 != 'AG--  ' or \
        ref2 != '    GA':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    #this test tests for breakpoint positions as well as preference for gaps to be at the end and the sequences compressed inward
    # e.g. not this: --ACG-TT--
    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                  'ACGTT',
                'GGACAA',
                   'CGTTTAA',
          gap_score=-3,
          mismatch_score=-3,
          jump_score=-2,
          ref1_cut_pos=None,
          ref2_cut_pos=None)
    if read !=  '--ACGTT---' or \
        ref1 != 'GGA       ' or \
        ref2 != '   CGTTTAA' or \
        breakpoints != ([1], [3], [0]) or \
        read_path != [1,2]:
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2 + '\nbreakpoints: ' + str(breakpoints))

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGGGGA',
                'AGGGGG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=1,
                ref2_cut_pos=0)
    if read !=  'AGGGGGA' or \
        ref1 != 'A      ' or \
        ref2 != ' GGGGGA' or \
        tx_info[0] != True:
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGGGGA',
                'AGGGGG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=2)
    if read !=  'AGGGGGA' or \
        ref1 != 'AGG    ' or \
        ref2 != '   GGGA' or \
        tx_info[0] != True:
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGGTGA',
                'AGGGTG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=2,
                cut_pos_tolerance_for_tx=0
                )
    if read !=  'AGGGTGA' or \
        ref1 != 'AGGGT  ' or \
        ref2 != '     GA' or \
        tx_info[0] != False:
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AGGGTGA',
                'AGGGTG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=2,
                cut_pos_tolerance_for_tx=2)
    if read !=  'AGGGTGA' or \
        ref1 != 'AGGGT  ' or \
        ref2 != '     GA' or \
        tx_info[0] != True:
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2 + '\ntxInfo: ' + str(tx_info))

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AAATTT',
                'AAAAAA',
                'TTTTTT',
		gap_score=-3,
		mismatch_score=-3,
		jump_score=-2)
#    print('read: ' + read)
#    print('ref1: ' + ref1)
#    print('ref2: ' + ref2)
    if read !=  'AAATTT' or \
        ref1 != 'AAA   ' or \
        ref2 != '   TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'TTTAAA',
                'AAAAAA',
                'TTTTTT',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'TTTAAA' or \
        ref1 != '   AAA' or \
        ref2 != 'TTT   ' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AAAATTT',
                'AAAAAA',
                'TTTTTTT',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'AAAATTT' or \
        ref1 != 'AAAA   ' or \
        ref2 != '    TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AAAATTT',
                'AAAAA',
                'TTTTTTT',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'AAAATTT' or \
        ref1 != 'AAAA   ' or \
        ref2 != '    TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'AACCTTGG',
                'AAA',
                'CCCTGG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'AACCTTGG' or \
        ref1 != 'AA      ' or \
        ref2 != '  CCCTGG':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'ACTGACTGACTG',
                'ACTGACTGACTG',
                'CCCTGG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'ACTGACTGACTG' or \
        ref1 != 'ACTGACTGACTG' or \
        ref2 != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'ACTGACTGACTG',
                'ACTGCTGACTG',
                'CCCTGG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'ACTGACTGACTG' or \
        ref1 != 'ACTG-CTGACTG' or \
        ref2 != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    #the last character must be a jump (-4 + 3) or a gap (-2)
    #in this case, the perimiter gap extension score doesn't apply here because it's the corner value
    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'ACTGACTGACTG',
                  'CTGCTGACT',
                      'CCCTGG',
	        gap_score=-2,
		mismatch_score=-1,
		jump_score=-4,
                      )
    if read !=  'ACTGACTGACTG' or \
        ref1 != '-CTG-CTGACT ' or \
        ref2 != '           G':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    #the last character can either be a jump (-5 + 3) or a gap (-2)
    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'ACTGACTGACTG',
                  'CTGCTGACT',
                      'CCCTGG',
	        gap_score=-2,
		mismatch_score=-1,
		jump_score=-5,
                      )
    if read !=  'ACTGACTGACTG' or \
        ref1 != '-CTG-CTGACT-' or \
        ref2 != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)

    read,ref1,ref2,breakpoints,score,read_path,tx_info = nw_breakpoint(
                'GTCGACTGACTG',
                'GTCAGTCAGTCA',
                'ACTGACTGACTG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if read !=  'GTCGACTGACTG' or \
        ref1 != 'GTC         ' or \
        ref2 != '   GACTGACTG':
            raise Exception('TEST DID NOT PASS\nread: ' + read + '\nref1: ' + ref1 + '\nref2: ' + ref2)



    print("Tests passed")
