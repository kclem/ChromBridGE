from ChromBridGE.ChromBridGE_aln import nw_breakpoint

if __name__ == "__main__":
    print('Performing tests..')

    aln_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=None,
                ref2_cut_pos=1)
    print('read: ' + aln_info['read_aln'])
    print('ref1: ' + aln_info['ref1_aln'])
    print('ref2: ' + aln_info['ref2_aln'])
    if aln_info['read_aln'] != 'AGGA' or \
       aln_info['ref1_aln'] != 'A   ' or \
       aln_info['ref2_aln'] != ' GGA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=1,
                ref2_cut_pos=None)
    if aln_info['read_aln'] != 'AGGA' or \
       aln_info['ref1_aln'] != 'A   ' or \
       aln_info['ref2_aln'] != ' GGA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=2,
                ref2_cut_pos=None)
    if aln_info['read_aln'] != 'AGGA' or \
       aln_info['ref1_aln'] != 'AG  ' or \
       aln_info['ref2_aln'] != '  GA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=None,
                ref2_cut_pos=2)
    if aln_info['read_aln'] != 'AGGA' or \
       aln_info['ref1_aln'] != 'AG  ' or \
       aln_info['ref2_aln'] != '  GA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=None)
    if aln_info['read_aln'] != 'AGGA' or \
       aln_info['ref1_aln'] != 'AGG ' or \
       aln_info['ref2_aln'] != '   A':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGA',
                'AGGG',
                'GGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=None,
                ref2_cut_pos=3)
    if aln_info['read_aln'] != 'AGGA' or \
       aln_info['ref1_aln'] != 'AGG ' or \
       aln_info['ref2_aln'] != '   A':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGTGGA',
                'AGCGG'.replace(" ",""),
                ' GTCGA'.replace(" ",""),
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-6,
                cut_pos_jump_incentive_score=2,
                ref1_cut_pos=2,
                ref2_cut_pos=2,
                )
    if aln_info['read_aln'] != 'AGTGGA' or \
       aln_info['ref1_aln'] != 'AG    ' or \
       aln_info['ref2_aln'] != '  TCGA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    #this test tests for breakpoint positions as well as preference for gaps to be at the end and the sequences compressed inward
    # e.g. not this: --ACG-TT--
    #this has multiple solutions...
    aln_info = nw_breakpoint(
                  'ACGTT',
                'GGACAA',
                   'CGTTTAA',
          gap_score=-3,
          mismatch_score=-3,
          jump_score=-2,
          ref1_cut_pos=None,
          ref2_cut_pos=None)
#    print('aln info: ' + str(aln_info))
    if aln_info['read_aln'] != '--ACGTT---' or \
       aln_info['ref1_aln'] != 'GGAC      ' or \
       aln_info['ref2_aln'] != '    GTTTAA' or \
       aln_info['breakpoints_read'] != [2] or \
       aln_info['breakpoints_ref1'] != [4] or \
       aln_info['breakpoints_ref2'] != [1] or \
       aln_info['read_path'] != [1,2]:
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGGGGA',
                'AGGGGG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=1,
                ref2_cut_pos=0)
    if aln_info['read_aln'] != 'AGGGGGA' or \
       aln_info['ref1_aln'] != 'A      ' or \
       aln_info['ref2_aln'] != ' GGGGGA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGGGGA',
                'AGGGGG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=2)
    if aln_info['read_aln'] != 'AGGGGGA' or \
       aln_info['ref1_aln'] != 'AGG    ' or \
       aln_info['ref2_aln'] != '   GGGA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGGTGA',
                'AGGGTG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=2,
                )
    if aln_info['read_aln'] != 'AGGGTGA' or \
       aln_info['ref1_aln'] != 'AGGGT  ' or \
       aln_info['ref2_aln'] != '     GA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AGGGTGA',
                'AGGGTG',
                'GGGGGA',
                gap_score=-3,
                mismatch_score=-3,
                jump_score=-2,
                cut_pos_jump_incentive_score=1,
                ref1_cut_pos=3,
                ref2_cut_pos=2,
                )
    if aln_info['read_aln'] != 'AGGGTGA' or \
       aln_info['ref1_aln'] != 'AGGGT  ' or \
       aln_info['ref2_aln'] != '     GA':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AAATTT',
                'AAAAAA',
                'TTTTTT',
		gap_score=-3,
		mismatch_score=-3,
		jump_score=-2)
#    print('read: ' + read)
#    print('ref1: ' + ref1)
#    print('ref2: ' + ref2)
    if aln_info['read_aln'] != 'AAATTT' or \
       aln_info['ref1_aln'] != 'AAA   ' or \
       aln_info['ref2_aln'] != '   TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'TTTAAA',
                'AAAAAA',
                'TTTTTT',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'TTTAAA' or \
       aln_info['ref1_aln'] != '   AAA' or \
       aln_info['ref2_aln'] != 'TTT   ' :
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AAAATTT',
                'AAAAAA',
                'TTTTTTT',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'AAAATTT' or \
       aln_info['ref1_aln'] != 'AAAA   ' or \
       aln_info['ref2_aln'] != '    TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AAAATTT',
                'AAAAA',
                'TTTTTTT',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'AAAATTT' or \
       aln_info['ref1_aln'] != 'AAAA   ' or \
       aln_info['ref2_aln'] != '    TTT' :
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'AACCTTGG',
                'AAA',
                'CCCTGG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'AACCTTGG' or \
       aln_info['ref1_aln'] != 'AA      ' or \
       aln_info['ref2_aln'] != '  CCCTGG':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'ACTGACTGACTG',
                'ACTGACTGACTG',
                'CCCTGG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'ACTGACTGACTG' or \
       aln_info['ref1_aln'] != 'ACTGACTGACTG' or \
       aln_info['ref2_aln'] != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'ACTGACTGACTG',
                'ACTGCTGACTG',
                'CCCTGG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'ACTGACTGACTG' or \
       aln_info['ref1_aln'] != 'ACTG-CTGACTG' or \
       aln_info['ref2_aln'] != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    #the last character must be a jump (-4 + 3) or a gap (-2)
    #in this case, the perimiter gap extension score doesn't apply here because it's the corner value
    aln_info = nw_breakpoint(
                'ACTGACTGACTG',
                  'CTGCTGACT',
                      'CCCTGG',
	        gap_score=-2,
		mismatch_score=-1,
		jump_score=-4,
                      )
    if aln_info['read_aln'] != 'ACTGACTGACTG' or \
       aln_info['ref1_aln'] != '-CTG-CTGACT ' or \
       aln_info['ref2_aln'] != '           G':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    #the last character can either be a jump (-5 + 3) or a gap (-2)
    aln_info = nw_breakpoint(
                'ACTGACTGACTG',
                  'CTGCTGACT',
                      'CCCTGG',
	        gap_score=-2,
		mismatch_score=-1,
		jump_score=-5,
                      )
    if aln_info['read_aln'] != 'ACTGACTGACTG' or \
       aln_info['ref1_aln'] != '-CTG-CTGACT-' or \
       aln_info['ref2_aln'] != '            ':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])

    aln_info = nw_breakpoint(
                'GTCGACTGACTG',
                'GTCAGTCAGTCA',
                'ACTGACTGACTG',
		gap_score=-2,
		mismatch_score=-1,
		jump_score=-3)
    if aln_info['read_aln'] != 'GTCGACTGACTG' or \
       aln_info['ref1_aln'] != 'GTC         ' or \
       aln_info['ref2_aln'] != '   GACTGACTG':
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'])



    #test for shifting jump to the left or right when not at cut site
    aln_info = nw_breakpoint(
                'AAATGGG'.replace(' ',''),
                'AAATG'.replace(' ',''),
                '  ATGGG'.replace(' ',''),
          gap_score=-3,
          mismatch_score=-3,
          jump_score=-2,
          ref1_cut_pos=0,
          ref2_cut_pos=0)
    if aln_info['read_aln'] != 'AAATGGG' or \
       aln_info['ref1_aln'] != 'AA     ' or \
       aln_info['ref2_aln'] != '  ATGGG' or \
       aln_info['breakpoints_read'] != [2] or \
       aln_info['breakpoints_ref1'] != [2] or \
       aln_info['breakpoints_ref2'] != [0] or \
       aln_info['read_path'] != [1,2]:
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: "' + aln_info['ref1_aln'] + '"\nref2: "' + aln_info['ref2_aln']+'"')

    aln_info = nw_breakpoint(
                'AAATGGG'.replace(' ',''),
                'AAATG'.replace(' ',''),
                '  ATGGG'.replace(' ',''),
          gap_score=-3,
          mismatch_score=-3,
          jump_score=-2,
          ref1_cut_pos=5,
          ref2_cut_pos=5)
    if aln_info['read_aln'] != 'AAATGGG' or \
       aln_info['ref1_aln'] != 'AAATG  ' or \
       aln_info['ref2_aln'] != '     GG' or \
       aln_info['breakpoints_read'] != [5] or \
       aln_info['breakpoints_ref1'] != [5] or \
       aln_info['breakpoints_ref2'] != [3] or \
       aln_info['read_path'] != [1,2]:
            raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: "' + aln_info['ref1_aln'] + '"\nref2: "' + aln_info['ref2_aln']+'"')



    print("Tests passed")
