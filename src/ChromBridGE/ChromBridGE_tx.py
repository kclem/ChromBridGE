import numpy as np
def get_first_matching_pos(read_aln_str, ref_aln_str,
        aln_idx_break, read_idx_break, ref_idx_break,
        num_bases_to_check, increment):
    """
    Identify bases that align exactly and return the coordinate where {num_bases_to_check} bases align exactly

    More precisely, given an alignment between a read and ref sequence (read_aln_str, ref_aln_str) start at the alignment position corresponding to the break/jump point in the read (aln_idx_break) and make sure that the last (num_bases_to_check) bases match. If they do: return True, if they don't: increment the index, and try again.
    If at any point a space is encountered (meaning that the alignment has jumped to the other reference, return False

    params:
        read_aln_str: Alignment of the read (including gaps)
        ref_aln_str: Alignment of the ref (including gaps)
        aln_idx_break: Index in the alignment at which the break occurs

        read_idx_break: location of break in the read, 1-based (0 means the break happens before the first base, 1 means it happens after the first base, etc)
        ref_idx_break: location of break in the ref sequence, 1-based (0 means the break happens before the first base, 1 means it happens after the first base, etc)
        num_bases_to_check: how many bases in a row must match exactly
        increment: which direction to move from the break

    returns:
        dict of (success, final_aln_ind, final_read_ind, final_ref_ind)
            success: boolean, whether a match of num_bases_to_check could be found
            aln_ind: the location of that match in the alignment (or -1 if success is False)
            read_ind: the location of that match in the read (or -1 if success is False)
            ref_ind: the location of that match in the ref (or -1 if success is False)
    """

    stop_idx = 0 + num_bases_to_check - 2
    left_right_adjustment = -1 # if I start at the right and move left, the base I'm looking at will be one base to the left of the cut
    if increment > 0:
        stop_idx = len(read_aln_str) - num_bases_to_check + 1
        left_right_adjustment = 0 #if I start at the left and move right, the base I'm looking at will be one base to the right of the cut

    curr_aln_idx = aln_idx_break + left_right_adjustment
    curr_read_idx = read_idx_break
    curr_ref_idx = ref_idx_break
    while curr_aln_idx != stop_idx:
        has_spaces = False
        bases_match = True
        for offset in range(num_bases_to_check):
            aln_idx_to_check = int(curr_aln_idx + (offset * increment))
            if ref_aln_str[aln_idx_to_check] == ' ':
                has_spaces = True
                break
            if ref_aln_str[aln_idx_to_check] != read_aln_str[aln_idx_to_check]:
                bases_match = False
                break

        if has_spaces:
            return({
                "success":False,
                "aln_ind":-1, 
                "read_ind":-1,
                "ref_ind":-1
                })

        if not has_spaces and bases_match:
            return({
                "success":True,
                "aln_ind":curr_aln_idx - left_right_adjustment,
                "read_ind":curr_read_idx,
                "ref_ind":curr_ref_idx,
                })

        curr_aln_idx += increment
        if ref_aln_str[curr_aln_idx] != '-':
            curr_ref_idx += increment
        if read_aln_str[curr_aln_idx] != '-':
            curr_read_idx += increment

    return({
        "success":False,
        "aln_ind":-1, 
        "read_ind":-1,
        "ref_ind":-1
        })

def get_last_matching_pos(read_aln_str, ref_aln_str,
        aln_idx_break, increment, mismatch_tolerance=0, gap_tolerance=0):
    """
    Identify bases that align exactly and return the last coordinate where bases align exactly from the break up to this point

    More precisely, given an alignment between a read and ref sequence (read_aln_str, ref_aln_str) start at the alignment position corresponding to the break/jump point in the read (aln_idx_break) and count the number of bases that match within the mismatch_tolerance and gap_tolerance. Return the index of the last match.
    If at any point a space is encountered (meaning that the alignment has jumped to the other reference), return the last position

    params:
        read_aln_str: Alignment of the read (including gaps)
        ref_aln_str: Alignment of the ref (including gaps)
        aln_idx_break: Index in the alignment at which the break occurs

        increment: which direction to move from the break
        mismatch_tolerance: int How many mismatches to tolerate before returning false
        gap_tolerance: int How many gaps to tolerate before returning false

    returns: tuple of (match_bases_count, match_read_bases_count, match_ref_bases_count)
        match_bases_count: number of bases that matched
        match_read_bases_count: number of bases that matched in read
        match_ref_bases_count: number of bases that matched in ref
    """


    stop_idx = -1 # continue until 0
    left_right_adjustment = -1 # if I start at the right and move left, the base I'm looking at will be one base to the left of the cut
    if increment > 0:
        stop_idx = len(read_aln_str)
        left_right_adjustment = 0 #if I start at the left and move right, the base I'm looking at will be one base to the right of the cut

    curr_aln_idx = aln_idx_break + left_right_adjustment
    gap_count = 0
    mismatch_count = 0
    match_bases_count = 0
    match_read_bases_count = 0
    match_ref_bases_count = 0
    while curr_aln_idx != stop_idx:
        if ref_aln_str[curr_aln_idx] == ' ':
            return (match_bases_count, match_read_bases_count, match_ref_bases_count)
       
        if ref_aln_str[curr_aln_idx] == '-' or read_aln_str[curr_aln_idx] == '-':
            gap_count += 1
            if gap_count > gap_tolerance:
                return (match_bases_count, match_read_bases_count, match_ref_bases_count)
           
        if ref_aln_str[curr_aln_idx] != read_aln_str[curr_aln_idx]:
            mismatch_count += 1
            if mismatch_count > mismatch_tolerance:
                return (match_bases_count, match_read_bases_count, match_ref_bases_count)
        if ref_aln_str[curr_aln_idx] != '-':

            match_ref_bases_count += 1
        if read_aln_str[curr_aln_idx] != '-':
            match_read_bases_count += 1

        match_bases_count += 1
        curr_aln_idx += increment

    return (match_bases_count, match_read_bases_count, match_ref_bases_count)

def get_ref_cut_pos_in_aln(read_aln_str, ref_aln_str, aln_idx_break, read_idx_break, ref_idx_break, ref_cut_pos):
    """
    Given two alignments, get the index within that alignment and the read of the reference cut position

    params:
        read_aln_str: Alignment of the read (including gaps)
        ref_aln_str: Alignment of the ref (including gaps)
        aln_idx_break: Index in the alignment at which the break occurs
        read_idx_break: location of break in the read, 1-based (0 means the break happens before the first base, 1 means it happens after the first base, etc)
        ref_idx_break: location of break in the ref sequence, 1-based (0 means the break happens before the first base, 1 means it happens after the first base, etc)
        ref_cut_pos: position in ref where the predicted cut site is (user input parameter)

    returns:
        ref_cut_pos_in_aln: index of cut position in the alignment
        ref_cut_pos_in_read: index of cut position in the read
    """
    if ref_cut_pos < ref_idx_break: #breakpoint to right of cut position
        curr_pos_in_ref = ref_idx_break
        curr_pos_in_aln = aln_idx_break
        curr_pos_in_read = read_idx_break
        while ref_cut_pos < curr_pos_in_ref:
            curr_pos_in_aln -= 1
            if ref_aln_str[curr_pos_in_aln] != '-':
                curr_pos_in_ref -= 1
            if read_aln_str[curr_pos_in_read] != '-':
                curr_pos_in_read -= 1
            if ref_aln_str[curr_pos_in_aln] == ' ':
                return None,None
        return (curr_pos_in_aln, curr_pos_in_read)

    elif ref_cut_pos > ref_idx_break: #breakpoint to left of cut position
        curr_pos_in_ref = ref_idx_break
        curr_pos_in_aln = aln_idx_break
        curr_pos_in_read = read_idx_break
        while ref_cut_pos > curr_pos_in_ref:
            curr_pos_in_aln += 1
            if ref_aln_str[curr_pos_in_aln] != '-':
                curr_pos_in_ref += 1
            if read_aln_str[curr_pos_in_read] != '-':
                curr_pos_in_read += 1
            if ref_aln_str[curr_pos_in_aln] == ' ':
                return None,None
        return (curr_pos_in_aln, curr_pos_in_read)
    else:
        return (aln_idx_break, read_idx_break)


def get_tx_offset(read_aln_str, ref_aln_str, aln_idx_break, read_idx_break, ref_idx_break, ref_cut_pos, direction, min_num_bases_beyond_cut=4, min_num_bases_before_cut=4, mismatch_tolerance=0, gap_tolerance=0):
    """
    Get possible offset for translocation (before or after given cut)
    If direction is negative, meaning that the alignment to ref happens before the breakpoint to the left:
        If ref_cut_pos is given, the number of perfectly-matched bases to the right(after) the ref_cut_pos are counted = num_bases_beyond_cut
        The index of the first run of {min_num_bases_before_cut} starting at the ref_aln_breakpoint and moving to the left is calculated = num_bases_before_cut
        A putative breakpoint is identified as:
            if ref_cut_pos is given:
                if num_bases_beyond_cut > min_num_bases_beyond_cut: ref_cut_pos+num_bases_beyond_cut
                otherwise: ref_cut_pos - num_bases_before_cut
            otherwise:
                read_idx_break - num_bases_before_cut

    Here, we distinguish between the cut point (ref_cut_pos) which is the predicted cut point provided by the user, and the breakpoint which is determined by the alignment algorithm

    params:
        read_aln_str: Alignment of the read (including gaps)
        ref_aln_str: Alignment of the ref (including gaps)
        aln_idx_break: Index in the alignment at which the break occurs
        read_idx_break: location of break in the read, 1-based (0 means the break happens before the first base, 1 means it happens after the first base, etc)
        ref_idx_break: location of break in the ref sequence, 1-based (0 means the break happens before the first base, 1 means it happens after the first base, etc)
        ref_cut_pos: position in ref where the predicted cut site is (user input parameter)
        direction: 1 or -1, whether alignment extends to the right (1) or the left(-1) of the breakpoint
        min_num_bases_beyond_cut: Min number of matching bases that must be seen beyond the cut site for a breakpoint to be reported beyond the cut site
        min_num_bases_before_cut: Min number of matching bases that must be seen before the cut site - the position of this first match will be reported
        mismatch_tolerance: int How many mismatches to tolerate beyond the cut
        gap_tolerance: int How many gaps to tolerate beyond the cut

    returns:
        dict of:
            num_bases_beyond_cut: number of reference bases beyond cut (or 0 if ref_cut_pos is not given)
            num_bases_before_cut: number of bases before cut (or 0 if {num_bases_to_check} matching bases couldn't be found)
            num_bases_offset: number of bases before (negative) or after (positive) the cut site
            found_break: boolean, whether breakpoint could be found with less than {num_bases_to_check} beyond the cut, but with more than num_bases_to_check bases after the cut
            breakpoint_aln: location of breakpoint in alignment
            breakpoint_ref: location of breakpoint in reference sequence
    """
    found_break = False
    breakpoint_in_aln = 0
    breakpoint_in_ref = 0
    breakpoint_offset = 0
    num_bases_beyond_cut = -1
    num_bases_before_cut = -1
    if ref_cut_pos is not None: #has breakpoint
        ref_cut_pos_in_aln, ref_cut_pos_in_read = get_ref_cut_pos_in_aln(
                read_aln_str = read_aln_str,
                ref_aln_str = ref_aln_str,
                aln_idx_break = aln_idx_break,
                read_idx_break = read_idx_break,
                ref_idx_break = ref_idx_break,
                ref_cut_pos = ref_cut_pos,
                )
        if ref_cut_pos_in_aln is not None: #is none if breakpoint before ref cut pos
            beyond_match_info = get_last_matching_pos(read_aln_str, ref_aln_str,
                ref_cut_pos_in_aln,-1*direction,
                mismatch_tolerance=mismatch_tolerance, gap_tolerance=gap_tolerance)
            before_match_info = get_first_matching_pos(read_aln_str, ref_aln_str, ref_cut_pos_in_aln, ref_cut_pos_in_read, ref_cut_pos, min_num_bases_before_cut, direction)

            num_bases_beyond_cut = beyond_match_info[2]
            if before_match_info['success']:
                num_bases_before_cut = before_match_info['aln_ind'] - ref_cut_pos_in_aln

                found_break = True
                # break has to have min_num_bases_before_cut matching at the cut site for the num_bases_beyond_cut to apply
                # e.g. with min_num_bases_before_cut=4 and min_num_bases_beyond_cut=4
                #  AAAAAAA|AAAA <ref (with cut at |)
                #  AAAAAAA|AAAA <this read would count 4 num_bases_beyond_cut
                #  AAAAA-A|AAAA <this read would not count 4 
                if num_bases_before_cut == 0 and num_bases_beyond_cut >= min_num_bases_beyond_cut:
                    breakpoint_in_ref = ref_cut_pos + (-1*direction*num_bases_beyond_cut)
                    breakpoint_in_aln = ref_cut_pos_in_aln + (-1*direction*num_bases_beyond_cut)
                else:
                    breakpoint_in_aln = before_match_info['aln_ind']
                    breakpoint_in_ref = before_match_info['ref_ind']
            else: #couldn't find bases before breakpoint
                found_break = False
                breakpoint_in_aln = aln_idx_break
                breakpoint_in_ref = ref_idx_break
        else: #breakpoint was before cut
            before_match_info = get_first_matching_pos(read_aln_str, ref_aln_str, aln_idx_break, read_idx_break, ref_idx_break, min_num_bases_before_cut, direction)
            if before_match_info['success']:
                found_break = True
                num_bases_beyond_cut = 0
                num_bases_before_cut = -1*direction*(ref_cut_pos - before_match_info['ref_ind'])
                breakpoint_in_aln = before_match_info['aln_ind']
                breakpoint_in_ref = before_match_info['ref_ind']
            else: #couldn't find bases before breakpoint
                found_break = False
                breakpoint_in_aln = aln_idx_break
                breakpoint_in_ref = ref_idx_break
            
    else: #does not have given cut_points
        before_match_info = get_first_matching_pos(read_aln_str, ref_aln_str, aln_idx_break, read_idx_break, ref_idx_break, min_num_bases_before_cut, direction)
        if before_match_info['success']:
            breakpoint_in_aln = before_match_info['aln_ind']
            num_bases_before_cut = direction*(breakpoint_in_aln - aln_idx_break)
            breakpoint_in_ref = before_match_info['ref_ind']
            found_break = True
    num_bases_offset = -1*num_bases_before_cut
    if num_bases_before_cut == 0 and num_bases_beyond_cut > 0:
        num_bases_offset = num_bases_beyond_cut


    return({
        "num_bases_beyond_cut":num_bases_beyond_cut, 
        "num_bases_before_cut":num_bases_before_cut, 
        "num_bases_offset":num_bases_offset,
        "found_break":found_break, 
        "breakpoint_in_aln":breakpoint_in_aln, 
        "breakpoint_in_ref":breakpoint_in_ref,
        })

def analyze_tx_alignment(read_aln_str, ref1_aln_str, ref2_aln_str,
        breakpoints_read,breakpoints_ref1,breakpoints_ref2,
        read_path, ref1_cut_pos=None, ref2_cut_pos=None,
        min_num_bases_beyond_cut=4, min_num_bases_before_cut=4,
        mismatch_tolerance=0, gap_tolerance=0):
    """
    Refine the alignment and possible translocation sites by trimming transition sites until they have at least min_num_bases_before_cut matching, (or min_num_bases_beyond_cut beyond the cut site)

    params:
        read_aln_str: Alignment of the read (including gaps)
        ref1_aln_str: Alignment of the ref1 (including gaps or spaces if the sequence doesn't align)
        ref2_aln_str: Alignment of the ref2 (including gaps or spaces if the sequence doesn't align)
        breakpoints_read: positions of breakpoints in read discovered by alignment
        breakpoints_ref1: positions of breakpoints in ref1 at which the optimal alignment switches references
        breakpoints_ref2: positions of breakpoints in ref2 at which the optimal alignment switches references
        read_path: index of ref that the read is aligned to, corresponding to the break points (there will be len(breakpoints)+1 items in read_path)
            e.g. if the read path is [1,2]
            len(breakpoints) is 1, and shows the location at which the alignment switches from ref1 to ref2
        ref1_cut_pos: position in ref1 where the predicted cut site is (user input parameter)
        ref2_cut_pos: position in ref2 where the predicted cut site is (user input parameter)
        min_num_bases_beyond_cut: Min number of matching bases that must be seen beyond the cut site for a breakpoint to be reported beyond the cut site
        min_num_bases_before_cut: Min number of matching bases that must be seen before the cut site - the position of this first match will be reported
        mismatch_tolerance: int How many mismatches to tolerate before returning false
        gap_tolerance: int How many gaps to tolerate before returning false

    returns:
        dict containing keys:
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
        tx_lucky_insertions: sum of the left- and right- distances if they extend beyond the cut. If the cut actually happened, these would be lucky insertions that happened to match the uncut reference sequence. If cuts are not given this value is 0.
    """

    #first, get the coordinates of the alignment with regard to the read
    #read_idx_in_aln: array giving the indices of the read in read_aln_str such that read_idx_in_aln[i] would yield the index of the ith base in the read in read_aln_str
    #at which position in the alignment can I find the read base x?
    #             aln ATT-CG
    # read_ind_in_aln 01245

    read_ind_in_aln = np.zeros(len(read_aln_str),dtype=np.intc)
    aln_idx = 0
    for idx in range(len(read_aln_str)):
        if read_aln_str[idx] != '-':
            read_ind_in_aln[aln_idx] = idx
            aln_idx += 1

    final_path = []
    final_breakpoint_ref1 = None
    final_breakpoint_ref2 = None

    left_tx_info = None
    left_dist = None #how many bp match beyond cut site (pos) or deleted (negative)
    left_bp_is_valid = False

    right_tx_info = None
    right_dist = None
    right_bp_is_valid = False
    if len(read_path) > 0:
        if read_path[0] == 1 and read_path[-1] == 2: # transition ref1 to ref2
            final_path = [1,2]
            left_tx_info = get_tx_offset(
                read_aln_str = read_aln_str,
                ref_aln_str =  ref1_aln_str,
                aln_idx_break = read_ind_in_aln[breakpoints_read[0]],
                read_idx_break = breakpoints_read[0],
                ref_idx_break = breakpoints_ref1[0],
                ref_cut_pos = ref1_cut_pos,
                direction = -1,
                min_num_bases_beyond_cut=min_num_bases_beyond_cut,
                min_num_bases_before_cut=min_num_bases_before_cut,
                mismatch_tolerance=mismatch_tolerance,
                gap_tolerance=gap_tolerance,
                )
            left_dist = left_tx_info['num_bases_offset']
            if ref1_cut_pos is not None:
                left_bp_is_valid = True
                final_breakpoint_ref1 = left_tx_info['breakpoint_in_ref']
            right_tx_info = get_tx_offset(
                read_aln_str = read_aln_str,
                ref_aln_str =  ref2_aln_str,
                aln_idx_break = read_ind_in_aln[breakpoints_read[-1]],
                read_idx_break = breakpoints_read[-1],
                ref_idx_break = breakpoints_ref2[-1],
                ref_cut_pos = ref2_cut_pos,
                direction = 1,
                min_num_bases_beyond_cut=min_num_bases_beyond_cut,
                min_num_bases_before_cut=min_num_bases_before_cut,
                mismatch_tolerance=mismatch_tolerance,
                gap_tolerance=gap_tolerance,
                )
            right_dist = right_tx_info['num_bases_offset']
            if ref2_cut_pos is not None:
                right_bp_is_valid = True
                final_breakpoint_ref2 = right_tx_info['breakpoint_in_ref']

        elif read_path[0] == 2 and read_path[-1] == 1: #transition from ref2 to ref1
            final_path = [2,1]
            left_tx_info = get_tx_offset(
                read_aln_str = read_aln_str,
                ref_aln_str =  ref2_aln_str,
                aln_idx_break = read_ind_in_aln[breakpoints_read[0]],
                read_idx_break = breakpoints_read[0],
                ref_idx_break = breakpoints_ref2[0],
                ref_cut_pos = ref2_cut_pos,
                direction = -1,
                min_num_bases_beyond_cut=min_num_bases_beyond_cut,
                min_num_bases_before_cut=min_num_bases_before_cut,
                mismatch_tolerance=mismatch_tolerance,
                gap_tolerance=gap_tolerance,
                )
            left_dist = left_tx_info['num_bases_offset']
            if ref2_cut_pos is not None:
                left_bp_is_valid = True
                final_breakpoint_ref2 = left_tx_info['breakpoint_in_ref']

            right_tx_info = get_tx_offset(
                read_aln_str = read_aln_str,
                ref_aln_str =  ref1_aln_str,
                aln_idx_break = read_ind_in_aln[breakpoints_read[-1]],
                read_idx_break = breakpoints_read[-1],
                ref_idx_break = breakpoints_ref1[-1],
                ref_cut_pos = ref1_cut_pos,
                direction = 1,
                min_num_bases_beyond_cut=min_num_bases_beyond_cut,
                min_num_bases_before_cut=min_num_bases_before_cut,
                mismatch_tolerance=mismatch_tolerance,
                gap_tolerance=gap_tolerance,
                )
            right_dist = right_tx_info['num_bases_offset']
            if ref1_cut_pos is not None:
                right_bp_is_valid = True
                final_breakpoint_ref1 = right_tx_info['breakpoint_in_ref']

    #add ~ for insertions at translocated regions
    read_aln_arr = list(read_aln_str)
    ref1_aln_arr = list(ref1_aln_str)
    ref2_aln_arr = list(ref2_aln_str)

    if left_tx_info is not None and right_tx_info is not None and \
            left_tx_info['found_break'] and right_tx_info['found_break']:
        for idx in range(left_tx_info['breakpoint_in_aln'], right_tx_info['breakpoint_in_aln']):
            ref1_aln_arr[idx] = '~'
            ref2_aln_arr[idx] = '~'

    #now get rid of positions where there is a gap but the alignment is removed
    final_read_aln_arr = []
    final_ref1_aln_arr = []
    final_ref2_aln_arr = []
    bp_match_ref1 = 0
    bp_match_ref2 = 0
    bp_insertion = 0
    for idx in range(len(read_aln_arr)):
        # these lines delete insertions in the alignments that align to gaps (~) but could make the alignment unstable
#        if read_aln_arr[idx] == '-' and ref1_aln_arr[idx] == '~' and ref2_aln_arr[idx] == '~':
#            pass
#        else:
            final_read_aln_arr.append(read_aln_arr[idx])
            final_ref1_aln_arr.append(ref1_aln_arr[idx])
            final_ref2_aln_arr.append(ref2_aln_arr[idx])

            if read_aln_arr[idx] == ref1_aln_arr[idx]:
                bp_match_ref1 += 1
            if read_aln_arr[idx] == ref2_aln_arr[idx]:
                bp_match_ref2 += 1
            if ref1_aln_arr[idx] == '~':
                bp_insertion += 0


    final_read_str = ''.join(final_read_aln_arr)
    final_ref1_str = ''.join(final_ref1_aln_arr)
    final_ref2_str = ''.join(final_ref2_aln_arr)


    is_tx = False
    tx_status = 'Unknown/breakpoints not given (' + str(ref1_cut_pos) + ' and ' + str(ref2_cut_pos) + ')'
    tx_lucky_insertions = 0

    if (ref1_cut_pos is not None) and (ref2_cut_pos is not None):
        if len(read_path) == 1:
            is_tx = False
            tx_status = 'No breakpoints detected'
        #alignments with multiple breakpoints are probably not tx
        elif len(read_path) > 2 and len(final_path) == 0:
            is_tx = False
            tx_status = 'Multiple breakpoints detected'
        elif not left_bp_is_valid or not right_bp_is_valid:
            is_tx = False
            tx_status = 'Valid breakpoints could not be identified'
        else:
            tx_lucky_insertions = max(left_dist,0) + max(right_dist,0)
            if (final_path[0] == 1 and final_path[1] == 2):
                if final_breakpoint_ref1 <= ref1_cut_pos and final_breakpoint_ref2 >= ref2_cut_pos:
                    is_tx = True
                    tx_status = 'Tx A>B'
                else:
                    tx_status = 'Breakpoints incompatible with given cuts'
            elif (final_path[0] == 2 and final_path[1] == 1):
                if final_breakpoint_ref2 <= ref2_cut_pos and final_breakpoint_ref1 >= ref1_cut_pos:
                    is_tx = True
                    tx_status = 'Tx B>A'
                else:
                    tx_status = 'Breakpoints incompatible with given cuts'

    return({"final_read_str":final_read_str,
            "final_ref1_str":final_ref1_str,
            "final_ref2_str":final_ref2_str,

            "final_path":final_path,

            "final_breakpoint_ref1":final_breakpoint_ref1,
            "final_breakpoint_ref2":final_breakpoint_ref2,

            "bp_match_ref1":bp_match_ref1,
            "bp_match_ref2":bp_match_ref2,
            "bp_insertion":bp_insertion,

            "is_tx":is_tx,
            "tx_status":tx_status,
            "left_dist":left_dist,
            "right_dist":right_dist,
            "tx_lucky_insertions":tx_lucky_insertions
            })

if __name__ == "__main__":
#    from ChromBridGE.ChromBridGE_aln import nw_breakpoint
    from ChromBridGE_aln import nw_breakpoint
    print('Performing tests..')

    print('Testing matching bases on both sides of breakpoint with cut given (left_dist, right_dist, and tx_lucky_insertions)')
    # cut is | breakpoint is ^
    seq_A =    "TGTCGCCCCCCAGCCAGCACGGT^CGCACG|CGTCTGCGCTGGGTGATTTGTA"
    read_seq = "TGTCGCCCCCCAGCCAGCACGGT^GTAACG|GACTAATGTGACATGTCCGGCA"
    seq_B =    "ACCGCGGTTTGTTCCTATGTCAG^GTAACG|GACTAATGTGACATGTCCGGCA"
    cut_A = seq_A.replace('^','').index('|') if '|' in seq_A else None
    cut_B = seq_B.replace('^','').index('|') if '|' in seq_B else None
    aln_info = nw_breakpoint(
                read_seq.replace('|','').replace("^",""),
                seq_A.replace('|','').replace("^",""),
                seq_B.replace('|','').replace("^",""),
                ref1_cut_pos=cut_A,
                ref2_cut_pos=cut_B
    )
    print('aln_info: ' + str(aln_info))
    tx_info = analyze_tx_alignment(
            read_aln_str = aln_info['read_aln'],
            ref1_aln_str = aln_info['ref1_aln'],
            ref2_aln_str = aln_info['ref2_aln'],
            breakpoints_read = aln_info['breakpoints_read'],
            breakpoints_ref1 = aln_info['breakpoints_ref1'],
            breakpoints_ref2 = aln_info['breakpoints_ref2'],
            read_path = aln_info['read_path'],
            ref1_cut_pos = cut_A,
            ref2_cut_pos = cut_B)

    print('tx_info: ' + str(tx_info))
    assert(tx_info['left_dist'] == -6) # left side has 6 deletions
    assert(tx_info['right_dist'] == 6) #  right side has 6 insertions
    assert(tx_info['tx_lucky_insertions'] == 6)

    seq_A =    "TGTCGCCCCCCAGCCAGCACGGTCGCACG|CGTCTG^CGCTGGGTGATTTGTA"
    read_seq = "TGTCGCCCCCCAGCCAGCACGGTCGCACG|CGTCTG^TGTGACATGTCCGGCA"
    seq_B =    "ACCGCGGTTTGTTCCTATGTCAGGTAACG|GACTAA^TGTGACATGTCCGGCA"
    cut_A = seq_A.replace('^','').index('|') if '|' in seq_A else None
    cut_B = seq_B.replace('^','').index('|') if '|' in seq_B else None
    aln_info = nw_breakpoint(
                read_seq.replace('|','').replace("^",""),
                seq_A.replace('|','').replace("^",""),
                seq_B.replace('|','').replace("^",""),
                ref1_cut_pos=cut_A,
                ref2_cut_pos=cut_B
    )
    print('aln_info: ' + str(aln_info))
    #print('read: ' + aln_info['read_aln'])
    #print('ref1: ' + aln_info['ref1_aln'])
    #print('ref2: ' + aln_info['ref2_aln'])
    #print('breakpoints_read: ' + str(aln_info['breakpoints_read']))
    #print('breakpoints_ref1: ' + str(aln_info['breakpoints_ref1']))
    #print('breakpoints_ref2: ' + str(aln_info['breakpoints_ref2']))
    #print('read_path: ' + str(aln_info['read_path']))
    tx_info = analyze_tx_alignment(
            read_aln_str = aln_info['read_aln'],
            ref1_aln_str = aln_info['ref1_aln'],
            ref2_aln_str = aln_info['ref2_aln'],
            breakpoints_read = aln_info['breakpoints_read'],
            breakpoints_ref1 = aln_info['breakpoints_ref1'],
            breakpoints_ref2 = aln_info['breakpoints_ref2'],
            read_path = aln_info['read_path'],
            ref1_cut_pos = cut_A,
            ref2_cut_pos = cut_B)

    print('tx_info: ' + str(tx_info))
    assert(tx_info['left_dist'] == 6) # left side has 6 insertions
    assert(tx_info['right_dist'] == -6) #  right side has 6 deletions
    assert(tx_info['tx_lucky_insertions'] == 6)

    print('Testing transitions B>A')
    seq_A =    "ACCGCGGTTTGTTCCTATGTCAG^GTAACG|GACTAATGTGACATGTCCGGCA"
    read_seq = "TGTCGCCCCCCAGCCAGCACGGT^GTAACG|GACTAATGTGACATGTCCGGCA"
    seq_B =    "TGTCGCCCCCCAGCCAGCACGGT^CGCACG|CGTCTGCGCTGGGTGATTTGTA"
    cut_A = seq_A.replace('^','').index('|') if '|' in seq_A else None
    cut_B = seq_B.replace('^','').index('|') if '|' in seq_B else None
    aln_info = nw_breakpoint(
                read_seq.replace('|','').replace("^",""),
                seq_A.replace('|','').replace("^",""),
                seq_B.replace('|','').replace("^",""),
                ref1_cut_pos=cut_A,
                ref2_cut_pos=cut_B
    )
    print('aln_info: ' + str(aln_info))
    tx_info = analyze_tx_alignment(
            read_aln_str = aln_info['read_aln'],
            ref1_aln_str = aln_info['ref1_aln'],
            ref2_aln_str = aln_info['ref2_aln'],
            breakpoints_read = aln_info['breakpoints_read'],
            breakpoints_ref1 = aln_info['breakpoints_ref1'],
            breakpoints_ref2 = aln_info['breakpoints_ref2'],
            read_path = aln_info['read_path'],
            ref1_cut_pos = cut_A,
            ref2_cut_pos = cut_B)

    print('tx_info: ' + str(tx_info))
    assert(tx_info['left_dist'] == -6) # left side has 6 deletions
    assert(tx_info['right_dist'] == 6) #  right side has 6 insertions
    assert(tx_info['tx_lucky_insertions'] == 6)

    #                                         LLLLLLL < 7 deletions from right
    seq_A =    "ACCGCGGTTTGTTCCTATGTCAGGTAACG|GACTAAG^TGTGACATGTCCGGCA"
    read_seq = "TGTCGCCCCCCAGCCAGCACGGTCGCACG|CGTCT^^^TGTGACATGTCCGGCA"
    seq_B =    "TGTCGCCCCCCAGCCAGCACGGTCGCACG|CGTCT^GCGCTGGGTGATTTGTA"
    #                                         RRRRR < 5 insertions from left
    cut_A = seq_A.replace('^','').index('|') if '|' in seq_A else None
    cut_B = seq_B.replace('^','').index('|') if '|' in seq_B else None
    aln_info = nw_breakpoint(
                read_seq.replace('|','').replace("^",""),
                seq_A.replace('|','').replace("^",""),
                seq_B.replace('|','').replace("^",""),
                ref1_cut_pos=cut_A,
                ref2_cut_pos=cut_B
    )
    print('aln_info: ' + str(aln_info))
    tx_info = analyze_tx_alignment(
            read_aln_str = aln_info['read_aln'],
            ref1_aln_str = aln_info['ref1_aln'],
            ref2_aln_str = aln_info['ref2_aln'],
            breakpoints_read = aln_info['breakpoints_read'],
            breakpoints_ref1 = aln_info['breakpoints_ref1'],
            breakpoints_ref2 = aln_info['breakpoints_ref2'],
            read_path = aln_info['read_path'],
            ref1_cut_pos = cut_A,
            ref2_cut_pos = cut_B)

    print('tx_info: ' + str(tx_info))
    assert(tx_info['left_dist'] == 5) # left side has 5 insertions
    assert(tx_info['right_dist'] == -7) #  right side has 7 deletions
    assert(tx_info['tx_lucky_insertions'] == 5)

    seq_A = "TACTGCGCCAATAAGTTACGGTACTGTCGCCCCCCAGCCAGCACGGGGGCACG|CGTCTGCGCTGGGTGATTTGTACATAGTAGCATTGTTAATCAACTCGGCGACGCAGGAGGAAACGGAGAG"
    seq_B = "CAGGAGACGCGCTCATGGTGATGCGCCGCGGTTTGTTCCTATGTCAGGTAACG|CGTCTGCGCTGGGTGATTTGCTCACATAAAAAACGAGAAGGCCAGTCCCACTAATGTGACATGTCCGGCA"
    read_seq = "TACTGCGCCAATAAGTTACGGTACTGTCGCCCCCCAGCCAGCACGGGG|CGCGTCTGCGCTGGGTGATTTGCTCACATAAAAAACGAGAAGGCCAGTCCCACTAATGTGACATGTCCGGCA"
    cut_A = seq_A.index('|') if '|' in seq_A else None
    cut_B = seq_B.index('|') if '|' in seq_B else None
    aln_info = nw_breakpoint(
                read_seq.replace('|',''),
                seq_A.replace('|',''),
                seq_B.replace('|',''),
                ref1_cut_pos=cut_A,
                ref2_cut_pos=cut_B
    )
#    print('read: ' + aln_info['read_aln'])
#    print('ref1: ' + aln_info['ref1_aln'])
#    print('ref2: ' + aln_info['ref2_aln'])
#    print('breakpoints_read: ' + str(aln_info['breakpoints_read']))
#    print('breakpoints_ref1: ' + str(aln_info['breakpoints_ref1']))
#    print('breakpoints_ref2: ' + str(aln_info['breakpoints_ref2']))
#    print('read_path: ' + str(aln_info['read_path']))
    tx_info = analyze_tx_alignment(
            read_aln_str = aln_info['read_aln'],
            ref1_aln_str = aln_info['ref1_aln'],
            ref2_aln_str = aln_info['ref2_aln'],
            breakpoints_read = aln_info['breakpoints_read'],
            breakpoints_ref1 = aln_info['breakpoints_ref1'],
            breakpoints_ref2 = aln_info['breakpoints_ref2'],
            read_path = aln_info['read_path'],
            ref1_cut_pos = cut_A,
            ref2_cut_pos = cut_B)

#    print('tx_info: ' + str(tx_info))
    if tx_info['final_breakpoint_ref1'] != 48 or \
       tx_info['final_breakpoint_ref2'] != 53 or \
       tx_info['tx_status'] != 'Tx A>B':
           raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'] + '\ntx info: ' + str(tx_info))

    seq_A = "CTGACTGACTCAAACATGTAACTTATAACGTCGTTAACCCTTAGTCCAAGGAA|CGGGAGACTAGACTAAGATCAGCATCGGACCGGGGACAATTTTCTATGGCGGTGCTCTATGTTTTCCAGT"
    seq_B = "CACATGTGTATGAAAATGTATGGTATAAACGTTTGTCTTCACAGCGCACGGAA|CGGGAGACTAGACTAAGATCGCTGCTACCCATGGTCAGCTAACGACAAGAGGGCTCTCTCTGAAGCTTAC"
    read_seq = "CTGACTGACTCAAACATGTAACTTATAACGTCGTTAACCCTTAGTCCAAGGA|AAACGCGGGAGACTAGACTAAGATCGCTGCTACCCATGGTCAGCTAACGACAAGAGGGCTCTCTCTGAAGCTTAC"
    cut_A = seq_A.index('|') if '|' in seq_A else None
    cut_B = seq_B.index('|') if '|' in seq_B else None
    aln_info = nw_breakpoint(
                read_seq.replace('|',''),
                seq_A.replace('|',''),
                seq_B.replace('|',''),
                ref1_cut_pos=cut_A,
                ref2_cut_pos=cut_B
    )
#    print('read: ' + aln_info['read_aln'])
#    print('ref1: ' + aln_info['ref1_aln'])
#    print('ref2: ' + aln_info['ref2_aln'])
#    print('breakpoints_read: ' + str(aln_info['breakpoints_read']))
#    print('breakpoints_ref1: ' + str(aln_info['breakpoints_ref1']))
#    print('breakpoints_ref2: ' + str(aln_info['breakpoints_ref2']))
#    print('read_path: ' + str(aln_info['read_path']))
    tx_info = analyze_tx_alignment(
            read_aln_str = aln_info['read_aln'],
            ref1_aln_str = aln_info['ref1_aln'],
            ref2_aln_str = aln_info['ref2_aln'],
            breakpoints_read = aln_info['breakpoints_read'],
            breakpoints_ref1 = aln_info['breakpoints_ref1'],
            breakpoints_ref2 = aln_info['breakpoints_ref2'],
            read_path = aln_info['read_path'],
            ref1_cut_pos = cut_A,
            ref2_cut_pos = cut_B)

#    print('tx_info: ' + str(tx_info))
    if tx_info['final_breakpoint_ref1'] != 53 or \
       tx_info['final_breakpoint_ref2'] != 53 or \
       tx_info['tx_status'] != 'Tx A>B' or \
       tx_info['final_read_str'] != 'CTGACTGACTCAAACATGTAACTTATAACGTCGTTAACCCTTAGTCCAAGGAAAAC-G--CGGGAGACTAGACTAAGATCGCTGCTACCCATGGTCAGCTAACGACAAGAGGGCTCTCTCTGAAGCTTAC' or \
       tx_info['final_ref1_str'] != 'CTGACTGACTCAAACATGTAACTTATAACGTCGTTAACCCTTAGTCCAAGGAA~~~~~~~                                                                      ' or \
       tx_info['final_ref2_str'] != '                                                     ~~~~~~~CGGGAGACTAGACTAAGATCGCTGCTACCCATGGTCAGCTAACGACAAGAGGGCTCTCTCTGAAGCTTAC':
           raise Exception('TEST DID NOT PASS\nread: ' + aln_info['read_aln'] + '\nref1: ' + aln_info['ref1_aln'] + '\nref2: ' + aln_info['ref2_aln'] + '\ntx info: ' + str(tx_info))

    val = get_first_matching_pos(
            read_aln_str = "AGACCCCCCTTATGGGG",
            ref_aln_str =  "      AGATTATGGGG",
            aln_idx_break = 6,
            read_idx_break = 6,
            ref_idx_break = 10,
            num_bases_to_check = 4,
            increment=1
            )
    print(val)
    assert(val['success'] == True)
    assert(val['aln_ind'] == 9)
    assert(val['read_ind'] == 9)
    assert(val['ref_ind'] == 13)


    #check right stop value
    val = get_first_matching_pos(
        read_aln_str = "AGACCCCCCTTAGGGGG",
        ref_aln_str =  "      AGATTATGGGG",
        aln_idx_break = 6,
        read_idx_break = 6,
        ref_idx_break = 10,
        num_bases_to_check = 4,
        increment=1
    )
    assert(val['success'] == True)
    assert(val['aln_ind'] == 13)
    assert(val['read_ind'] == 13)
    assert(val['ref_ind'] == 17)

    val = get_first_matching_pos(
        read_aln_str = "AGACCCC-CTTAGGGGG",
        ref_aln_str =  "      AGATTATGGGG",
        aln_idx_break = 6,
        read_idx_break = 6,
        ref_idx_break = 10,
        num_bases_to_check = 4,
        increment=1
    )
    print(val)
    assert(val['success'] == True)
    assert(val['aln_ind'] == 13)
    assert(val['read_ind'] == 12)
    assert(val['ref_ind'] == 17)

    val = get_first_matching_pos(
        read_aln_str = "AGACCCCTCTTAGGGGG",
        ref_aln_str =  "     CA-ATTATGGGG",
        aln_idx_break = 5,
        read_idx_break = 5,
        ref_idx_break = 10,
        num_bases_to_check = 4,
        increment=1
    )
    assert(val['success'] == True)
    assert(val['aln_ind'] == 13)
    assert(val['read_ind'] == 13)
    assert(val['ref_ind'] == 17)

    val = get_first_matching_pos(
        read_aln_str = "AG-ACCCCTCTTAGGGGG",
        ref_aln_str =  "      CA-ATTATGGGG",
        aln_idx_break = 6,
        read_idx_break = 5,
        ref_idx_break = 10,
        num_bases_to_check = 4,
        increment=1
    )
    assert(val['success'] == True)
    assert(val['aln_ind'] == 14)
    assert(val['read_ind'] == 13)
    assert(val['ref_ind'] == 17)

    val = get_first_matching_pos(
        read_aln_str = "AAAACCCCTTTTGGGG",
        ref_aln_str =  "AAAACGGG        ",
        aln_idx_break = 8,
        read_idx_break = 8,
        ref_idx_break = 8,
        num_bases_to_check = 4,
        increment=-1
    )
    print(val)
    assert(val['success'] == True)
    assert(val['aln_ind'] == 5)
    assert(val['read_ind'] == 5)
    assert(val['ref_ind'] == 5)

    val = get_first_matching_pos(
        read_aln_str = "AAAACC-CCTTTTGGGG",
        ref_aln_str =  "AAAACGGGG        ",
        aln_idx_break = 9,
        read_idx_break = 8,
        ref_idx_break = 9,
        num_bases_to_check = 4,
        increment=-1
    )
    print(val)
    assert(val['success'] == True)
    assert(val['aln_ind'] == 5)
    assert(val['read_ind'] == 5)
    assert(val['ref_ind'] == 5)

    val = get_first_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "AAAAC-GG        ",
        aln_idx_break = 8,
        read_idx_break = 8,
        ref_idx_break = 7,
        num_bases_to_check = 4,
        increment = -1
    )
    assert(val['success'] == True)
    assert(val['aln_ind'] == 4)
    assert(val['read_ind'] == 4)
    assert(val['ref_ind'] == 4)

    val = get_first_matching_pos(
        read_aln_str = "-AAAAACCCTTTTGGGG",
        ref_aln_str =  "AAAAAC-GG        ",
        aln_idx_break = 9,
        read_idx_break = 8,
        ref_idx_break = 8,
        num_bases_to_check = 4,
        increment = -1
    )
    assert(val['success'] == True)
    assert(val['aln_ind'] == 5)
    assert(val['read_ind'] == 4)
    assert(val['ref_ind'] == 5)

    val = get_first_matching_pos(
        read_aln_str = "-AAAAACCCTTTTGGGG",
        ref_aln_str =  "AAAAA            ",
        aln_idx_break = 5,
        read_idx_break = 4,
        ref_idx_break = 5,
        num_bases_to_check = 4,
        increment = -1
    )
    assert(val['success'] == True)
    assert(val['aln_ind'] == 5)
    assert(val['read_ind'] == 4)
    assert(val['ref_ind'] == 5)

    val = get_first_matching_pos(
        read_aln_str = "-AAAAACCCTTTTGGGG",
        ref_aln_str =  "AAAAA            ",
        aln_idx_break = 5,
        read_idx_break = 4,
        ref_idx_break = 5,
        num_bases_to_check = 5,
        increment = -1
    )
    print(val)
    assert(val['success'] == False)
    assert(val['aln_ind'] == -1)
    assert(val['read_ind'] == -1)
    assert(val['ref_ind'] == -1)



    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "   AACCC        ",
        aln_idx_break = 8,
        increment = -1
    )
    print(val)
    assert(val == (5,5,5))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "   AACCC        ",
        aln_idx_break = 8,
        increment = -1,
        mismatch_tolerance = 0,
        gap_tolerance = 1
    )
    print(val)
    assert(val == (5,5,5))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "   AAGCC        ",
        aln_idx_break = 8,
        increment=-1,
        mismatch_tolerance = 1,
        gap_tolerance = 1
    )
    print(val)
    assert(val == (5,5,5))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "  AAAGGC        ",
        aln_idx_break = 8,
        increment=-1,
        mismatch_tolerance = 1,
        gap_tolerance = 1
    )
    print(val)
    assert(val == (2,2,2))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "  AAA-CC        ",
        aln_idx_break = 8,
        increment=-1,
        mismatch_tolerance = 1,
        gap_tolerance = 1
    )
    print(val)
    assert(val == (6,6,5))

    val = get_last_matching_pos(
        read_aln_str = "AAAAA-CCTTTTGGGG",
        ref_aln_str =  "AAAAACCC        ",
        aln_idx_break = 8,
        increment=-1,
        mismatch_tolerance = 1,
        gap_tolerance = 1
    )
    print(val)
    assert(val == (8,7,8))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "        TTAT    ",
        aln_idx_break = 8,
        increment=1
    )
    print(val)
    assert(val == (2,2,2))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "        TTAT    ",
        aln_idx_break = 8,
        increment=1,
        mismatch_tolerance = 1
    )
    print(val)
    assert(val == (4,4,4))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "        TTATGGGG",
        aln_idx_break = 8,
        increment=1,
        mismatch_tolerance = 1
    )
    print(val)
    assert(val == (8,8,8))

    val = get_last_matching_pos(
        read_aln_str = "AAAAACCCTTTTGGGG",
        ref_aln_str =  "        T ATGGGG",
        aln_idx_break = 8,
        increment=1,
        mismatch_tolerance = 1,
        gap_tolerance = 1
    )
    print(val)
    assert(val == (1,1,1))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCC        ',
            aln_idx_break=8,
            read_idx_break=8,
            ref_idx_break=8,
            ref_cut_pos = 4
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (4,4))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCC        ',
            aln_idx_break=8,
            read_idx_break=8,
            ref_idx_break=16,
            ref_cut_pos = 12
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (4,4))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAA-CCC        ',
            aln_idx_break=8,
            read_idx_break=8,
            ref_idx_break=16,
            ref_cut_pos = 12
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (3,3))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = '-AAACCCCTTTTGGGG',
            ref_aln_str =  'AAAA-CCC        ',
            aln_idx_break=8,
            read_idx_break=7,
            ref_idx_break=16,
            ref_cut_pos = 12
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (3,2))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '        TTTTGGG ',
            aln_idx_break=8,
            read_idx_break=8,
            ref_idx_break=8,
            ref_cut_pos = 10
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (10,10))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '        TTTTGGG ',
            aln_idx_break=8,
            read_idx_break=8,
            ref_idx_break=4,
            ref_cut_pos = 6
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (10,10))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = 'AAAACCCCTT-TGGGG',
            ref_aln_str =  '        T-TTGGGG',
            aln_idx_break=8,
            read_idx_break=8,
            ref_idx_break=4,
            ref_cut_pos = 6
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (11,10))

    ref_cut_pos_in_aln = get_ref_cut_pos_in_aln(
            read_aln_str = '-AAAACCCCTT-TGBBB',
            ref_aln_str =  '         T-TTGBB ',
            aln_idx_break=9,
            read_idx_break=8,
            ref_idx_break=4,
            ref_cut_pos = 8
            )
    print(ref_cut_pos_in_aln)
    assert(ref_cut_pos_in_aln == (14,11))


    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCC        ',
            aln_idx_break = 8,
            read_idx_break = 8,
            ref_idx_break = 8,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '--AACCCC        ',
            aln_idx_break = 8,
            read_idx_break = 8,
            ref_idx_break = 6,
            ref_cut_pos = 6,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 6)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCCT       ',
            aln_idx_break = 9,
            read_idx_break = 9,
            ref_idx_break = 9,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCCC       ',
            aln_idx_break = 9,
            read_idx_break = 9,
            ref_idx_break = 9,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCCC       ',
            aln_idx_break = 9,
            read_idx_break = 9,
            ref_idx_break = 9,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCCCC       ',
            aln_idx_break = 9,
            read_idx_break = 9,
            ref_idx_break = 9,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=1,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 9)
    assert(bp['breakpoint_in_ref'] == 9)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACC          ',
            aln_idx_break = 6,
            read_idx_break = 6,
            ref_idx_break = 6,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 2)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 6)
    assert(bp['breakpoint_in_ref'] == 6)

    bp = get_tx_offset(
            read_aln_str = '-AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAAACC          ',
            aln_idx_break = 7,
            read_idx_break = 6,
            ref_idx_break = 7,
            ref_cut_pos = 9,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 2)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 7)
    assert(bp['breakpoint_in_ref'] == 7)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCT         ',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 2)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 6)
    assert(bp['breakpoint_in_ref'] == 6)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACCT         ',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = None,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == -1)
    assert(bp['num_bases_before_cut'] == 1)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 6)
    assert(bp['breakpoint_in_ref'] == 6)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAACC          ',
            aln_idx_break = 6,
            read_idx_break = 6,
            ref_idx_break = 6,
            ref_cut_pos = None,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == -1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 6)
    assert(bp['breakpoint_in_ref'] == 6)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '-AAACC          ',
            aln_idx_break = 6,
            read_idx_break = 6,
            ref_idx_break = 5,
            ref_cut_pos = None,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == -1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 6)
    assert(bp['breakpoint_in_ref'] == 5)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAGCCT         ',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == -1)
    assert(bp['num_bases_before_cut'] == -1)
    assert(bp['found_break'] == False)
    assert(bp['breakpoint_in_aln'] == 7)
    assert(bp['breakpoint_in_ref'] == 7)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  'AAAATCT         ',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = 8,
            direction = -1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 4)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 4)
    assert(bp['breakpoint_in_ref'] == 4)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '        TTTTGGG-',
            aln_idx_break = 8,
            read_idx_break = 8,
            ref_idx_break = 8,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '        GTTTGGG-',
            aln_idx_break = 8,
            read_idx_break = 8,
            ref_idx_break = 8,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 0)
    assert(bp['num_bases_before_cut'] == 1)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 9)
    assert(bp['breakpoint_in_ref'] == 9)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '       CGTTTGGG-',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 1)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 9)
    assert(bp['breakpoint_in_ref'] == 9)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '       CTTTTGGG-',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=1,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'AAAACCCCTTTTGGGG',
            ref_aln_str =  '       CTTTTGGG-',
            aln_idx_break = 7,
            read_idx_break = 7,
            ref_idx_break = 7,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=1,
            min_num_bases_before_cut=0,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 7)
    assert(bp['breakpoint_in_ref'] == 7)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '    GGGCTTTTGGG-',
            aln_idx_break = 4,
            read_idx_break = 3,
            ref_idx_break = 4,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 8)
    assert(bp['breakpoint_in_ref'] == 8)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '    GGGAGTTTGGG-',
            aln_idx_break = 4,
            read_idx_break = 3,
            ref_idx_break = 4,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 9)
    assert(bp['breakpoint_in_ref'] == 9)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '    GGGAGTTTGGG-',
            aln_idx_break = 4,
            read_idx_break = 3,
            ref_idx_break = 6,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 9)
    assert(bp['breakpoint_in_ref'] == 11)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '    GGGCCCTTGGG-',
            aln_idx_break = 4,
            read_idx_break = 3,
            ref_idx_break = 4,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 1)
    assert(bp['num_bases_before_cut'] == 2)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 10)
    assert(bp['breakpoint_in_ref'] == 10)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '   TCCCCTTTTGGG-',
            aln_idx_break = 3,
            read_idx_break = 2,
            ref_idx_break = 3,
            ref_cut_pos = 8,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == 4)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 4)
    assert(bp['breakpoint_in_ref'] == 4)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '   TCCCCTTTTGGG-',
            aln_idx_break = 3,
            read_idx_break = 2,
            ref_idx_break = 3,
            ref_cut_pos = None,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == -1)
    assert(bp['num_bases_before_cut'] == 1)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 4)
    assert(bp['breakpoint_in_ref'] == 4)

    bp = get_tx_offset(
            read_aln_str = 'A-AACCCCTTTTGGGG',
            ref_aln_str =  '    CCCCTTTTGGG-',
            aln_idx_break = 4,
            read_idx_break = 3,
            ref_idx_break = 5,
            ref_cut_pos = None,
            direction = 1,
            min_num_bases_beyond_cut=4,
            min_num_bases_before_cut=4,
            mismatch_tolerance=0,
            gap_tolerance=0)
    print(bp)
    assert(bp['num_bases_beyond_cut'] == -1)
    assert(bp['num_bases_before_cut'] == 0)
    assert(bp['found_break'] == True)
    assert(bp['breakpoint_in_aln'] == 4)
    assert(bp['breakpoint_in_ref'] == 5)


    print('Finished tests')
