/** 
 * mismatch_and_alignment
 *
 * This function performs two operations:
 * - Computes the Levenshtein distance between a specified substring of s1,
 *    and s2, up to a maximum bound on edit operations
 * - Determines the best alignment of s2 to the full s1, also bounded by a
 *    different number of edit operations.
 *
 * The Levenshtein distance is computed up to length barcode_length characters
 * of s2. If there is a match within barcode_length, continue processing for
 * all of the length of s1 and determine the aligned length of s2 (this is 
 * used to know how much to trim).
 *
 * If there are more characters left in s2 after all characters in s1 are used,
 * there is no penalty for these extra characters.
 *
 * This algorithm is faster than a naive approach, as it applies the max_val
 * bounds at every step, excluding a lot of computations which are guaranteed
 * to yield a distance above or equal to max_val. It aborts as soon as there
 * is no chance of getting a value below max_val.
 *
 * Parameters:
 *  max_val_barcode: The max value to return for Levenshtein distance of 
 *                   barcode (all higher values are clipped). Such a
 *                   threshold is relevant only for performance reasons.
 *  max_val_align:   Limit of edit operations allowed for alignment of the
 *                   whole sequence of s1. This parameter works similarly
 *                   to max_val_barcode, but the distance of the full string
 *                   alignment is not actually returned. Instead, if alignment
 *                   reaches this value, the aligned distance of zero is
 *                   returned.
 *  barcode_length:  Length of substring of s1 to use for barcode comparison.
 *                   The bound max_val_barcode is only applied to this sequence.
 *  s1:              The first string, representing the known index and spacer
 *                   sequence, concatenated.
 *  s2:              The data to examine for alignment and mismatches. The full
 *                   length read can be passed without performance impact.
 *                   
 *
 * Returns:
 *  first:  The edit distance between s2 and substring(s1, 0, barcode_length),
 *          using the best matching substring of s2 (starting at 0). If the
 *          distance exceeds max_val_barcode, that value is returned instead.
 *  second: The length of s2 when aligned to the full sequence of s1. The length
 *          of the match with fewest edit operations is returned. If several
 *          matches have equal distance, the shortest match is returned (please
 *          see loop before return statement, this could change). If the edit
 *          distance of the barcode is greater or equal to max_val_barcode, the
 *          value of second is undefined. If the total number of edit operations is
 *          equal or greater than max_val_align, the value zero is returned.
 *
 */

pair<int,int> mismatch_and_alignment(int max_val_barcode, int max_val_align,
        int barcode_length, const string& s1, const string& s2)
{
    size_t s1len = s1.length();
    size_t s2len = min(s2.length(), s1len + (size_t)max_val_align);

    // Buffer for two columns of mismatch value results
    int buf1[s2len+1], buf2[s2len+1];

    int* column = buf1, * prevcol = buf2;

    // Columns include a dummy cell at index 0. Every time an element in
    // the column is indexed, one has to add 1 to the index i.

    // Initialise "prev" column before the actual data. This encodes the 
    // penalties for insertions at the beginning of the aligned string.
    // As a special case, prevcol[0] is used in case there is a match at
    // the first position of s1 and s2.
    for (int j=0; j<s2len+1; j++) {
        prevcol[j] = j;
    }
    
    // - first points to the first row to analyse. It is greater than 0 when
    // rows are excluded by the max_val bound. first is the actual position
    // in s2, not the index into the column array (which is first+1).
    // - last is the index of the last row to analyse (position in s2). Points
    // to the last element in column with value below max_val. The initial value
    // represents max_val-1 insertions before s2. Each column opens up another
    // row for last.
    int first = 0, last = min(max_val_align-1, (int)s2len-1);
    int index_mismatches = index_mismatches;
    for (int j = 0; j < s1len; j++) { // Column!

        int i; // i is an index of the row / position in s2, does not include the
               // offset of 1, as in the column arrays.

        int minval = max_val_barcode;
        
        if (first > 0) {
            // If first is nonzero, some rows have already been deemed excluded
            // Row first here correspond to position first - 1 (due to one-based
            // offset in the array)
            column[first] = max_val_align;
        }
        else {
            // The first element of column represents deletions at the start of the
            // aligned string. Not an actual position in s2.
            column[0] = j+1;
        }

        int last_below = -1, first_tmp = first;
        first = s2len; // unless found, set to position past end of list

        // Second loop: find new values in range; keep track of last row below
        // threshold, so we can possibly constrain the "last"
        for (i=first_tmp; i <= last && i < s2len; i++) {
            auto possibilities = {
                prevcol[1 + i] + 1,
                column[i] + 1,
                prevcol[i] + (s1[j] == s2[i] ? 0 : 1)
            };
            int best_possibility = std::min(possibilities);
            column[i+1] = best_possibility;
            if (best_possibility < max_val_align) {
                if (first == s2len) first = i;
                last_below = i;
                minval = min(best_possibility, minval);
            }
        }
        // Possibly, we allow to go down one diagonally in case we're still within
        // the threshold
        if (i == last + 1 && i < s2len) {
            auto possibilities = {
                column[i] + 1,
                prevcol[i] + (s1[j] == s2[i] ? 0 : 1)
            };
            int best_possibility = std::min(possibilities);
            column[i + 1] = best_possibility;
            if (best_possibility < max_val_align) {
                if (first == s2len) first = i;
                last_below = i;
                minval = min(best_possibility, minval);
            }
        }

        // Next round 
        last = last_below;

        /*
        cerr << " j is " << j << " | max barcode= " << max_val_barcode
            << " | max align= " << max_val_align << endl;
        for (int k=-1; k<(int)s2len; ++k) {
           cerr << "prevcol[" << k << " + 1] = " << prevcol[k + 1];
           cerr << " \tcolumn[" << k << " + 1]  = " << column[k + 1];
           if (k == first) cerr << " <- first";
           if (k == last) cerr << " <- last";
           cerr << endl;
        }
        */
        
        if (j < barcode_length && minval >= max_val_barcode) {
            //cout << "RET: too high minval" << endl;
            return make_pair(max_val_barcode, 0);
        }
        
        if (j == barcode_length-1) {
            index_mismatches = minval;
        }
        swap(column, prevcol);
        if (first > last) {
            //cout << "RET: last > first" << endl;
            return make_pair(index_mismatches, 0);
        }
    }
    int best_len = 0, minval = max_val_align;
    //for (int i=first; i<=last; ++i) { // Prefer long match when equal score
    for (int i=last; i>=first; --i) { // Prefer short match when equal score
        if (prevcol[i+1] < minval) {
            minval = prevcol[i+1];
            best_len = i+1;
        }
    }
    //cout << "RET: loop terminated normally" << endl;
    return make_pair(index_mismatches, best_len);
}

