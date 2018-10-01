int bounded_levenshtein_distance(
    int max_val,
    int s1len, const string& s1,
    int s2len, const string& s2
    )
{
    int buf1[s2len+1], buf2[s2len+1];

    int* column = buf1, * prevcol = buf2;

    for (int i=0; i<s2len+1; i++) {
        prevcol[i] = i;
        column[i] = i;
    }

    int first = 0, last = min(max_val, s2len);
    for (int i = 0; i < s1len; i++) { // Column!
        int j;
        
        if (i >= max_val-1) {
            // First loop: Eliminate leading rows above threshold. Always move one
            // down, so first will be at least incremented by one.
            for (j = first; j < last; j++) { // Row
                // Only need to check diagonally and left, not from top, as it is already 
                // at the max.
                int current_val = min(
                    prevcol[1 + j] + 1,
                    prevcol[j] + (s1[i] == s2[j] ? 0 : 1)
                    );
                if (current_val < max_val) {
                    column[j+1] = current_val;
                    break;
                }
            }
            // Should now point to first row which is below max_val
            first = j+1;
        }
        else {
            column[0] = i+1;
            j = 0;
        }

        // Second loop: find new values in range; keep track of last row below
        // threshold, so we can possibly constrain the "last"
        int last_below = j;
        for (; j <= last && j < s2len; j++) {
            auto possibilities = {
                prevcol[1 + j] + 1,
                column[j] + 1,
                prevcol[j] + (s1[i] == s2[j] ? 0 : 1)
            };
            column[j+1] = std::min(possibilities);
            if (column[j+1] < max_val) {
                last_below = j+1;
            }
        }
        last = last_below;

        swap(column, prevcol);
        if (first > last)
            return max_val;
    }
    return prevcol[s2len];
}



/** aligned_s2_length
 * Finds the length of the best alignment of s2 to s1. s1 is considered
 * to be a template of fixed length.
 * s2 will be aligned to it, with a penalty for insertions, deletions and
 * substitutions. The length of the best match of s2 is returned. If there
 * are several equally scored alignments, the length of the longest one is
 * returned.
 * 
 * For the purpose of reducing computational time, this function takes
 * an argument maxval indicating the maximum number of edit operations
 * performed. This is easily computable in our case.
 */
int aligned_s2_length(int max_val, const string& s1, const string& s2)
{
    const int s1len = s1.length();
    const int s2len = s2.length();

    int buf1[s2len+1], buf2[s2len+1];

    int* column = buf1, * prevcol = buf2;

    for (int i=0; i<s2len+1; i++) {
        prevcol[i] = i;
        column[i] = i;
    }

    int first = 0, last = min(max_val, s2len);
    for (int i = 0; i < s1len; i++) { // Column!
        int j;
        
        if (i >= max_val-1) {
            // First loop: Eliminate leading rows above threshold. Always move one
            // down, so first will be at least incremented by one.
            for (j = first; j < last; j++) { // Row
                // Only need to check diagonally and left, not from top, as it is already 
                // at the max.
                int current_val = min(
                    prevcol[1 + j] + 1,
                    prevcol[j] + (s1[i] == s2[j] ? 0 : 1)
                    );
                if (current_val < max_val) {
                    column[j+1] = current_val;
                    break;
                }
            }
            // Should now point to first row which is below max_val
            first = j+1;
        }
        else {
            column[0] = i+1;
            j = 0;
        }

        // Second loop: find new values in range; keep track of last row below
        // threshold, so we can possibly constrain the "last"
        int last_below = j;
        for (; j <= last && j < s2len; j++) {
            auto possibilities = {
                prevcol[1 + j] + 1,
                column[j] + 1,
                prevcol[j] + (s1[i] == s2[j] ? 0 : 1)
            };
            column[j+1] = std::min(possibilities);
            if (column[j+1] < max_val) {
                last_below = j+1;
            }
        }
        last = last_below;

        swap(column, prevcol);
        if (first > last)
            return max_val;
    }
    // prevcol now contains the best alignments of s2 to s1:
    // The value of prevcol[j] is the alignment score of s1 to s2, with j
    // being the position in s2 (TBC).

    // The longest alignment with the best score will be 
    int best = max_val;
    int bestj = 0;
    for (int al=first; al <= last; ++al) {
        if (prevcol[al] <= best) {
            bestj = al;
            best = prevcol[al];
        }
    }
    return bestj;
}

