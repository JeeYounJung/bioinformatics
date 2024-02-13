def overlap_length(left, right):
    """Returns the length of the longest suffix of left that is a prefix of right
    
    Args:
        left: a string
        right: a string
    Returns:
        An integer length of the longest overlap (0 if there is no overlap)
    """
    ans = 0
    for i in left:
        if right.startswith(i):
            left = left[left.index(i):]
            if left == right[:len(left)]:
                ans = len(left)
            else:
                left = left[left.index(i)+1:]
        else:
            continue
    return ans

def merge_ordered_reads(reads):
    """Returns the shortest superstring resulting from
    merging a list of ordered reads.
    
    Args:
        reads: a list of strings
    Returns:
        A string that is a shortest superstring of the ordered input read strings.
    """
    ans = ""
    oldI = ""
    if len(reads) == 0:
        return ans
    for i in reads:
        num = overlap_length(oldI, i)
        ans += i[num:]
        oldI = i
    return ans


def queue(reads):
    ans = []
    f = ''
    for i in reads:
        length = overlap(f, i)
        ans.append(length)
        f = i
    return sorted(ans, reverse=True)
            
