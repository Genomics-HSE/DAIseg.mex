import numpy as np

from functools import reduce
from operator import xor
from itertools import chain


# two set of intervals intersection 
def intersections(a,b):
    if b==[[]]:
        return a
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges
    
    


def exclude_gaps(set1, set2): #set2 is set of gaps [[],[],[]]

    l = sorted((reduce(xor, map(set, chain(set1 , set2)))))
    XOR=[l[i:i + 2] for i in range(0, len(l), 2)]

    return intersections(XOR, set1)
    
    

def point_in_set(p, m_a):
    if m_a==[[]]:
        return False
    def point_in(p, a):
        if p>=a[0] and p<=a[1]:
            return True
        else:
            return False
    f=False
    for j in m_a:
        f=point_in(p,j)
        if f==True:
            return f
    return f
