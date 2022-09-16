## This is a general library of useful py functions

def rangeR(s,e,n):
    out = []
    for i in range(n):
        out.append((e-s)/(n-1)*i + s)
    return out
