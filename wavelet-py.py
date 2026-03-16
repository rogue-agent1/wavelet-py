#!/usr/bin/env python3
"""Haar wavelet transform (1D)."""
import sys,math

def haar_forward(x):
    n=len(x);assert n&(n-1)==0,"Length must be power of 2"
    a=list(x);h=n
    while h>1:
        h//=2
        for i in range(h):
            s=(a[2*i]+a[2*i+1])/math.sqrt(2)
            d=(a[2*i]-a[2*i+1])/math.sqrt(2)
            a[i],a[h+i]=s,d
    return a

def haar_inverse(c):
    n=len(c);a=list(c);h=1
    while h<n:
        tmp=list(a)
        for i in range(h):
            a[2*i]=(tmp[i]+tmp[h+i])/math.sqrt(2)
            a[2*i+1]=(tmp[i]-tmp[h+i])/math.sqrt(2)
        h*=2
    return a

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--test":
        x=[1,2,3,4,5,6,7,8]
        c=haar_forward(x)
        r=haar_inverse(c)
        assert all(abs(a-b)<1e-10 for a,b in zip(x,r))
        # Energy preservation
        e_orig=sum(v**2 for v in x);e_coeff=sum(v**2 for v in c)
        assert abs(e_orig-e_coeff)<1e-10
        # Simple case
        c2=haar_forward([1,1]);assert abs(c2[0]-math.sqrt(2))<1e-10;assert abs(c2[1])<1e-10
        print("All tests passed!")
    else:
        x=[1,4,2,3];c=haar_forward(x)
        print(f"Signal:  {x}");print(f"Coeffs:  {[f'{v:.3f}' for v in c]}")
if __name__=="__main__":main()
