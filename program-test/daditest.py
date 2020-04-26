import numpy
import dadi
from dadi import *
import sys

def IM(params, ns, pts):
    print("params: ",end="",file=sys.stderr)
    print(params, file=sys.stderr)
    T, nu1, nu2, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12 = m12, m21 = m21)

    fs = Spectrum.from_phi(phi, ns, (xx, xx))
    return fs

res = IM([float(x) for x in sys.argv[3:8]], (int(sys.argv[1]),int(sys.argv[1])), 200)

for i in range(0,int(sys.argv[1])+1):
    for j in range(0,int(sys.argv[2])+1):
        if isinstance(res[i,j], float):
            print(round(res[i,j],3), end="\t")
        else:
            print(0, end="\t")
    print("")
