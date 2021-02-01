#!/bin/python3
import sys
import getopt
import math
import numpy as np

def env(fc, s):
    if s ==  "home":
        return max_range(fc, 0.231,5)
    elif s == "mall":
        return max_range(fc, 0.047,12)
    elif s == "los":
        return max_range(fc, 0,0)

def max_range(fc, alpha, Lw):
    for x in range(1,10000):
        r = loss(x,fc,alpha,Lw) 
        if r > 102:
            return x

def get_loss_for_env(distance, fc, s):
    if s ==  "home":
        return loss(distance, fc, 0.231, 5)
    elif s == "mall":
        return loss(distance, fc, 0.047, 12)
    elif s == "los":
        return loss(distance, fc, 0, 0)

def loss(d,fc,alpha,Lw):
    np.seterr(over='ignore')
    freeatt=20*np.log10(d)
    wallatt=alpha*d*Lw
    fratt=20*math.log(fc/5,10)
    return( 46.4 + freeatt + fratt + wallatt)
loss_v=np.vectorize(loss)

def main(argv):
    fc=2.44
    alpha=0
    Lw=0
    Pt=20
    outputfile=''
    try:
        opts, args = getopt.getopt(argv, "hf:a:L:P:o:")
    except getopt.GetoptError:
        print ('range.py -f <Center Frequency> -a <Wall probability> -L <Wall attenuation> -P <Transmission Power> -o <Output File>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print ('help range.py -f <Center Frequency> -a <Wall probability> -L <Wall attenuation> -P <Transmission Power> -o <Output File>')
            sys.exit()
        elif opt in ("-f"):
            fc=float(arg)#float(args.pop(0))
        elif opt in ("-a"):
            alpha=float(arg)#float(args.pop(0))
            print(alpha)
        elif opt in ("-L"):
            Lw=float(arg)#float(args.pop(0))
        elif opt in ("-P"):
            Pt=float(arg)#float(args.pop(0))
        elif opt in ("-o"):
            outputfile=arg
    print("fc=",fc)
    print("alpha=",alpha)
    print("Lw=",Lw)
    print("Pt=",Pt)
    distance=np.arange(1,101)
    rp=Pt-loss_v(distance,fc,alpha,Lw)
    out=np.column_stack((distance,rp))
    print(out)
    fmt = '%d','%1.3f'
    try:
        np.savetxt(outputfile, out, delimiter=' ', fmt=fmt)
    except Exception as e:
        print(traceback.format_exception(*sys.exc_info()))
        sys.exit(-1)
if __name__ == "__main__":
    main(sys.argv[1:])
