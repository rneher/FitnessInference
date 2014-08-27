#!/ebio/ag-neher/share/programs/bin/python2.7
import sys,os
import subprocess
import time

dt=int(sys.argv[-1])
print dt, sys.argv[-3]
for ni in range(int(sys.argv[-2])):
    gen = int(sys.argv[-3])+ni*dt
    if gen<25000:
        print gen
        args = sys.argv[1:-3] + [gen]
        cmd = '../src/analyze_toy_data.py '+' '.join(map(str, args))
        print cmd
        os.system(cmd)
    
    else:
        print "generation out of range", gen
