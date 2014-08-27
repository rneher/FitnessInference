import glob,os

dirlist = glob.glob('../data/N_*')

for dirname in dirlist:
    flist = glob.glob(dirname+'/20140418_seqs*_d_0.5.dat')
    for fname in flist:
        print fname
        os.remove(fname)
