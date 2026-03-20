#Tool used to loop through the directory structure of my analysis and combine the files produced by the grid submission scrips, when broken down into multiple smaller files

import ROOT
import os
import sys, argparse

parser = argparse.ArgumentParser( prog='combineTargets',
                description='automates the submission of runEventLoopTargetsNew to fermigrid')
parser.add_argument('indir', nargs=1, help="the data input directory")
parser.add_argument('-n', "--num", nargs=1, help="Number of subruns expected")
args = parser.parse_args()

indir = args.indir[0]
numSplits = int(args.num[0])
bases = []

for (root, dirs, file) in os.walk(indir):
    for f in file:
        if '_p' in f:
            fullbase = root+"/"+f[:(f.find("_p"))]
            if fullbase not in bases:
                bases.append(fullbase)
vettedBases = []
rejectedBases = []
for base in bases:
    count = 0
    for i in range(numSplits):
        probepath = base + "_n"+str(i)+".root"
        if os.path.exists(probepath):
            count+=1
        else:
            print(probepath, " doesn't exist")
    if count == numSplits:
        vettedBases.append(base)
    else:
        print(base, count)
        rejectedBases.append(base)

print("some files missing from following sets:", rejectedBases)

for base in vettedBases:
    print("Now working on set: ", base)
    outfile = base+".root"
    # Create the merger
    merger = ROOT.TFileMerger()
    merger.OutputFile(outfile)
    for i in range(numSplits):
        path = base + "_p"+str(i)+".root"
        merger.AddFile(path)
    # Loop over all ROOT files in the directory

    success = merger.Merge()

    f = ROOT.TFile(outfile, "UPDATE")
    mcpot = f.Get("MCPOT")
    datapot = f.Get("DataPOT")
    if mcpot:
        mcpot.SetVal(mcpot.GetVal()/numSplits)
        mcpot.Write("", ROOT.TObject.kOverwrite)
    if datapot:
        datapot.SetVal(datapot.GetVal()/numSplits)
        datapot.Write("", ROOT.TObject.kOverwrite)
    f.Close()

    if success:
        print(f"Successfully merged files into {outfile}")
    else:
        print("Merging failed")
