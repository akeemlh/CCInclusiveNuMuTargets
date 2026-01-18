#Tool used to loop through the directory structure of my analysis and combine the files produced by the grid submission scrips, when broken down into multiple smaller files

import ROOT
import os
import sys, argparse, math

parser = argparse.ArgumentParser( prog='combineSplitRuns',
                description='automates the submission of runEventLoopTargetsNew to fermigrid')
parser.add_argument('indir', nargs=1, help="the data input directory")
parser.add_argument('-n', "--num", nargs=1, help="Number of subruns expected")
parser.add_argument('-o', "--outdir", help="Output parent directory")
parser.add_argument('-f', "--force", action='store_true', help="Overwrite existing files at output directory")
args = parser.parse_args()

indir = args.indir[0]
numSplits = int(args.num[0])
outdir = args.outdir
if outdir is not None and outdir[:-1]=="/":
    outdir = outdir[:-1]
bases = []
force = (args.force!=None)

for (root, dirs, file) in os.walk(indir):
    for f in file:
        idx = f.find("_n")
        if idx!=-1 and f[idx+2]!="e":
            print("f: ", f)
            fullbase = root+"/"+f[:(f.find("_n"))]
            if fullbase not in bases:
                print(fullbase, root, f)
                bases.append(fullbase)
            
print("Testing:", bases)
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
        #print(base, count)
        rejectedBases.append(base)

print("some files missing from following sets:", rejectedBases)

for base in vettedBases:
    print("Now working on set: ", base)
    lastslash = base.rfind("/")
    secondlastslash = base.rfind("/", 0, lastslash)
    playlist = base[secondlastslash+1:lastslash]
    filebase = base[lastslash+1:]
    """ print(lastslash)
    print(secondlastslash)
    print(base[secondlastslash+1:lastslash])
    continue """
    outbase = ""
    if outdir is None:
        outfile = base+".root"
    else:
        out_dir = f"{outdir}/{playlist}/"
        if not os.path.isdir(out_dir):
            print("Couldn't find directory %s, trying to make it " % out_dir )
            os.mkdir(out_dir)
        outfile = f"{out_dir}{filebase}.root"

    if os.path.exists(outfile):
        if force:
            print("File %s already exists. Overwriting it" % outfile)
        else:
            print("File %s already exists. Skipping. To overwrite existing files rerun with -f/--force" % outfile)
            continue
    print("Saving to %s" % outfile)
    # Create the merger
    merger = ROOT.TFileMerger()
    merger.OutputFile(outfile)
    if numSplits > 20:
        divs = math.ceil(numSplits / 20)
        for i in range(divs):
            start = i *20
            if i == divs:
                end = numSplits
            else:
                end = ((i+1) *20)-1
    if numSplits < 20:
        for i in range(numSplits):
            path = base + "_n"+str(i)+".root"
            merger.AddFile(path)
    # Loop over all ROOT files in the directory

    success = merger.Merge()
    if success:
        print(f"Successfully merged files into {outfile}")
    else:
        print("Merging failed")

    print(f"Fixing POTs")
    """ f = ROOT.TFile(outfile, "UPDATE")
    mcpot = f.Get("MCPOT")
    datapot = f.Get("DataPOT")
    if mcpot:
        mcpot.SetVal(mcpot.GetVal()/numSplits)
        mcpot.Write("", ROOT.TObject.kOverwrite)
    if datapot:
        datapot.SetVal(datapot.GetVal()/numSplits)
        datapot.Write("", ROOT.TObject.kOverwrite)
    f.Close() """
    print(f"Fixing POTs finished")
