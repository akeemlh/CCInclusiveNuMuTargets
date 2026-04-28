#Tool used to loop through the directory structure of my analysis and combine the files produced by the grid submission scrips, when broken down into multiple smaller files

import ROOT
import os
import sys, argparse, math
import PlotUtils 

VERBOSITY = True
def merge_rootfile(target_dir, source_files, current_path=""):
    """
    Recursively merges directories, trees, and PlotUtils histograms 
    from a list of source files into a target directory.
    """
    first_source = source_files[0]
    first_source.cd(current_path)

    # Prevent ROOT from globally tracking histograms in memory
    status = ROOT.TH1.AddDirectoryStatus()
    ROOT.TH1.AddDirectory(False)

    current_sourcedir = ROOT.gDirectory
    keys = current_sourcedir.GetListOfKeys()

    # Loop over all keys in this directory
    for key in keys:
        first_source.cd(current_path)
        obj = key.ReadObj()
        
        if not obj:
            continue

        obj_name = key.GetName()
        obj_class = obj.ClassName()
        print("obj_name: ", obj_name, " obj_class: ", obj_class)

        # Merge Histograms (including MnvH1D and MnvH2D)
        if "MnvH1D" in obj_class or "MnvH2D" in obj_class or obj.InheritsFrom("TH1"):
            if VERBOSITY:
                print(f"Merging {obj_class}: {obj_name}")
            
            for next_source in source_files[1:]:
                next_source.cd(current_path)
                key2 = ROOT.gDirectory.GetListOfKeys().FindObject(obj_name)
                if key2:
                    h2 = key2.ReadObj()
                    if h2:
                        obj.Add(h2)
                        # FIX 1: Explicitly delete the temporary histogram from C++ memory
                        h2.Delete()
            
            target_dir.cd()
            obj.Write(obj_name)
            # FIX 2: Delete the main accumulated object after writing to disk
            obj.Delete()

        # Merge TParameter objects
        elif obj_class.startswith("TParameter<"):
            if VERBOSITY:
                print(f"Merging {obj_class}: {obj_name}")
            
            # Extract the initial value
            total_val = obj.GetVal()
            
            for next_source in source_files[1:]:
                next_source.cd(current_path)
                key2 = ROOT.gDirectory.GetListOfKeys().FindObject(obj_name)
                if key2:
                    param2 = key2.ReadObj()
                    if param2:
                        # Add the value from the current file
                        total_val += param2.GetVal()
                        param2.Delete()
            
            # Set the accumulated value and save
            obj.SetVal(total_val)
            target_dir.cd()
            obj.Write(obj_name)
            obj.Delete()

        # Merge TTrees
        elif obj.InheritsFrom("TTree"):
            if VERBOSITY:
                print(f"Merging TTree: {obj_name}")
            glob_chain = ROOT.TChain(obj_name)
            
            for source in source_files:
                tree_path = f"{source.GetName()}/{current_path}/{obj_name}" if current_path else f"{source.GetName()}/{obj_name}"
                glob_chain.Add(tree_path)

            target_dir.cd()
            glob_chain.Merge(target_dir.GetFile(), 0, "keep")
            # FIX 3: Free the TChain
            glob_chain.Delete()

        # Merge Subdirectories (Recursion)
        elif obj.InheritsFrom("TDirectory"):
            print(f"Found subdirectory: {obj_name}")
            target_dir.cd()
            newdir = target_dir.mkdir(obj_name, obj.GetTitle())

            new_path = f"{current_path}/{obj_name}" if current_path else obj_name
            merge_rootfile(newdir, source_files, new_path)

        # Unhandled object types
        else:
            if VERBOSITY:
                print(f"Unhandled object type written directly, name: {obj_name}, class: {obj_class}")
            target_dir.cd()
            obj.Write(obj_name)
            obj.Delete()

    target_dir.SaveSelf(True)
    ROOT.TH1.AddDirectory(status)


playlist_divs = {"1A" : 20, #10?
                 "1B" : 5, 
                 "1C" : 5, 
                 "1D" : 20, 
                 "1E" : 15, 
                 "1F" : 30, #24
                 "1G" : 25, 
                 "1L" : 5, 
                 "1M" : 30, 
                 "1N" : 15, 
                 "1O" : 5, 
                 "1P" : 5, 
                 "5A" : 10, 
                 "6A" : 20, 
                 "6B" : 15, 
                 "6C" : 15, 
                 "6D" : 15, 
                 "6E" : 10, 
                 "6F" : 20, 
                 "6G" : 10, 
                 "6H" : 15, 
                 "6I" : 10, 
                 "6I_WaterIN" : 5, 
                 "6I_WaterOUT" : 5, 
                 "6J" : 15, 
                 }

#Double-time
#playlist_divs = {"1A" : 25, #10?
#                 "1B" : 10, 
#                 "1C" : 15, 
#                 "1D" : 35, 
#                 "1E" : 30, 
#                 "1F" : 35, #24
#                 "1G" : 35, 
#                 "1L" : 5, 
#                 "1M" : 50, 
#                 "1N" : 20, 
#                 "1O" : 10, 
#                 "1P" : 15, 
#                 "5A" : 15, 
#                 "6A" : 35, 
#                 "6B" : 25, 
#                 "6C" : 25, 
#                 "6D" : 25, 
#                 "6E" : 25, 
#                 "6F" : 20, 
#                 "6G" : 25, 
#                 "6H" : 39, 
#                 "6I" : 25, 
#                 "6I_WaterIN" : 15, 
#                 "6I_WaterOUT" : 10, 
#                 "6J" : 25, 
#                 }

parser = argparse.ArgumentParser( prog='combineSplitRuns',
                description='automates the submission of runEventLoopTargetsNew to fermigrid')
parser.add_argument('indir', nargs=1, help="the data input directory")
parser.add_argument('-n', "--num", nargs=1, help="Number of subruns expected. Enter -1 to use my preset binnings")
parser.add_argument('-o', "--outdir", help="Output parent directory")
parser.add_argument('-f', "--force", action='store_true', help="Overwrite existing files at output directory")
parser.add_argument('-v', "--verbose", action='store_true', help="Verbose output")
parser.add_argument('-p', "--playlist", help="Specify single playlist")
parser.add_argument('-d', "--dry_run", action='store_true', help="Don't actually perform merging")
args = parser.parse_args()

indir = args.indir[0]
numSplitsarg = int(args.num[0])
outdir = args.outdir
if outdir is not None and outdir[:-1]=="/":
    outdir = outdir[:-1]

if not (os.path.isdir(outdir)):
    print(f"Could not find target output directory {outdir}, so attempting to make it")
    os.mkdir(outdir)

bases = []
force = args.force
specified_playlist = args.playlist
VERBOSITY = args.verbose
dry_run = args.dry_run
ROOT.EnableImplicitMT()  # Enable threading - NOTE: If files are over slow network (XRootD) it might cause thread contention waiting for data.

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
    lastslash = base.rfind("/")
    print("lastslash: ", lastslash)
    secondlastslash = base.rfind("/", 0, lastslash)
    print("secondlastslash: ", secondlastslash)
    playlist = base[secondlastslash+1:lastslash]
    #playlist="1L"
    if specified_playlist is not None:
        playlist = specified_playlist
    print("playlist: ", playlist)
    if numSplitsarg==-1:
        numSplits=playlist_divs[playlist]
    else:
        numSplits = numSplitsarg
    print("playlist: ", playlist, "numSplits=", numSplits)
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
if dry_run:
    exit(0)
for base in vettedBases:
    print("Now working on set: ", base)
    lastslash = base.rfind("/")
    secondlastslash = base.rfind("/", 0, lastslash)
    playlist = base[secondlastslash+1:lastslash]
    #playlist="1L"
    if specified_playlist is not None:
        playlist = specified_playlist
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
    filesToMerge = []
    if numSplitsarg==-1:
        numSplits=playlist_divs[playlist]
    else:
        numSplits = numSplitsarg
    for i in range(numSplits):
        path = base + "_n"+str(i)+".root"
        filesToMerge.append(path)
    """ merger = ROOT.TFileMerger()
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
            filesToMerge.append(path)
    # Loop over all ROOT files in the directory

    success = merger.Merge()
    if success:
        print(f"Successfully merged files into {outfile}")
    else:
        print("Merging failed") """

    target_file = ROOT.TFile.Open(outfile, "RECREATE")
    if not target_file or target_file.IsZombie():
        print(f"Fatal: Could not open target file {outfile}")
        sys.exit(1)

    source_files = []
    for sf in filesToMerge:
        f = ROOT.TFile.Open(sf, "READ")
        if not f or f.IsZombie():
            print(f"Fatal: Could not open source file {sf}")
            sys.exit(1)
        source_files.append(f)
    print("merging into: ", target_file)
    merge_rootfile(target_file, source_files,)
    # Cleanup
    target_file.Close()
    for f in source_files:
        f.Close()

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

print("some files missing from following sets:", rejectedBases)