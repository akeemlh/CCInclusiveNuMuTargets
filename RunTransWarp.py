#Tool used to loop through the directory structure of my analysis and combine the files produced by the grid submission scrips, when broken down into multiple smaller files

import ROOT
import os
import sys, argparse, math

parser = argparse.ArgumentParser( prog='RunTransWarp.py',
                description='automates the submission of runEventLoopTargetsNew to fermigrid')
parser.add_argument('inbase', nargs=1, help="the data input base location directory")
parser.add_argument('-n', "--nums", nargs='*', help="Numbers of iterations to run")
parser.add_argument('-w', "--warps", nargs='*', help="Warps to run")
parser.add_argument('-v', "--variables", nargs='*', help="Variables")
parser.add_argument("--keep", help="Don't remove the hadded files after")
parser.add_argument('-o', "--outdir", nargs=1, help="Output directory")
parser.add_argument('-l', "--playlists", nargs='*', help="Playlists to work on, leave empty to default to all neutrino mode (FHC) playlists",)
parser.add_argument('-t', "--targets", nargs='*', help="targes to work on, leave empty to default to all targets")
parser.add_argument('-p', "--petals", nargs='*', help="petals to work on, leave empty to default to all petals as well as sum")
args = parser.parse_args()

playlists = ["1A", "1B", "1C", "1D", "1E", "1F", "1G", "1L", "1M", "1N", "1O", "1P"]
if args.playlists is not None:
  playlists = args.playlists
petals = ["SUM", "-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]
if args.petals is not None:
  petals = args.petals
targets = ["SUM", "2026", "2082", "3006", "3026", "3082", "4082", "5026", "5082", "6000", "Iron", "Lead", "Carbon", "Water"]
if args.targets is not None:
  targets = args.targets
variables = args.variables
warps = args.warps
nums = args.nums
outdir = args.outdir[0]
inbase = args.inbase[0]
keep = False
if args.keep is not None:
  keep = True

targetsToSum = {}
for target in targets:
  if target == "Iron":
    targetsToSum["Iron"]=["2026", "3026", "5026"]
  elif target == "Lead":
    targetsToSum["Lead"]=["2082", "3082", "4082", "5082"]
  elif target == "Carbon":
    targetsToSum["Carbon"]=["3006"]
  elif target == "Water":
    targetsToSum["Water"]=["6000"]
  elif target == "SUM": #If we're doing everything
    targetsToSum["Iron"]=["2026", "3026", "5026"]
    targetsToSum["Lead"]=["2082", "3082", "4082", "5082"]
    targetsToSum["Carbon"]=["3006"]
    targetsToSum["Water"]=["6000"]
    targetsToSum["2026"]=["2026"]
    targetsToSum["2082"]=["2082"]
    targetsToSum["3006"]=["3006"]
    targetsToSum["3026"]=["3026"]
    targetsToSum["3082"]=["3082"]
    targetsToSum["4082"]=["4082"]
    targetsToSum["5026"]=["5026"]
    targetsToSum["5082"]=["5082"]
    targetsToSum["6000"]=["6000"]
  elif target.isDigit(): #Maybe add more checks, make sure it's in the right range (>2026 <6000 eg?) #What about pseudotargets?
    targetsToSum[target]=[target]

petalsToSum = {}
for petal in petals:
  if petal == "SUM": #If we're doing everything
    for i in range(-1, 11):
      print(i) #Just to check we're doing the full range
      i_str = str(i)
      petalsToSum[i_str]=[i_str]
  elif petal.isDigit(): #Maybe add more checks, make sure it's in the right numerical range
    petalsToSum[petal]=[petal]

tmpfolder = outdir+"/tmpSummed/"

mergeFiles = {}
for key, value in petalsToSum.items():
    for val in value:
        datapath = inbase + "/runEventLoopTrackerData_petal_"+val+".root"
        mcpath = inbase + "/runEventLoopTrackerMC_petal_"+val+".root"
        if not os.path.exists(datapath):
          print(datapath, " doesn't exist")
        if not os.path.exists(mcpath):
          print(mcpath, " doesn't exist")
        if os.path.exists(datapath) and os.path.exists(mcpath):
          outfiledata = tmpfolder+"/runEventLoopTrackerData_petal_"+key+".root"
          outfilemc = tmpfolder+"/runEventLoopTrackerMC_petal_"+key+".root"
          if outfiledata not in mergeFiles:
            mergeFiles[outfiledata]=[datapath]
            mergeFiles[outfilemc]=[mcpath]
          else:
            mergeFiles[outfiledata].append(datapath)
            mergeFiles[outfilemc].append(mcpath)

for key, value in targetsToSum.items():
    for val in value:
        datapath = inbase + "/runEventLoopTargetsData"+val+".root"
        mcpath = inbase + "/runEventLoopTargetsMC"+val+".root"
        if not os.path.exists(datapath):
          print(datapath, " doesn't exist")
        if not os.path.exists(mcpath):
          print(mcpath, " doesn't exist")
        if os.path.exists(datapath) and os.path.exists(mcpath):
          outfiledata = tmpfolder+"/runEventLoopTargetsData"+key+".root"
          outfilemc = tmpfolder+"/runEventLoopTargetsMC"+key+".root"
          if outfiledata not in mergeFiles and outfilemc not in mergeFiles:
            mergeFiles[outfiledata]=[datapath]
            mergeFiles[outfilemc]=[mcpath]
          else:
            mergeFiles[outfiledata].append(datapath)
            mergeFiles[outfilemc].append(mcpath)
          
for key, value in mergeFiles.items():
    merger = ROOT.TFileMerger()
    merger.OutputFile(key)
    for toAdd in value:
        merger.AddFile(toAdd)
    success = merger.Merge()
    if success:
        print(f"Successfully merged files into {key}")
    else:
        print("Merging failed")