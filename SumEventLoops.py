#Tool used to loop through the directory structure of my analysis and combine the files produced by the grid submission scrips, when broken down into multiple smaller files

import ROOT
import os
import sys, argparse, math

parser = argparse.ArgumentParser( prog='SumEventLoops.py',
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

base_dir_prefix = os.path.join(os.getenv("HOME"),"data","UnfoldingStudies","MigrationMatrix")
variations = ["_CV","_piontune","_2p2h","_rpa","_default","_rpapiontune","_2p2hrpa"]
#variatons = ["_CV","rpa"]
#_infile = "Migration_MakeFlux-1_ApplyFlux-1_Multiplicity-0_CombinedPlaylists-me1A.root"
#_infile = "Merged_MigrationMatrix_optimized_minervame5A_6A_6B.root"
#target_dir = "/minerva/data/users/bashyal8/UnfoldingStudies/pzpt_unfolding/"
#target_dir = "/minerva/data/users/"+os.getenv("USER")+"/UnfoldingStudies"
target_dir = os.path.join(os.getenv("HOME"),"data","UnfoldingStudies_hms")
reco_hist = "h_q2_enu_qelike_reco"  #"h_pzmu_ptmu_qelike_reco"
truth_hist = "h_q2_enu_qelike_truth"#"h_pzmu_ptmu_qelike_truth"
migration_hist = "h_q2_enu_qelike_migration"#"h_pzmu_ptmu_qelike_migration"
_output_prefix = reco_hist.replace("qelike_reco","")
output_prefix = _output_prefix.replace("h","stat_universe")
#variation = variations[2]
_infile = "MigrationMatrix_enu.root"
num_dim = str(2)
num_uni = str(25)
num_iter = "1,2,3,4,5,6,10"
#num_iter = "20"
#num_iter = "1,3,4,5,10"
max_chi2 = str(1000)
step_chi2 = str(10)


#needed hists
##########
# --data --> Data to be unfolded
# --data_truth -->True Data Distribution
# --reco -->Reconstructed MC distribution
# --truth -->Truth MC distribution
# --migration -->Migration matrix for Reconstructed MC to True MC
#########

for variation in variations:
  output_file = target_dir+"/"+output_prefix+variation+"_"+num_uni+".root"
  dat_file = base_dir_prefix+variation+"/"+_infile
  mc_file = base_dir_prefix+variations[0]+"/"+_infile #should be 0 mostly....
  mc_file = base_dir_prefix+variations[0]+"/"+_infile #should be 0 mostly....
  migration_file = mc_file #I think migration matrix is only constructed from MC info...

  cmd = "../TransWarpExtraction -o "+output_file\
    + " --data_file "+dat_file \
    + " --data "+reco_hist \
    + " --data_truth_file "+dat_file \
    + " --data_truth "+truth_hist \
    + " --reco_file "+mc_file \
    + " --reco "+reco_hist \
    + " --truth_file "+mc_file \
    + " --truth "+truth_hist \
    + " --migration_file "+migration_file \
    + " --migration "+migration_hist \
    + " --num_dim "+num_dim \
    + " --num_uni "+num_uni \
    + " --num_iter "+num_iter \
    + " --max_chi2 "+max_chi2 \
    + " --step_chi2 "+step_chi2 
#    + " --exclude_bin "

  print cmd

  os.system(cmd)
