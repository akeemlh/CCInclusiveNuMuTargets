import os

memmap = {"Tracker":5000, "Targets":5000}
lifetime = 24 #hours
outdir_logs = ""
#dirname  = "output11July0"
dirname  = "10OctTemp"
inputdir = "output18July0"
numIter = 5

#Tarring MAT opts folder
cmd = "tar -cvzf /exp/minerva/app/users/alhart/opt.tar.gz -C /exp/minerva/app/users/alhart/MAT_AL9/opt/ ."
os.system(cmd)

#If we already have a tarred opts folder, remove it so I can copy over new one (/persistent/ doesn't like being overwritten)
if (os.path.isfile("/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz")):
    cmd = "rm /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz"
    os.system(cmd)

#Copy opts tar
cmd = "cp /exp/minerva/app/users/alhart/opt.tar.gz /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz"
os.system(cmd)

playlists = ["1A", "1B", "1C", "1D", "1E", "1F", "1G", "1L", "1M", "1N", "1O", "1P", "Test"]
#playlists = ["Test1L", "Test1C"] #Used for rapid testing only, 1L is water filled, 1C is unfilled
#playlists = ["Test"]

cmd = "rmdir /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/Extracted1DXSecs"
os.system(cmd)
cmd = "mkdir /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/Extracted1DXSecs"
os.system(cmd)
#Which event loop(s) to run
#sets=["Tracker", "Targets"]
sets=["Tracker"]
for runType in sets:
    memory = memmap[runType]
    for playlist in playlists:
        cmd = "mkdir /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/Extracted1DXSecs/"+playlist
        os.system(cmd)
        # Create wrapper
        wrapper_name = "wrapper-extractxsec-"+runType+playlist+".sh"
        my_wrapper = open(wrapper_name,"w")
        my_wrapper.write("#!/bin/bash\n")
        my_wrapper.write("cd $CONDOR_DIR_INPUT\n")
        my_wrapper.write("mkdir opt\n")
        my_wrapper.write("echo Untarring\n")
        my_wrapper.write("tar -xvzf opt.tar.gz -C opt\n")
        my_wrapper.write("echo Untarring - DONE\n")
        my_wrapper.write("echo Setting up environment\n")
        my_wrapper.write("export XRD_NETWORKSTACK=IPv4\n")
        my_wrapper.write("export MINERVA_PREFIX=${CONDOR_DIR_INPUT}/opt\n")
        my_wrapper.write("source ${MINERVA_PREFIX}/bin/setup.sh\n")
        my_wrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh\n")
        my_wrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
        my_wrapper.write("spack load root@6.28.12\n")
        my_wrapper.write("spack load cmake\n")
        my_wrapper.write("spack load gcc\n")
        my_wrapper.write("spack load fife-utils\n")
        my_wrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${MINERVA_PREFIX}/lib:${LD_LIBRARY_PATH}\n")
        my_wrapper.write("echo Setting up environment - DONE\n")
        my_wrapper.write("echo Running xec loop\n")

        my_wrapper.write("${MINERVA_PREFIX}/bin/Extract1DCrossSection" +runType+ " "+str(numIter)+" ${CONDOR_DIR_INPUT}/"+playlist+"-runEventLoopDataTracker.root " + " ${CONDOR_DIR_INPUT}/"+playlist+"-runEventLoopMCTracker.root\n")
        my_wrapper.write("echo Running xsec loop - DONE\n")
        my_wrapper.write("echo Copying files back to persistent\n")
        ##Very temporary, find some other way that doesn invovle hardcoding tracker vvvvv
        my_wrapper.write("ifdh cp  ./GENIEXSECEXTRACT_"+playlist+"-MC.txt.root //pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/"+playlist+"-GENIEXSEC_"+runType+".root\n")
        my_wrapper.write("ifdh cp  ./tracker_* /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/Extracted1DXSecs/"+playlist+ "/\n")
        my_wrapper.write("echo Copying files back to persistent - DONE\n")
        my_wrapper.write("echo SUCCESS\n")
        my_wrapper.close()
        cmd = "jobsub_submit --group=minerva --cmtconfig=x86_64-slc7-gcc49-opt --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest --expected-lifetime %sh --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --mail_always --memory %dMB  -f /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/%s/%s-runEventLoopDataTracker.root -f /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/%s/%s-runEventLoopMCTracker.root -f /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz  file://%s/%s" % ( lifetime, memory, inputdir, playlist, inputdir, playlist, os.environ["PWD"] ,wrapper_name )    
        print(cmd)
        os.system(cmd)
        print("Done")