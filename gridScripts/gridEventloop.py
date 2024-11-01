import os
memmap = {"Tracker":15000, "Targets":15000, "FullDet":15000, "TrackerDSWater":15000, "TrackerUSWater":15000}
lifetime = 24 #hours
outdir_logs = ""
#dirname  = "output11July0"
#dirname  = "output18July0"
dirname = "output28Oct_2"

playlistDir = "/exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/PlaylistFiles/me-playlists"

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

#Which event loop(s) to run
sets=["Targets"]
#sets=["Tracker", "Targets", "FullDet", "TrackerDSWater", "TrackerUSWater"]
for runType in sets:
    memory = memmap[runType]
    for playlist in playlists:
        # Create wrapper
        wrapper_name = "wrapper-"+runType+playlist+".sh"
        wrapper_path = "/nashome/a/alhart/gridWrappers/"+wrapper_name
        my_wrapper = open(wrapper_path,"w")
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
        my_wrapper.write("spack load fife-utils@3.7.4\n")
        my_wrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${MINERVA_PREFIX}/lib:${LD_LIBRARY_PATH}\n")
        my_wrapper.write("echo Setting up environment - DONE\n")
        my_wrapper.write("echo Running event loop\n")
        my_wrapper.write("${MINERVA_PREFIX}/bin/runEventLoop" +runType+" ${CONDOR_DIR_INPUT}/"+playlist+"-Data.txt "+"${CONDOR_DIR_INPUT}/"+playlist+"-MC.txt\n")
        my_wrapper.write("echo Running event loop - DONE\n")
        my_wrapper.write("echo Copying files back to persistent\n")
        my_wrapper.write("ifdh cp ./runEventLoop"+runType+"Data.root //pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/"+playlist+"-runEventLoopData"+runType+".root\n")
        my_wrapper.write("ifdh cp ./runEventLoop"+runType+"MC.root //pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/"+dirname+"/"+playlist+"-runEventLoopMC"+runType+".root\n")
        my_wrapper.write("echo Copying files back to persistent - DONE\n")
        my_wrapper.write("echo SUCCESS\n")
        my_wrapper.close()
        cmd = "jobsub_submit --group=minerva --cmtconfig=x86_64-slc7-gcc49-opt --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest --expected-lifetime %sh --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --mail_always --memory %dMB --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=1024' --lines '+FERMIHTC_GraceLifetime=1800' -f -f dropbox://%s/Data/%s-Data.txt -f -f dropbox://%s/MC/%s-MC.txt -f /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz  file://%s" % ( lifetime, memory, playlistDir, playlist, playlistDir, playlist ,wrapper_path )    
        print(cmd)
        os.system(cmd)
        if os.path.exists(wrapper_path):
            os.remove(wrapper_path)
        else:
            print("Wrapper file not created successfully")
        print("Done")
