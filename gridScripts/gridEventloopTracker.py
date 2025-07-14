import os, sys, argparse, time
memory = 6000
lifetime = 12 #hours

if __name__ == '__main__':
    parser = argparse.ArgumentParser( prog='gridEventLoopTracker',
                    description='automates the submission of runEventLoopTracker to fermigrid')
    parser.add_argument('data', nargs=1, help="the data input directory")
    parser.add_argument('mc', nargs=1, help="the data input directory")
    parser.add_argument('out', nargs=1, help="the output directory")
    parser.add_argument('-p', '--playlists', nargs='*', help="the playlists to run")
    parser.add_argument('-t', '--petals', nargs='*', help="the petals to run")
    parser.add_argument('-s', '--skipsys', action='store_true', help="skip systematics? Used in warping studies")
    parser.add_argument('--no2p2hwarp', action='store_true', help="Turning off LowRecoil2p2hReweighter")
    parser.add_argument('--amudiswarp', action='store_true', help="Turning on AMUDISReweighter")
    parser.add_argument('--lowq2warp', action='store_true', help="Turning on LowQ2PiReweighter")
    args = parser.parse_args()
    #Tarring MAT opts folder
    epoch_time = int(time.time())
    tarballfolder = "/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/"+str(epoch_time)
    os.mkdir(tarballfolder)
    tarballpath = tarballfolder+"/opt.tar.gz"
    print("tarring into ", tarballpath)
    cmd = "tar -cvzf "+tarballpath+" -C /exp/minerva/app/users/alhart/MAT_AL9/opt/ ."
    os.system(cmd)

    playlists=args.playlists
    dataInDir=args.data[0]
    mcInDir=args.mc[0]
    outDir=args.out[0]
    skipSys = args.skipsys
    no2p2hwarp = args.no2p2hwarp
    amudiswarp = args.amudiswarp
    lowq2warp = args.lowq2warp

    #Checking all input directories exist
    if not (os.path.isdir(dataInDir)):
        print("Cannot access data input directory at %s" % dataInDir)
        exit(1)
    if not (os.path.isdir(mcInDir)):
        print("Cannot access mc input directory at %s" % mcInDir)
        exit(1)
    #Checking if the target output directory exists and if not attempting to create it
    if not (os.path.isdir(outDir)):
        print("Cannot access target output directory at %s\nAttempting to create it" % outDir)
        os.mkdir(outDir)
    if not playlists:
        playlists=["1A", "1B", "1C", "1D", "1E", "1F", "1G", "1L", "1M", "1N", "1O", "1P"]
        print("No playlists specified, using default set")

    vettedPlaylists = []
    #Checking all playlists actually correspond to a file within the data and mc input directories
    for playlist in playlists:
        datafilename = dataInDir+"/"+playlist+"-Data.txt"
        datafilefound = os.path.exists(datafilename)
        mcfilename = mcInDir+"/"+playlist+"-MC.txt"
        mcfilefound = os.path.exists(mcfilename)
        if not (datafilefound):
            print("Cannot access data playlist file at %s. Continuing without it" % datafilename)
        if not (mcfilefound):
            print("Cannot access MC playlist file at %s. Continuing without it"% mcfilename)
        if datafilefound and mcfilefound:
            vettedPlaylists.append(playlist)
            print("Data and MC files found for playlist %s. Using it"% playlist)


    print("Submitting playlists: ", vettedPlaylists)
    petals = args.petals
    if not petals:
        petals = ["-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]
        print("No target specified defaulting to all")
    print("Submitting targets: ", " ".join(petals))
    #exit()

    for petal in petals:
        for playlist in vettedPlaylists:
            outdirplaylist = outDir+"/"+playlist
            #Checking if the target output directory exists and if not attempting to create it
            if not (os.path.isdir(outdirplaylist)):
                print("Cannot access target output directory at %s\nAttempting to create it" % outdirplaylist)
                os.mkdir(outdirplaylist)
            else:
                #We cannot overwrite on persistent so I'm making sure the files we want to write out dont exist
                migrationpath = outdirplaylist+"/runEventLoopTrackerMigration_petal_"+petal+".root"
                datapath = outdirplaylist+"/runEventLoopTrackerData_petal_"+petal+".root"
                mcpath = outdirplaylist+"/runEventLoopTrackerMC_petal_"+petal+".root"
                """ if (os.path.exists(migrationpath)):
                    os.remove(migrationpath)
                if (os.path.exists(datapath)):
                    os.remove(datapath)
                if (os.path.exists(mcpath)):
                    os.remove(mcpath) """
            # Create wrapper
            wrapper_name = "wrapper-EvLoopTracker"+playlist+petal+".sh"
            wrapper_path = "/nashome/a/alhart/gridWrappers/"+wrapper_name
            my_wrapper = open(wrapper_path,"w")
            my_wrapper.write("#!/bin/bash\n")
            my_wrapper.write("cd $CONDOR_DIR_INPUT\n")
            my_wrapper.write("mkdir opt\n")
            my_wrapper.write("echo Untarring\n")
            my_wrapper.write("tar -xvzf opt.tar.gz -C opt\n")
            my_wrapper.write("echo Untarring - DONE\n")
            my_wrapper.write("echo Setting up environment\n")
            if skipSys:
                my_wrapper.write("export MNV101_SKIP_SYST=True\n")
            if no2p2hwarp:
                my_wrapper.write("export NO_2P2H_WARP=True\n")
            if amudiswarp:
                my_wrapper.write("export AMU_DIS_WARP=True\n")
            if lowq2warp:
                my_wrapper.write("export LOW_Q2_PION_WARP=True\n")
            my_wrapper.write("export XRD_NETWORKSTACK=IPv4\n")
            my_wrapper.write("export MINERVA_PREFIX=${CONDOR_DIR_INPUT}/opt\n")
            my_wrapper.write("source ${MINERVA_PREFIX}/bin/setup.sh\n")
            my_wrapper.write("source /cvmfs/minerva.opensciencegrid.org/minerva/setup/setup_minerva_products.sh\n")
            my_wrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
            my_wrapper.write("spack load root@6.28.12\n")
            #my_wrapper.write("spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
            my_wrapper.write("spack load gcc\n")
            my_wrapper.write("spack load fife-utils@3.7.4\n")
            my_wrapper.write("htgettoken -a htvaultprod.fnal.gov -i minerva\n")
            my_wrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${MINERVA_PREFIX}/lib:${LD_LIBRARY_PATH}\n")
            my_wrapper.write("echo Setting up environment - DONE\n")
            my_wrapper.write("ls -la\n")
            my_wrapper.write("echo $(ls -la)\n")
            my_wrapper.write("echo Running event loop\n")
            my_wrapper.write("${MINERVA_PREFIX}/bin/runEventLoopTracker ${CONDOR_DIR_INPUT}/"+playlist+"-Data.txt "+"${CONDOR_DIR_INPUT}/"+playlist+"-MC.txt "+petal+"\n")
            my_wrapper.write("echo Running event loop - DONE\n")
            my_wrapper.write("echo Copying files back to persistent\n")
            my_wrapper.write("ifdh cp ./runEventLoopTracker* "+outdirplaylist+"\n")
            my_wrapper.write("echo Copying files back to persistent - DONE\n")
            my_wrapper.write("echo SUCCESS\n")
            my_wrapper.close()
            cmd = "jobsub_submit --group=minerva --cmtconfig=x86_64-slc7-gcc49-opt  -c has_avx2==True --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest --expected-lifetime %sh --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --mail_always --memory %dMB --disk 15GB --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=1024' --lines '+FERMIHTC_GraceLifetime=1800' -f dropbox://%s/%s-Data.txt -f dropbox://%s/%s-MC.txt -f dropbox://%s  file://%s " % ( lifetime, memory, dataInDir, playlist, mcInDir, playlist , tarballpath, wrapper_path )    
            print(cmd)
            os.system(cmd)
            #if os.path.exists(wrapper_path):
                #os.remove(wrapper_path)
            #else:
            #    print("Wrapper file not created successfully")
            print("Done")

    cmd = "rm -r "+tarballfolder
    os.system(cmd)
