import os, sys, argparse
memory = 6000
lifetime = 24 #hours

if __name__ == '__main__':
    parser = argparse.ArgumentParser( prog='gridEventLoopTargets',
                    description='automates the submission of runEventLoopTargetsNew to fermigrid')
    parser.add_argument('data', nargs=1, help="the data input directory")
    parser.add_argument('mc', nargs=1, help="the data input directory")
    parser.add_argument('out', nargs=1, help="the output directory")
    parser.add_argument('-p', '--playlists', nargs='*', help="the playlists to run")
    parser.add_argument('-t', '--targets', nargs='*', help="the target materials to run")
    args = parser.parse_args()
    #Tarring MAT opts folder
    cmd = "tar -cvzf /exp/minerva/app/users/alhart/opt.tar.gz -C /exp/minerva/app/users/alhart/MAT_AL9/opt/ ."
    os.system(cmd)

    #If we already have a tarred opts folder, remove it so I can copy over new one (/persistent/ doesn't like being overwritten)
    if (os.path.isfile("/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz")):
        cmd = "rm /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz"
        os.system(cmd)

    #Copy opts tar
    cmd = "mv /exp/minerva/app/users/alhart/opt.tar.gz /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz"
    os.system(cmd)

    playlists=args.playlists
    dataInDir=args.data[0]
    mcInDir=args.mc[0]
    outDir=args.out[0]

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
    targets = args.targets
    if not targets:
        targets = [""]
        print("No target specified")
    print("Submitting targets: ", " ".join(targets))
    #exit()

    for target in targets:
        for playlist in vettedPlaylists:
            outdirplaylist = outDir+"/"+playlist
            #Checking if the target output directory exists and if not attempting to create it
            if not (os.path.isdir(outdirplaylist)):
                print("Cannot access target output directory at %s\nAttempting to create it" % outdirplaylist)
                os.mkdir(outdirplaylist)
            else:
                #We cannot overwrite on persistent so I'm making sure the files we want to write out dont exist
                migrationpath = outdirplaylist+"/runEventLoopTracker2DMigration.root"
                datapath = outdirplaylist+"/runEventLoopTrackerData.root"
                mcpath = outdirplaylist+"/runEventLoopTrackerMC.root"
                if (os.path.exists(migrationpath)):
                    os.remove(migrationpath)
                if (os.path.exists(datapath)):
                    os.remove(datapath)
                if (os.path.exists(mcpath)):
                    os.remove(mcpath)
            # Create wrapper
            wrapper_name = "wrapper-EvLoopTargets"+playlist+target+".sh"
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
            #my_wrapper.write("spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
            my_wrapper.write("spack load gcc\n")
            my_wrapper.write("spack load fife-utils@3.7.4\n")
            my_wrapper.write("htgettoken -a htvaultprod.fnal.gov -i minerva\n")
            my_wrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${MINERVA_PREFIX}/lib:${LD_LIBRARY_PATH}\n")
            my_wrapper.write("echo Setting up environment - DONE\n")
            my_wrapper.write("ls -la\n")
            my_wrapper.write("echo $(ls -la)\n")
            my_wrapper.write("echo Running event loop\n")
            my_wrapper.write("${MINERVA_PREFIX}/bin/runEventLoopTracker ${CONDOR_DIR_INPUT}/"+playlist+"-Data.txt "+"${CONDOR_DIR_INPUT}/"+playlist+"-MC.txt\n")
            my_wrapper.write("echo Running event loop - DONE\n")
            my_wrapper.write("echo Copying files back to persistent\n")
            my_wrapper.write("ifdh cp ./runEventLoopTargets* "+outdirplaylist+"\n")
            my_wrapper.write("echo Copying files back to persistent - DONE\n")
            my_wrapper.write("echo SUCCESS\n")
            my_wrapper.close()
            cmd = "jobsub_submit --group=minerva --cmtconfig=x86_64-slc7-gcc49-opt  -c has_avx2==True --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest --expected-lifetime %sh --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --mail_always --memory %dMB --lines '+FERMIHTC_AutoRelease=True' --lines '+FERMIHTC_GraceMemory=1024' --lines '+FERMIHTC_GraceLifetime=1800' -f dropbox://%s/%s-Data.txt -f dropbox://%s/%s-MC.txt -f /pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/opt.tar.gz  file://%s " % ( lifetime, memory, dataInDir, playlist, mcInDir, playlist ,wrapper_path )    
            print(cmd)
            os.system(cmd)
            #if os.path.exists(wrapper_path):
                #os.remove(wrapper_path)
            #else:
            #    print("Wrapper file not created successfully")
            print("Done")
