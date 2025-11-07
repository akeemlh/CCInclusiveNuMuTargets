import os, sys, argparse, time, subprocess, re
memory = 4000
lifetime = 36 #hours

epoch_time = int(time.time())

with open(f"failingJobsTracker{epoch_time}.txt", "w") as outfile:
    outfile.write("Failing jobs")

def submit_with_retries(command_in, max_retries, delay):
    #Sometimes, inexplicably, a perfectly valid jobsub_submit command will fail, this function loops the submission until it succeeds for up to max_retries times
    """
    Executes a command and retries on failure, parsing for a Job ID on success.

    Args:
        command (list): The command and its arguments as a list of strings.
        max_retries (int): The maximum number of times to retry the command.
        delay (int): The delay in seconds between retries.
    
    Returns:
        bool: True if the command succeeded, False otherwise.
    """
    command = command_in.split(" ")
    for attempt in range(max_retries):
        print(f"--- Attempt {attempt + 1} of {max_retries} ---")
        print(f"▶️  Executing: {' '.join(command)}")

        try:
            # Run the command, capturing output and not raising exceptions for failure
            result = subprocess.run(
                command,
                #capture_output=True,
                text=True,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT
            )

            # Check if the command was successful (exit code 0)
            if result.returncode == 0 and not ("Error: condor_submit exited with failed status code" in result.stdout):
                print("\n✅ Success! The job was submitted.")
                print("\n--- STDOUT ---")
                print(result.stdout)
                print("--------------")

                # Use regex to find the Job ID in the standard output
                job_id_match = re.search(r"JobsubJobId of first job: (\S+)", result.stdout)
                if job_id_match:
                    print(f"\n🚀 Job ID found: {job_id_match.group(1)}")
                
                return True # Exit function on success
            
            # If the command failed
            else:
                print(f"\n❌ Command failed with exit code {result.returncode}")
                print("\n--- STDERR ---")
                # Print stderr if it exists, otherwise provide a default message
                print(result.stdout or "No standard error output captured.")
                print("--------------")
                
                if attempt < max_retries - 1:
                    print(f"\nRetrying in {delay} seconds...")
                    time.sleep(delay)
                else:
                    print("\n🚫 All submission attempts have failed. Giving up.")
                    with open(f"failingJobsTracker{epoch_time}.txt", "a") as outfile:
                        outfile.write(command_in)
                    return False

        except FileNotFoundError:
            print(f"Fatal Error: The command '{command[0]}' was not found.")
            print("Please ensure 'jobsub_lite' is installed and accessible in your system's PATH.")
            sys.exit(1)
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            return False

    return False

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
    parser.add_argument('--susa2p2hwarp', action='store_true', help="Replacing LowRecoil2p2hReweighter with SuSAFromValencia2p2hReweighter")
    parser.add_argument('-N', '--NumSubruns', help="Break down this run into N smaller runs - Runs will then need to be merged (using hadd for example) afterwards. Used to speed up long jobs.") 
    parser.add_argument('-sN', '--specifySubrun', help="Specify which number of subrun to run, used when rerunning failed submissions, in conjuction with the option -N above. This will submit a job broken down by N, as above, but only run the subjob denoted by this value")
    parser.add_argument("-a", "--tarfile", help="Specify the tar archive - can save submission time. Tar using \"tar -cvzf {tarball output path}.tar.gz -C /exp/minerva/app/users/alhart/MAT_AL9/opt/ .\". Otherwise it will tar the current state of my analysis opt folder")
    parser.add_argument('-m', '--mnvTune', help="Specify the tune to apply formatted as 3 digit integer. Eg -m 431 is mnvTune v4.3.1")
    args = parser.parse_args()

    tarpath =args.tarfile
    print("tarpath: ", tarpath)
    if tarpath is None:
        #Tarring opts folder
        tarballfolder = "/pnfs/minerva/persistent/users/alhart/NuMuNukeIncl/TarredMATFramework/"+str(epoch_time)
        os.mkdir(tarballfolder)
        tarballpath = tarballfolder+"/opt.tar.gz"
        print("tarring into ", tarballpath)
        cmd = "tar -cvzf "+tarballpath+" -C /exp/minerva/app/users/alhart/MAT_AL9/opt/ ."
        os.system(cmd)
    else: 
        print("Using supplied tarball: ", tarpath)
        tarballpath = tarpath

    playlists=args.playlists
    dataInDir=args.data[0]
    mcInDir=args.mc[0]
    outDir=args.out[0]
    skipSys = args.skipsys
    no2p2hwarp = args.no2p2hwarp
    amudiswarp = args.amudiswarp
    lowq2warp = args.lowq2warp
    susa2p2hwarp=args.susa2p2hwarp
    tune = 100 #default to mnvtune100
    if args.mnvTune is not None:
        tune = args.mnvTune
    print("Using mnvTune", tune)
    if args.NumSubruns is not None and args.specifySubrun is None:
        numOfSubRuns = args.NumSubruns
    else:
        numOfSubRuns = 0
    if args.specifySubrun is not None:
        specifySubrun = args.specifySubrun
        specifySubrunTotal = args.NumSubruns
    else:
        specifySubrun = None
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
        print("No target petals specified - defaulting to all")
    print("Submitting targets: ", " ".join(petals))
    #exit()
    failingPlaylist = []
    failingPetals = []
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
            my_wrapper.write("echo Setting up environment\n")
            my_wrapper.write("export MINERVA_PREFIX=${INPUT_TAR_DIR_LOCAL}\n")
            if skipSys:
                my_wrapper.write("export MNV101_SKIP_SYST=True\n")
                memory=2000
            if no2p2hwarp:
                my_wrapper.write("export NO_2P2H_WARP=True\n")
            if amudiswarp:
                my_wrapper.write("export AMU_DIS_WARP=True\n")
            if lowq2warp:
                my_wrapper.write("export LOW_Q2_PION_WARP=True\n")
            if susa2p2hwarp:
                my_wrapper.write("export SUSA_2P2H_WARP=True\n")
            tuneStr ="export MnvTune="+tune+"\n"
            my_wrapper.write(tuneStr)
            if numOfSubRuns != 0:
                subrunstr = "export NumGridSubruns="+str(numOfSubRuns)+"n"
                my_wrapper.write(subrunstr)
            if specifySubrun!= None:
                subrunstr = "export NumGridSubruns="+str(specifySubrunTotal)+"n"
                my_wrapper.write(subrunstr)
                subrunstr = "export PROCESS="+str(specifySubrun)+"n"
                my_wrapper.write(subrunstr)
            my_wrapper.write("export XRD_NETWORKSTACK=IPv4\n")
            my_wrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh\n")
            #my_wrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh\n")
            my_wrapper.write("spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3\n")
            my_wrapper.write("spack load cmake@3.27.9%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
            my_wrapper.write("spack load gcc\n")
            my_wrapper.write("spack load python@3.9.15\n")
            my_wrapper.write("spack load ifdhc@2.8.0%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3\n")
            my_wrapper.write("spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
            my_wrapper.write("spack load py-numpy@1.24.3%gcc@12.2.0\n")
            my_wrapper.write("export JOBSUB_GROUP=minerva\n")
            #my_wrapper.write("htgettoken -a htvaultprod.fnal.gov -i minerva\n")
            #my_wrapper.write("export BEARER_TOKEN_FILE=/run/user/`id -u`/bt_u`id -u`\n")
            my_wrapper.write("source ${MINERVA_PREFIX}/bin/setup.sh\n")
            my_wrapper.write("export LD_LIBRARY_PATH=${ROOTSYS}/lib:${ROOTSYS}/lib/root:${MINERVA_PREFIX}/lib:${LD_LIBRARY_PATH}\n")
            my_wrapper.write("echo Setting up environment - DONE\n")
            my_wrapper.write("echo Running event loop\n")

            my_wrapper.write("MAX_RETRIES=3\n")
            my_wrapper.write("SLEEP=10\n")
            my_wrapper.write("attempt=1\n")

            my_wrapper.write("while true; do\n")
            my_wrapper.write("    echo \"Attempt ${attempt} of ${MAX_RETRIES}\"\n")
            my_wrapper.write("    ${MINERVA_PREFIX}/bin/runEventLoopTracker ${CONDOR_DIR_INPUT}/"+playlist+"-Data.txt "+"${CONDOR_DIR_INPUT}/"+playlist+"-MC.txt "+petal+"\n") #The "|| exit $?" is so that if this command fails the bash script will also exit which will trigger the condor job retry
            my_wrapper.write("    status=$?\n")
            my_wrapper.write("    if [ $status -eq 0 ]; then\n")
            my_wrapper.write("        echo \"Success on attempt $attempt\"\n")
            my_wrapper.write("        break\n")
            my_wrapper.write("    elif [ $attempt -ge $MAX_RETRIES ]; then\n")
            my_wrapper.write("        echo \"Failed after $MAX_RETRIES attempts (last exit code $status)\"\n")
            my_wrapper.write("        exit $status\n")
            my_wrapper.write("    fi\n")
            my_wrapper.write("\n")
            my_wrapper.write("    attempt=$((attempt+1))\n")
            my_wrapper.write("    echo \"Retrying in $SLEEP seconds...\"\n")
            my_wrapper.write("    sleep $SLEEP\n")
            my_wrapper.write("done\n")

            my_wrapper.write("echo Running event loop - DONE\n")
            my_wrapper.write("echo Copying files back to persistent\n")
            my_wrapper.write("ifdh cp ./runEventLoopTracker* "+outdirplaylist+"\n")
            my_wrapper.write("echo Copying files back to persistent - DONE\n")
            my_wrapper.write("echo SUCCESS\n")
            my_wrapper.close()
            logpath = outdirplaylist+"/LogPetal"+petal+".log"
            cmd = "jobsub_submit --group=minerva"
            cmd+= " --cmtconfig=x86_64-slc7-gcc49-opt -c has_avx2==True"
            cmd+= " --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest"
            cmd+= f" -L {logpath}"
            cmd+= f" --expected-lifetime {lifetime}h"
            cmd+= " --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --mail_always"
            cmd+= f" --memory {memory}MB"
            cmd+= f" --disk 15GB"
            cmd+= " --lines=+FERMIHTC_AutoRelease=True"
            cmd+= " --lines=+FERMIHTC_GraceMemory=1024"
            cmd+= " --lines=+FERMIHTC_GraceLifetime=1800"
            cmd+= " --lines=+max_retries=3" #Retry if job fails
            cmd+= f" -f dropbox://{dataInDir}/{playlist}-Data.txt -f dropbox://{mcInDir}/{playlist}-MC.txt"
            cmd+= f" --tar_file_name dropbox://{tarballpath}"
            #cmd+= f" -f /{dataInDir}/{playlist}-Data.txt -f /{mcInDir}/{playlist}-MC.txt"
            #cmd+= f" --tar_file_name /{tarballpath}"
            
            if numOfSubRuns != 0:
                subrunstr = " -N "+str(numOfSubRuns)
                cmd += subrunstr
            cmd+=(" file://" + wrapper_path)
            print(cmd)
            success = submit_with_retries(cmd, 5, 10)
            #os.system(cmd)
            #if os.path.exists(wrapper_path):
                #os.remove(wrapper_path)
            #else:
            #    print("Wrapper file not created successfully")
            print("Done")
            if not success:
                failingPlaylist.append(playlist)
                failingPetals.append(petal)

    if len(failingPetals)>0:
        print("Failed to submit jobs:")
        for i in range(len(failingPlaylist)):
            print("Playlist: ", failingPlaylist[i])
            print("Petal: ", failingPetals[i])
    else:
        print("All jobs submitted successfully")

    cmd = "rm -r "+tarballfolder
    os.system(cmd)
