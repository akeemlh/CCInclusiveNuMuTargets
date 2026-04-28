import os, sys, argparse, time, subprocess, re
#memory = 50000
memory = 2000
#memory = 85000
#memory = 500000
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
    parser = argparse.ArgumentParser( prog='gridTrainCalibrator',
                    description='automates the submission of SMSCalibrator to fermigrid')

    parser.add_argument('-i', '--input', nargs=1, help="the mc input file")
    parser.add_argument('-o', '--outdir', nargs=1, help="the output directory")
    parser.add_argument('-s', '--sampling', default=0.1, help="The subsampling rate")
    parser.add_argument('-c', '--calibrator', default="sms", help="The calibrator type")
    parser.add_argument("-f", "--fillmode", type=int, help="Fill mode for filling probability space for .root inputs")
    args = parser.parse_args()

    inputfile = args.input[0]
    outDir=args.outdir[0]
    samplerate = args.sampling
    calibrator = args.calibrator
    fillmode = args.fillmode

    #Checking if the target output directory exists and if not attempting to create it
    if not (os.path.isdir(outDir)):
        print("Cannot access target output directory at %s\nAttempting to create it" % outDir)
        os.mkdir(outDir)


    # Create wrapper
    wrapper_name = "wrapper-TrainCalibrator.sh"
    wrapper_path = "/nashome/a/alhart/gridWrappers/"+wrapper_name
    my_wrapper = open(wrapper_path,"w")
    my_wrapper.write("#!/bin/bash\n")
    my_wrapper.write("date\n")
    my_wrapper.write("cd $CONDOR_DIR_INPUT\n")
    my_wrapper.write("echo Setting up environment\n")
    my_wrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh\n")
    my_wrapper.write("spack load python@3.9.15\n")
    my_wrapper.write("spack load ifdhc@2.8.0%gcc@12.2.0 arch=linux-almalinux9-x86_64_v3\n")
    my_wrapper.write("spack load ifdhc-config@2.6.20%gcc@11.4.1 arch=linux-almalinux9-x86_64_v2\n")
    my_wrapper.write("python -m venv venv\n")
    my_wrapper.write("source venv/bin/activate\n")
    my_wrapper.write("pip install -r requirements.txt\n")
    my_wrapper.write("date\n")
    fname = os.path.basename(inputfile)
    if ".gz" in fname:
        my_wrapper.write(f"gunzip {fname}\n")
        fname = fname[:-3]
    my_wrapper.write("date\n")
    my_wrapper.write(f"python SMSCalibrator.py -i {fname} -s {samplerate} -c {calibrator} -f {fillmode}\n")
    my_wrapper.write("echo Copying files back\n")
    my_wrapper.write("ifdh cp ./calibration_results/* "+outDir+"\n")
    my_wrapper.write("echo Copying files back - DONE\n")
    my_wrapper.write("echo SUCCESS\n")
    my_wrapper.close()
    logpath = outDir+"/LogTraining.log"
    cmd = "jobsub_submit --group=minerva"
    cmd+= " --cmtconfig=x86_64-slc7-gcc49-opt -c has_avx2==True"
    cmd += " --append_condor_requirements='(regexp(\".*fnpc7.*\",Machine)==FALSE)'"
    cmd+= " --singularity-image /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-el9:latest"
    #cmd+= f" -L {logpath}"
    cmd+= f" --expected-lifetime {lifetime}h"
    cmd+= " --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --role=Analysis --mail_always"
    cmd+= f" --memory {memory}MB"
    cmd+= f" --disk 60GB"
    cmd+= " --lines=+FERMIHTC_AutoRelease=True"
    cmd+= " --lines=+FERMIHTC_GraceMemory=1024"
    cmd+= " --lines=+FERMIHTC_GraceLifetime=1800"
    cmd+= " --lines=+max_retries=3" #Retry if job fails
    cmd+= " --use-pnfs-dropbox"
    cmd+= f" -f dropbox:///exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/tools/SMSCalibrator.py"
    cmd+= f" -f dropbox:///exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/tools/requirements.txt"
    cmd+= f" -f dropbox://{inputfile}"
    cmd+=(" file://" + wrapper_path)
    print(cmd)
    success = submit_with_retries(cmd, 5, 10)
    #os.system(cmd)
    #if os.path.exists(wrapper_path):
        #os.remove(wrapper_path)
    #else:
    #    print("Wrapper file not created successfully")
    print("Done")
