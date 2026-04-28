import os, sys, argparse, time, subprocess, re
#memory = 50000

#memory = 100000
#memory = 30000
#memory = 8000
memory = 5000

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
    parser.add_argument('-t', '--test', nargs=1, help="the data test input file")
    parser.add_argument('-o', '--outdir', nargs=1, help="the output directory")
    parser.add_argument('-s', '--sampling', default=50000, help="The subsampling rate")
    parser.add_argument('-c', '--calibrator', default="sms", help="The calibrator type")
    parser.add_argument("--cpu", default = 1, help="Number of cpus or cores to request - this will then be the number of cores over which the calibration is split up")
    parser.add_argument('--lambdas', nargs='+', help='The values to sample for the lambda variables')
    parser.add_argument('--rhos', nargs='+', help='The values to sample for the rho variables')
    parser.add_argument('--taus', nargs='+', help='The values to sample for the tau variables')
    parser.add_argument("--rhc",action='store_true', help="Flag determining if we're looking at a RHC/anti-nu playlist'")

    args = parser.parse_args()

    inputfile = args.input[0]
    outDir=args.outdir[0]
    samplerate = args.sampling
    calibrator = args.calibrator
    cpus = args.cpu
    testfile = args.test[0]
    lambdas = []
    if not args.lambdas or len(args.lambdas)==0:
        lambdas = [0.01, 0.1, 1, 10, 100]
    else:
        lambdas = args.lambdas
    rhos = []
    if not args.rhos or len(args.rhos)==0:
        rhos = [-1, -0.5, 0, 0.5, 1]
    else:
        rhos = args.rhos
    taus = []
    if not args.taus or len(args.taus)==0:
        taus = [0.5, 1, 1.5, 2]
    else:
        taus = args.taus

    print("Searching lamba values over range: ", lambdas)
    print("Searching rho values over range: ", rhos)
    print("Searching tau values over range: ", taus)
    
    lambda_str = " ".join(map(str, lambdas))
    rho_str = " ".join(map(str, rhos))
    tau_str = " ".join(map(str, taus))

    #Checking if the target output directory exists and if not attempting to create it
    if not (os.path.isdir(outDir)):
        print("Cannot access target output directory at %s\nAttempting to create it" % outDir)
        os.mkdir(outDir)


    #Lets break down the grid jobs by rho and tau - scratch that, lets do it by the lambdas
    for lambda_int in lambdas:
        for lambda_diag in lambdas:
            for lambda_offdiag in lambdas:
                for rho in rhos:
                    for tau in taus:
                        # Create wrapper
                        #wrapper_name = f"wrapper-TrainOptimiser{lambda_int}_{lambda_diag}_{lambda_offdiag}.sh"
                        wrapper_name = f"wrapper-TrainOptimiser{lambda_int}_{lambda_diag}_{lambda_offdiag}_{rho}_{tau}.sh"
                        #wrapper_name = f"wrapper-TrainOptimiser{rho}_{tau}.sh"
                        wrapper_path = "/nashome/a/alhart/gridWrappers/"+wrapper_name
                        my_wrapper = open(wrapper_path,"w")
                        my_wrapper.write("#!/bin/bash\n")
                        my_wrapper.write("date\n")
                        my_wrapper.write("cd $CONDOR_DIR_INPUT\n")
                        my_wrapper.write("echo Setting up environment\n")
                        my_wrapper.write("source /cvmfs/larsoft.opensciencegrid.org/spack-v0.22.0-fermi/setup-env.sh\n")
                        my_wrapper.write("spack load root@6.28.12 arch=linux-almalinux9-x86_64_v3\n")
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
                        tname = os.path.basename(testfile)
                        my_wrapper.write("date\n")
                        #my_wrapper.write(f"python SMSCalibrator.py -i {fname} -d {tname} -s {samplerate} -c {calibrator} --lambda_int {lambda_int} --lambda_diag {lambda_diag} --lambda_offd {lambda_offdiag} --rho {rho_str} --tau {tau_str} -o {outDir} -t {cpus}")
                        my_wrapper.write(f"python SMSCalibrator.py -i {fname} -d {tname} -s {samplerate} -c {calibrator} --lambda_int {lambda_int} --lambda_diag {lambda_diag} --lambda_offd {lambda_offdiag} --rho {rho} --tau {tau} -o {outDir} -t {cpus}")
                        if args.rhc:
                            my_wrapper.write(" --rhc")
                        my_wrapper.write("\n")
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
                        cmd+= f" --cpu {cpus}"
                        cmd+= " --lines=+FERMIHTC_AutoRelease=True"
                        cmd+= " --lines=+FERMIHTC_GraceMemory=1024"
                        cmd+= " --lines=+FERMIHTC_GraceLifetime=1800"
                        cmd+= " --lines=+max_retries=3" #Retry if job fails
                        cmd+= " --use-pnfs-dropbox"
                        cmd+= f" -f dropbox:///exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/tools/SMSCalibrator.py"
                        cmd+= f" -f dropbox:///exp/minerva/app/users/alhart/MAT_AL9/MINERvA-101-Cross-Section/tools/requirements.txt"
                        cmd+= f" -f dropbox://{inputfile}"
                        cmd+= f" -f dropbox://{testfile}"
                        cmd+=(" file://" + wrapper_path)
                        print(cmd)
                        success = submit_with_retries(cmd, 5, 10)
                        #os.system(cmd)
                        #if os.path.exists(wrapper_path):
                            #os.remove(wrapper_path)
                        #else:
                        #    print("Wrapper file not created successfully")
                        print("Done")
