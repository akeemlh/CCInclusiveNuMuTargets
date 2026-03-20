import os, sys, argparse, time, subprocess, re

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
    epoch_time = int(time.time())
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
