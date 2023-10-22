The container for running jobs on AWS Braket.
Clone this folder onto your local machine, install the pyqubo package into your current directory by opening Command Prompt or any terminal on the directory you cloned these two files into and run the command _pip install pyqubo -t ._
Similarly you can install any package into this container that is not present in the AWS Braket instance which you have imported for the code.
When running a job, make sure you carefully choose which QPU to run the job in. Make sure to check the gates used in the quantum circuit are compatible with that device.
