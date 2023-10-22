from braket.aws import AwsQuantumJob

job = AwsQuantumJob.create(
    "arn:aws:braket:::device/quantum-simulator/amazon/sv1",
    source_module="algorithm_script.py",
    entry_point="algorithm_script:start_here",
    wait_until_complete=True
)