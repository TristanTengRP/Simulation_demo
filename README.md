# Simulation_demo
This repository contains simulation scripts for validating the Linear Weighted Combination (LWC) and Two-Threshold Combined (TTC) schemes. It includes estimated feature datasets for antenna-specific Mutual Coupling (MC) and spatial Angle of Arrival (AoA) to ensure the reproducibility of the authentication performance.
The provided scripts and datasets allow for the reproduction of the theoretical and simulation results concerning the Linear Weighted Combination (LWC) and Two-Threshold Combined (TTC) authentication schemes.

**File Structure:**

1) **LWC_Model_Validation.m:** Validates the theoretical detection and false alarm probabilities for the LWC scheme.

2) **TTC_Model_Validation.m:** Verifies the authentication performance of the TTC-G schemes across various SNR regimes.

3) **Fingerprint Datasets (.mat files):** The repository includes 5 .mat files containing estimated feature data used in our performance evaluation.

**System Requirements**

1) Software: MATLAB R2021b or later.

2) Toolboxes: Signal Processing Toolbox, Optimization Toolbox (for SQP weight allocation).

**Usage**

1) Clone this repository to your local machine. 

2) Ensure all .mat files are in the same directory as the .m scripts. 

3) Run LWC_Model_Validation.m or TTC_Model_Validation.m to generate the performance curves.
