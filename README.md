# Simulation_demo
This repository contains simulation scripts for validating the Linear Weighted Combination (LWC) and Two-Threshold Combined (TTC) schemes. It includes estimated feature datasets for antenna-specific Mutual Coupling (MC) and spatial Angle of Arrival (AoA) to ensure the reproducibility of the authentication performance.
The provided scripts and datasets allow for the reproduction of the theoretical and simulation results concerning the Linear Weighted Combination (LWC) and Two-Threshold Combined (TTC) authentication schemes.

**File Structure:**

**LWC_Model_Validation.m:** Validates the theoretical detection and false alarm probabilities for the LWC scheme.
**TTC_Model_Validation.m:** Verifies the authentication performance of the TTC-G schemes across various SNR regimes.

Fingerprint Datasets (.mat files): The repository includes 5 .mat files containing estimated feature data used in our performance evaluation.

**System Requirements**
Software: MATLAB R2021b or later.
Toolboxes: Signal Processing Toolbox, Optimization Toolbox (for SQP weight allocation).

**Usage**
Clone this repository to your local machine. 
Ensure all .mat files are in the same directory as the .m scripts. 
Run LWC_Model_Validation.m or TTC_Model_Validation.m to generate the performance curves.
