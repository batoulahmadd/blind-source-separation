# blind-source-separation

## Project Overview
This project focuses on Blind Source Separation (BSS) for communication systems. Three Independent Component Analysis (ICA) algorithms — SOBI, FastICA, and AMUSE — were evaluated using MATLAB to determine the most efficient approach for separating mixed radio signals corrupted by noise and interference.

Based on benchmarking results (correlation & execution time), SOBI was selected for real-time implementation on GNU Radio using Python, demonstrating practical applicability for Software-Defined Radio (SDR) systems.


## Motivation
In modern wireless communications, received signals are often corrupted or mixed with other signals. BSS provides a robust solution without prior knowledge of the sources or the mixing process. This improves:
  - Communication quality
  - Signal recovery accuracy
  - Error rate reduction
  - SDR flexibility for dynamic environments


## Methods & System Architecture

1- Simulation in MATLAB:
  - Generated synthetic signal mixtures (sin, square, sawtooth)
  - Two scenarios tested:
      - Clean signals
      - Noisy signals under varying SNR ∈ [-20 dB → 20 dB]
  - Monte-Carlo evaluation across 1000 iterations
  - Algorithms compared using:
      - Correlation coefficient (signal recovery quality)
      - Execution time (computational efficiency)

2- Real Implementation in GNU Radio:
- SOBI executed on live mixed signals
- Real-time correlation calculated
- Three experiments:
      - Clean signals
      - Noisy signals
      - Noise separation (AWGN suppression)


## Technologies Used

| Category | Tools / Frameworks |
|---------|-------------------|
| Simulation | MATLAB (R2021a) |
| Practical SDR Implementation | GNU Radio (Python-based) |
| Programming Languages | MATLAB, Python |
| Hardware Target | Software-Defined Radio (SDR) compatible |



## Key Results Summary

| Scenario | Best Algorithm | Criteria | Notes |
|---------|----------------|---------|------|
| 3 signals, no noise | AMUSE | Fastest execution | Correlation close to 1 |
| 5 signals, no noise | AMUSE | Best accuracy | Very high correlation |
| Low SNR (high noise) | SOBI | Most stable | Robust correlation results |
| Increasing number of sources | SOBI | Most scalable | Minimal performance degradation |
| Practical SDR (GNU Radio) | SOBI | Best overall | Selected for deployment |

With Noise (SNR variation):
- SOBI outperformed others especially at low SNR
- Maintained stable execution time
- Most practical for communication environments

Conclusion: SOBI provides the best balance of accuracy + speed + stability → selected for SDR implementation.


## Repository Structure

blind-source-separation-ica/
│
├── matlab_simulation/
│   ├── fastica/
│   ├── sobi/
│   ├── amuse/
│   └── results/
│
├── gnu_radio_implementation/
│   ├── sobi_block.py
│   ├── correlation_block.py
│   └── flowgraphs/
│
├── docs/
│   ├── block_diagram.png
│   ├── sdr_architecture.png
│   └── algorithm_comparison_plots/
│
├── data/
│   └── sample_signals/
│
├── requirements.txt
└── README.md


## How to Run
SOBI in GNU Radio

1- Install dependencies
    pip install -r requirements.txt

2- Open GNU Radio and load the flowgraph:
    gnu_radio/flowgraphs/sobi_separation.grc

3- Start the flow — view separated signals and correlation results in real time
