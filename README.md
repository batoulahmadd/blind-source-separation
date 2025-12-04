# blind-source-separation

## Project Overview

This project implements Blind Source Separation (BSS) techniques to separate mixed communication signals from interference and noise.

Two major stages were developed:

1- MATLAB Simulation:

Comparing three ICA algorithms for accuracy + computation time:
  - FastICA
  - SOBI
  - AMUSE

2- Practical SDR Implementation:
The selected algorithm (SOBI) was implemented on GNU Radio using Python blocks, with real-time signal separation tests:
  - Clean signals
  - Noisy signals (with / without filtering)
  - Single-signal vs noise separation


## Motivation
Blind source separation plays an important role in communication systems, biomedical signal processing, and audio enhancement. This work demonstrates the efficiency of ICA techniques in extracting independent components from mixed observations, even in noisy environments.

## Technologies Used
| Category | Details |
|---------|---------|
| Programming Language | MATLAB |
| Algorithms | FastICA, AMUSE, SOBI |
| Evaluation Metrics | Correlation coefficient, Execution time |
| Experiment Type | Synthetic signals with varying SNR levels |

## Key Results
| Algorithm | Correlation Performance | Execution Time | Best Use Case |
|----------|------------------------|----------------|---------------|
| **SOBI** | ⭐ Highest accuracy (best correlation across SNRs) | Medium | Signals with temporal structure |
| **AMUSE** | Good accuracy | ⭐ Fastest | Low-latency processing |
| **FastICA** | Moderate accuracy | Slowest | Highly independent sources |

- SOBI provides the **best signal recovery**
- AMUSE is the **fastest algorithm**  
- ICA works reliably even with **added Gaussian noise**


## Repository Structure

bss-ica-gnu-radio/
│
├── matlab_simulation/
│   ├── scripts/
│   │   ├── fastica_sim.m
│   │   ├── sobi_sim.m
│   │   ├── amuse_sim.m
│   │   └── utils/ (correlation, SNR, mixing, etc.)
│   ├── plots/ (correlation graphs, timing charts)
│   └── README.md  <-- short note for running MATLAB scripts
│
├── gnu_radio/
│   ├── experiment_1_clean_signals/
│   │   ├── sobi_block.py
│   │   ├── correlation_block.py
│   │   └── flowgraph.png
│   │
│   ├── experiment_2_noisy_signals/
│   │   ├── without_filter/
│   │   │   ├── sobi_block.py
│   │   │   ├── correlation_block.py
│   │   │   └── flowgraph.png
│   │   └── with_filter/
│   │       ├── sobi_block.py
│   │       ├── correlation_block.py
│   │       └── flowgraph.png
│   │
│   ├── experiment_3_signal_vs_noise/
│   │   ├── without_filter/
│   │   │   ├── sobi_block.py
│   │   │   └── flowgraph.png
│   │   └── with_filter/
│   │       ├── sobi_block.py
│   │       └── flowgraph.png
│
│
├── LICENSE
└── README.md

## How to Run MATLAB Simulations
Open MATLAB → run desired script in matlab_simulation/scripts/

Results are saved automatically (correlation + execution time).

## How to Run GNU Radio Experiments
Open .grc files in GNU Radio Companion

Connect Python blocks

Execute flowgraph

Correlation values will appear inside GNU Radio as live output.


## Publications & Documentation

Full research report (Arabic PDF):
🔗 https://drive.google.com/file/d/1_SUB_-elds2g_sSnEgMkjOzltCwgjL4D/view


## Skills Demonstrated

- Blind Source Separation (ICA methods)
- DSP for wireless systems
- GNU Radio SDR implementation
- Experimental design + performance analytics
- Python + MATLAB development


## Author

Albatoul Ahmad

Telecom Engineer | Wireless & Signal Processing

📩 batoulahmad292@gmail.com
