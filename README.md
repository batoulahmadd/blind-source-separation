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
