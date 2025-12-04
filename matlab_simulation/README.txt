
Overview:

	This folder contains MATLAB simulations to evaluate and compare three ICA algorithms: SOBI, FastICA, and AMUSE. The comparison is based on execution time and correlation between original and separated signals.

	The simulations help determine the most efficient algorithm, which is later implemented in GNU Radio using Python for SDR-based communication.


Scripts:

	1- Demo scripts (main experiments):

		demo_three_signals_conversion.m – Three signals using conversion methods

		demo_five_signals_var_noise.m – Five signals, with variable noise conditions

		demo_single_signal_noisy.m – Single signal with noise

	2- Algorithm scripts:

		fastICA.m – FastICA algorithm

		sobi.m – SOBI algorithm

		amuse.m – AMUSE algorithm

	3- Utility scripts

		whitened_rows.m – Whitening function for preprocessing

		center_rows.m – Centering function for preprocessing


How to Run:

	1- Select the demo script depending on your scenario:

		- Three signals → demo_three_signals_conversion.m

		- Five signals with/without noise → demo_five_signals_var_noise.m

		- Single noisy signal → demo_single_signal_noisy.m

	2- Make sure all scripts (algorithms and utilities) are in the MATLAB path.

	3- Run the script. It will automatically:

		Generate signals (clean or noisy)

		Apply the three ICA algorithms

		Compute correlation with original signals

		Measure execution time

		Save plots in the plots/ folder

Outputs:

	Plots:

		Correlation vs. SNR

		Execution time vs. SNR