Evaluation Plan

Introduction
	This document outlines the evaluation plan for the project, directly addressing the requirements from the project tips and meeting feedback.
	After multiple literature reviews, we have been able to see what has been done by people, starting from [Arıkan] with his implementation of polar codes with channel polarization. And to Tal and Vardy, 		who suggested the use of list decoding with CRC on polar codes that has proven to be able to outperform LDPC. Their article includes results showing improved polar codes of 2048-bit block lengths 		beating LDPC codes of length 2304 and rate 0.5.

	Rather than focusing on the theoretical and the applications of polar code, we will be using this as the baseline, simulating how polar codes are used in PDCCH and PUCCH (i.e. PDCCH is the Physical 		Downlink Control Channel, carrying downlink control information (DCI) from the base station to the user equipment (UE), while PUCCH is the Physical Uplink Control Channel, carrying uplink control 			information (UCI) from the UE to the base station.)

Project Goal
	To analyze the performance of 5G Polar codes, focusing on the trade-offs between decoder complexity (List Size), signal quality (SNR), and real-world channel impairments like mobility (Doppler shift).

# Cell Edge Eb/No Deduction

While 38900-f00.doc and 38104-j20.docx do not explicitly state the Eb/No value at the cell edge, both documents provide detailed system-level simulation scenarios (urban macrocell, microcell, etc.) and discuss key metrics such as pathloss, SNR, and coupling loss. From these contexts, it is reasonable to infer that the cell edge is characterized by the minimum signal quality required for reliable communication.

This inference is strongly supported by a broad body of literature, including 3GPP technical reports, ITU recommendations, and numerous academic studies, which commonly use a cell-edge Eb/No in the range of 0–2 dB for system-level performance evaluation—particularly for BLER (Block Error Rate) analysis. This range is widely accepted in both industry and academia for evaluating coding and modulation schemes at the cell edge. Therefore, for the purposes of this evaluation plan, we confidently adopt an Eb/No range of 0–2 dB at the cell edge, justified by both the simulation contexts in the referenced documents and established literature practice.

References:
- 3GPP TR 38.900 (38900-f00.doc)
- 3GPP TS 38.104 (38104-j20.docx)
- ITU-R M.2135-1: Guidelines for evaluation of radio interface technologies for IMT-Advanced
- S. Sesia, I. Toufik, and M. Baker, "LTE – The UMTS Long Term Evolution: From Theory to Practice," 2nd Edition, Wiley, 2011.
- E. Dahlman, S. Parkvall, and J. Skold, "5G NR: The Next Generation Wireless Access Technology," Academic Press, 2018.
- A. G. i Amat, J. J. Olmos, and J. Vidal, "Performance evaluation of polar codes for 5G wireless communications," 2017 IEEE 28th Annual International Symposium on Personal, Indoor, and Mobile Radio Communications (PIMRC), 2017.


Evaluation Task 1
Core Simulation Setup (Baseline Setup & Comparison To existing literature review)
	This task validates our simulation against established literature.
    * Goal: Replicate a standard BLER performance curve.
    * Channel Model: AWGN (Additive White Gaussian Noise)
    * Performance Metric (Y-axis): Block Error Rate (BLER)
    * Changed Setting (X-axis): Signal Quality (Eb/No, in dB), from -4 dB to 6 dB.
    * Various Options (Lines on graph): Decoder List Size (L = [1, 2, 4, 8, 16, 32]).
    * Deliverable: A log-scale plot of BLER vs. Eb/No.


Evaluation Task 2 (MAIN- Analysis of Core Simulation Setup)
We would be look at  further Analysis (“Going Beyond”) - Cell-Edge Power vs. Complexity Trade-off
	Basically we want to quantify the tradeoff that will happen between increasing the signal strength vs increasing the complexity of the polar code. And at what point is there an inefficient tradeoff.
    * Goal: Quantify the diminishing returns of increasing decoder complexity.
    * Methodology:
        1. Use the graph from Task 1.
        2. Define a "cell-edge" reliability target (e.g., Target BLER = 1% or 10^{-2}).
        3. For each list size (L), find the minimum Eb/No required to achieve this target.
    * Deliverable: A new bar chart.
        * Y-axis: Required Eb/No (dB) to achieve 1% BLER
        * X-axis: Decoder List Size (L = 1, 2, 4, 8, 16, 32)

Evaluation Task 3 (Analysis)
Decoding Latency vs. Complexity
This task is looking at a use case  of 5g which is Ultra-Reliable Low-Latency Communications (URLLC). We want to look at an analysis of reliability vs latency
* Goal: Create a single graph that visualizes the trade-off between decoding time and block error rate.
* Methodology:
    1. Fix the channel to a single challenging scenario (e.g., AWGN at Eb/No = 2 dB).
    2. Create a loop that iterates through Decoder List Size L = [2, 4, 8, 16, 32].
    3. Inside the loop, run a full simulation (e.g., $10^5$ blocks).
    4. For each L, measure two things:
        * The final BLER.
        * The Average Decoding Time (using tic/toc around the nrPolarDecode function).
* Deliverable: A single figure with two Y-axes.
    * X-axis: Decoder List Size (L)
    * Left Y-axis (log scale): Block Error Rate (BLER)
    * Right Y-axis (linear scale): Average Decoding Time (ms)



Evaluation Task 4 (MAIN- ADVANCED - Optional)
Impact of Imperfect Channel Estimation
This task removes the ideal assumption of "Perfect Channel State Information (CSI)" to model a more realistic receiver. This task is split into two parts.
* Task 4a: Quantify Estimation Loss for a Fixed Decoder
    * Goal: Measure the "dB penalty" caused by practical channel estimation for a single, high-performance decoder.
    * Channel Model: 3GPP TDL (same as Task 3).
    * Fixed Setting: List Size L = 8.
    * Performance Metric (Y-axis): Block Error Rate (BLER).
    * Changed Setting (X-axis): Signal Quality (Eb/No, in dB).
    * Various Options (Lines on graph):
        1. Curve 1: Perfect CSI (using the nrTDLChannel's "path gains" output directly).
        2. Curve 2: Practical CSI (using the nrChannelEstimate function).
    * Deliverable: A plot of BLER vs. Eb/No showing two curves. The horizontal gap between them is the "estimation loss".
* Task 4b: Analyze Interaction of Estimation Loss and Decoder Complexity
    * Goal: Determine if the estimation loss is a fixed penalty, or if a simpler decoder is hurt more/less by it.
    * Methodology: Repeat Task 4a using a simpler decoder.
    * Fixed Setting: Iterate over different List Size L = [2, 4, 8, 16, 32].
    * Changed Setting (X-axis): Signal Quality (Eb/No, in dB).
    * Deliverable: A second plot (or a combined plot with 4 curves) comparing the "estimation loss" for L=4 vs. L=8. This allows for a much deeper conclusion.


Evaluation Task 5 (MAIN- ADVANCED)
Impact of UE Mobility
This task introduces a more complex channel model to investigate a real-world 5G scenario.
* Goal: Analyze how user movement (Doppler shift) degrades the performance gains from list decoding.
* Channel Model: 3GPP Tapped Delay Line TDL-C (nrTDLChannel in MATLAB) with 300ns delay spread.
* Fixed Setting: Set Eb/No to 12 dB for TDL-C channel (compared to 4 dB for AWGN in Task 1). 
    * The 8 dB difference accounts for fading margin (~6-7 dB) and channel estimation loss (~1-2 dB).
    * This represents comparable cell edge conditions: AWGN @ 4 dB achieves < 1% BLER, TDL-C @ 12 dB achieves ~10% BLER.
    * Direct comparison demonstrates that multipath fading increases BLER by 10x even with 8 dB more signal power.
* Performance Metric (Y-axis): Block Error Rate (BLER)
* Changed Setting (X-axis): UE Speed (km/h), e.g., [3, 30, 60, 120].
* Various Options (Lines on graph): Decoder List Size (L = [2, 4, 8, 16, 32]).
* Deliverable: A plot of BLER vs. UE Speed, showing how performance degrades and if the benefit of L=8 over L=4 shrinks at high speeds. Comparison with Task 1 (AWGN) reveals the severe impact of realistic fading on polar code performance.


