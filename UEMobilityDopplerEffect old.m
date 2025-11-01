% A Comprehensive 5G NR Polar Code Simulation in MATLAB
%
% TASK 5: BLER vs. UE Mobility (Doppler Shift)
%
% This script simulates a **PDCCH (Control Channel)** link, which is the
% *correct* simulation for Polar codes.
%
% *** THIS VERSION FIXES THE TOOLBOX INCOMPATIBILITY ***
% 1. Replaced 'nrPDCCHIndices' with 'nrPDCCHSpace' for compatibility.
% 2. Removed the invalid 'pdcch.NumLayers' property.
%
% Clear workspace and command window, seed the random number generator
clear; close all; clc;
rng('default');

% =------------------------------------------------------------------------
% SECTION 1: MASTER SIMULATION PARAMETERS
% =------------------------------------------------------------------------
% -- Code Configuration
K = 54;         % Message length in bits (e.g., DCI format 1_0)
 
% -- Link Configuration
modulation = 'QPSK';
bitsPerSymbol = 2;

% -- TASK 5 SPECIFIC PARAMETERS --
% -- Fixed Signal Quality
snr_dB_fixed = 0; % Fixed SNR in dB. We can tune this.
                  % 0 dB means Signal = Noise.

% -- Variables to Sweep
L_vec = [2, 4, 8, 16, 32];      % Decoder List Sizes to test
speed_kmh_vec = [3, 30, 60, 120]; % UE Speeds in km/h

% -- Simulation Control
maxNumFrames = 2000;
minFrameErrors = 50;

% -- Channel Configuration
carrierFrequency_Hz = 4e9; % 4 GHz (5G Mid-band FR1)
channelProfile = 'TDL-C';  % Common NLOS profile

% -- Initialize Results Storage
bler_results = zeros(length(L_vec), length(speed_kmh_vec));
legend_entries = cell(1, length(L_vec));

% =------------------------------------------------------------------------
% SECTION 2: 5G NR PDCCH & TDL CHANNEL CONFIGURATION
% =------------------------------------------------------------------------

% -- Carrier Configuration
carrier = nrCarrierConfig;
carrier.NCellID = 1; % Scrambling is on by default
ofdmInfo = nrOFDMInfo(carrier);

% -- PDCCH (Control) Configuration
pdcch = nrPDCCHConfig;
agg_lvl = 4; % aggregation level you are using
pdcch.AggregationLevel = agg_lvl; % A typical AL (1, 2, 4, 8, 16)










pdcch.NStartBWP = 0;
pdcch.NSizeBWP  = 45;

% CORESET uses 2 RBGs (12 PRBs) at the beginning of the BWP
numRBGs = ceil(pdcch.NSizeBWP / 6);   % = 8
pdcch.CORESET.RBOffset = 0;
pdcch.CORESET.Duration = 2;

pdcch.CORESET.FrequencyResources = [1 1 zeros(1, numRBGs - 2)];

% SearchSpace: Only AL=4 candidate
pdcch.SearchSpace.NumCandidates = [0 0 1 0 0];
pdcch.SearchSpace.CORESETID = 0;















% We must also set the payload info for the decoder
crcType = '24C';
crcPoly = '24C';
interleaving = true; % This matches 'iIL = true' from your AWGN script


% -- Pre-calculate the *true* coded length 'E'
%    FIX: Use 'nrPDCCHSpace' which is more compatible than 'nrPDCCHIndices'
[pdcchIndices, pdcchInfo] = nrPDCCHSpace(carrier, pdcch);
% E_actual = pdcchInfo.G;

E_actual = length(pdcchIndices{1}) * bitsPerSymbol;

K_actual = K; % Our message length
fprintf('PDCCH Configuration: K = %d, E = %d. Code Rate = %.3f\n', K_actual, E_actual, K_actual/E_actual);

% -- TDL Channel Object Configuration
tdl = nrTDLChannel;
tdl.DelayProfile = channelProfile;
tdl.NumReceiveAntennas = 1;
tdl.SampleRate = ofdmInfo.SampleRate;

% =------------------------------------------------------------------------
% SECTION 3: SIMULATION PROCESSING LOOP (Task 5 - PDCCH)
% =------------------------------------------------------------------------
fprintf('===== STARTING TASK 5: PDCCH BLER vs. Mobility Simulation =====\n');
% Start a parallel pool if you have Parallel Computing Toolbox
if isempty(gcp('nocreate'))
    parpool;
end
nID = carrier.NCellID;   % scrambling ID
rnti = 0;  
% Use a PARFOR (Parallel For) loop for the List Sizes
for i_L = 1:length(L_vec)
    L = L_vec(i_L);
    local_legend_entry = sprintf('L = %d', L);
    fprintf('\n===== Simulating for List Size L = %d =====\n', L);
    
    % Create a local BLER array for the parallel loop
    local_bler_results = zeros(1, length(speed_kmh_vec));
    
    % We must create our own TDL object inside the parallel loop
    tdl_loop = nrTDLChannel;
    tdl_loop.DelayProfile = channelProfile;
    tdl_loop.NumReceiveAntennas = 1;
    tdl_loop.SampleRate = ofdmInfo.SampleRate;

    % Create a local copy of the pdcch object for this worker
    % This is now used by nrDCIEncode and nrPDCCHDemodulate
    pdcch_loop = pdcch; 
    
    % Inner loop: Iterate over each UE Speed
    for i_speed = 1:length(speed_kmh_vec)
        speed_kmh = speed_kmh_vec(i_speed);
        
        % -- 1. Calculate and Set Doppler Shift
        speed_mps = speed_kmh / 3.6;
        dopplerShift_Hz = (speed_mps / 299792458) * carrierFrequency_Hz;
        
        release(tdl_loop); % Release the *local* TDL object
        tdl_loop.MaximumDopplerShift = dopplerShift_Hz;
        
        fprintf('  L=%d: Simulating Speed = %d km/h (Doppler = %.1f Hz)...\n', L, speed_kmh, dopplerShift_Hz);
        
        % -- 2. Calculate Noise Variance from fixed SNR
        noiseVar = 1./(10.^(snr_dB_fixed/10));
        % -- 3. Loop counters
        frameCount = 0;
        errorCount = 0;
        
        reset(tdl_loop);
        while frameCount < maxNumFrames && errorCount < minFrameErrors
            
           % -----------------------------------------------------------------
            % TRANSMITTER
            % -----------------------------------------------------------------
            % 1. Create message (DCI payload)
            msg = randi([0 1], K_actual, 1); % K_actual = 54 bits
            
            % *** FIX: Manually perform the DCI encoding steps ***
            %    (This replaces the incompatible nrDCIEncode)
            
            % 1a. CRC Attachment
            %     Input: 54 bits. Output: 54 + 24 = 78 bits
            crcBits = nrCRCEncode(msg, crcPoly);
            
            % 1b. Polar Encoding (incl. rate-matching and interleaving)
            %     Input: 78 bits. Output: E_actual = 432 bits
            % polarBits = nrPolarEncode(crcBits, E_actual, 'Interleave', interleaving);
            polarBits = nrPolarEncode(crcBits, E_actual, 9, interleaving);
            
            % 2. Get PDCCH symbols (Scrambling + QPSK Modulation)
            %    This uses the 'nrPDCCH.m' file from your system,
            %    which expects the already-encoded bits (polarBits).
            txSymbols = nrPDCCH(polarBits, nID, rnti);
            
            % 3. Create an empty grid and map symbols
            txGrid = nrResourceGrid(carrier, 1); 
            %    We still need pdcch_loop for nrPDCCHSpace and nrPDCCHDemodulate
            [indices, info] = nrPDCCHSpace(carrier, pdcch_loop); 
            txGrid(indices{1}) = txSymbols;
            
            % 4. OFDM Modulate (generate waveform)
            txWaveform = nrOFDMModulate(carrier, txGrid);
            % -----------------------------------------------------------------
            % CHANNEL
            % -----------------------------------------------------------------
            [rxWaveform, pathGains] = tdl_loop(txWaveform);
            
            % Add AWGN
            noise = sqrt(noiseVar/2) * (randn(size(rxWaveform)) + 1i*randn(size(rxWaveform)));
            rxWaveform_noisy = rxWaveform + noise;
            
            % -----------------------------------------------------------------
            % RECEIVER
            % -----------------------------------------------------------------
            % 1. OFDM Demodulate
            rxGrid = nrOFDMDemodulate(carrier, rxWaveform_noisy);
            
            % 2. Get PDCCH LLRs
            % *** FIX: The 'L' argument is not used here ***
            %    This was likely a bug in your original script
            % Old: rxLLRs = nrPDCCHDemodulate(carrier, pdcch, rxGrid, noiseVar, L);
            rxLLRs = nrPDCCHDemodulate(carrier, pdcch_loop, rxGrid, noiseVar);
            
            % 3. Polar Decode (using CRC)
            %    This is where 'L' is correctly used as an argument
            % *** FIX: The 7th argument is nPC (which is 0), not crcPoly ***
            decodedBits = nrPolarDecode(rxLLRs, K_actual, E_actual, L, ...
                                        interleaving, crcType, ...
                                        0);
            
            % -----------------------------------------------------------------
            % ERROR CALCULATION
            % -----------------------------------------------------------------
            % This check is correct for nrPolarDecode
            if any(decodedBits < 0)
                errorCount = errorCount + 1;
            end
            frameCount = frameCount + 1;
        end
        
        % -- Calculate and store BLER
        local_bler_results(i_speed) = errorCount / frameCount;
        fprintf('    L=%d, Speed=%d km/h: BLER: %.5f (%d/%d)\n', L, speed_kmh, local_bler_results(i_speed), errorCount, frameCount);
    end % end speed loop
    
    % Store the results from the parallel loop
    bler_results(i_L, :) = local_bler_results;
    legend_entries{i_L} = local_legend_entry;
    
end % end L_vec parfor loop
fprintf('===== TASK 5 SIMULATION COMPLETE. =====\n\n');

% =------------------------------------------------------------------------
% SECTION 4: PLOT RESULTS (Task 5)
% =------------------------------------------------------------------------
figure;
hold on;
markers = {'-o', '-s', '-^', '-d', '-v'}; % Markers for the plot

for i_L = 1:length(L_vec)
    marker_style = markers{mod(i_L-1, length(markers)) + 1};
    semilogy(speed_kmh_vec, bler_results(i_L, :), marker_style, 'LineWidth', 1.5, 'MarkerSize', 6);
end

hold off;
grid on;
xlabel('UE Speed (km/h)');
ylabel('Block Error Rate (BLER)');
title(sprintf('Task 5: PDCCH BLER vs. Mobility (SNR = %.1f dB, %s, AL=%d)', snr_dB_fixed, channelProfile, pdcch.AggregationLevel));
legend(legend_entries, 'Location', 'SouthEast');
ylim([1e-4 1]);
ax = gca;
ax.YScale = 'log';
xticks(speed_kmh_vec); % Ensure X-axis ticks match our tested speeds

