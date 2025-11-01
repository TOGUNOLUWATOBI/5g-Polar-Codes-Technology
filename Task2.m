% A Comprehensive 5G NR Polar Code Simulation in MATLAB
% This script models a baseline PDCCH link over an AWGN channel
% to evaluate Block Error Rate (BLER) performance.
%
% This script runs BOTH Task 1 (BLER vs. Eb/No) and Task 2 (Analysis)
%
% Clear workspace and command window, seed the random number generator
clear; close all; clc;
rng('default');

% =========================================================================
% SECTION 1: MASTER SIMULATION PARAMETERS
% =========================================================================
% -- Code Configuration
K = 54;         % Message length in bits, including CRC (e.g., 30 info + 24 CRC)
E = 124;        % Rate-matched output length in bits
% -- List Sizes to Sweep
L_vec = [1, 2, 4, 8, 16, 32]; % Vector of SCL decoder list sizes to test

% -- Link Configuration
linkDir = 'DL'; % 'DL' for Downlink (PDCCH) or 'UL' for Uplink (PUCCH)
modulation = 'QPSK'; % Modulation scheme for control channels
EbNo = -3:0.5:5;     % Eb/No range in dB to sweep for BLER curve

% -- Simulation Control
maxNumFrames = 5000;   % Maximum number of frames to simulate at each Eb/No point
minFrameErrors = 100;   % Minimum number of frame errors to collect for statistical significance

% -- Initialize Results Storage
bler_results = zeros(length(L_vec), length(EbNo)); % 2D array for BLER
legend_entries = cell(1, length(L_vec)); % For storing plot legend strings

% =========================================================================
% SECTION 2: DERIVED PARAMETERS (based on Link Direction)
% =========================================================================
% These parameters are set according to the 3GPP standard for Polar codes
if strcmpi(linkDir,'DL')
    % Downlink (PDCCH) specific parameters
    crcPoly = '24C';      % CRC polynomial for DCI
    crcLen = 24;          % CRC length for DCI
    nMax = 9;             % Max value of n for mother code length 2^n
    iIL = true;           % Input interleaving is enabled for downlink
    iBIL = false;         % Coded-bit interleaving is disabled for downlink
else % Uplink
    % Uplink (PUCCH) specific parameters
    crcPoly = '11';       % CRC polynomial for UCI (for K > 30)
    crcLen = 11;          % CRC length for UCI
    nMax = 10;            % Max value of n for mother code length 2^n
    iIL = false;          % Input interleaving is disabled for uplink
    iBIL = true;          % Coded-bit interleaving is enabled for uplink
end

% =========================================================================
% SECTION 3: SIMULATION PROCESSING LOOP (Task 1)
% =========================================================================
fprintf('===== STARTING TASK 1: BLER vs. Eb/No Simulation =====\n\n');

% Outer loop to iterate over each List Size
for i_L = 1:length(L_vec)
    
    L = L_vec(i_L); % Get current list size
    legend_entries{i_L} = sprintf('L = %d', L); % Create legend string
    
    fprintf('===== Simulating for List Size L = %d =====\n', L);
    
    % Inner loop to iterate over each Eb/No point
    for i_ebno = 1:length(EbNo)
        
        % -- Calculate Noise Variance from Eb/No
        codeRate = K/E;
        bitsPerSymbol = 2; % For QPSK
        EsNo = EbNo(i_ebno) + 10*log10(bitsPerSymbol);
        snr_dB = EsNo + 10*log10(codeRate);
        noiseVar = 1./(10.^(snr_dB/10));
        
        % -- Instantiate AWGN Channel with the calculated noise variance
        channel = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);
        
        % -- Loop counters for each Eb/No point
        frameCount = 0;
        errorCount = 0;
        
        fprintf('  Simulating Eb/No = %.1f dB... \n', EbNo(i_ebno));
        
        while frameCount < maxNumFrames && errorCount < minFrameErrors
            
            % -----------------------------------------------------------------
            % TRANSMITTER
            % -----------------------------------------------------------------
            msg = randi([0 1], K-crcLen, 1);
            msg_crc = nrCRCEncode(msg, crcPoly);
            rateMatchedBits = nrPolarEncode(msg_crc, E, nMax, iIL);
            modulatedSymbols = nrSymbolModulate(rateMatchedBits, modulation);
            
            % -----------------------------------------------------------------
            % CHANNEL
            % -----------------------------------------------------------------
            receivedSymbols = channel(modulatedSymbols);
            
            % -----------------------------------------------------------------
            % RECEIVER
            % -----------------------------------------------------------------
            rxLLRs = nrSymbolDemodulate(receivedSymbols, modulation, noiseVar);
            decodedBits = nrPolarDecode(rxLLRs, K, E, L, nMax, iIL, crcLen);
            
            % -----------------------------------------------------------------
            % ERROR CALCULATION
            % -----------------------------------------------------------------
            if any(decodedBits ~= msg_crc)
                errorCount = errorCount + 1;
            end
            
            frameCount = frameCount + 1;
        end
        
        % -- Calculate and store BLER for this Eb/No point and L value
        bler_results(i_L, i_ebno) = errorCount / frameCount;
        fprintf('    BLER: %.5f (%d/%d)\n', bler_results(i_L, i_ebno), errorCount, frameCount);
        
    end % end EbNo loop
    
    fprintf('\n'); % Add space between L-value simulations
    
end % end L_vec loop


% =========================================================================
% SECTION 4: TASK 2 ANALYSIS & PLOT (Cell-Edge Trade-off)
% =========================================================================
fprintf('===== STARTING TASK 2: Power vs. Complexity Analysis =====\n');

% -- Analysis Parameters
target_bler = 0.01; % 1% BLER reliability target
required_ebno = zeros(1, length(L_vec)); % Array to store results

% -- Analysis Loop
for i_L = 1:length(L_vec)
    % Find the first index (lowest Eb/No) where BLER is <= target_bler
    % This finds the *first* point where the line crosses *below* the target
    idx = find(bler_results(i_L, :) <= target_bler, 1, 'first');
    
    if ~isempty(idx)
        % Found an intersection.
        % A simple approximation is to just take the Eb/No value
        % A more precise method (linear interpolation) is also possible
        required_ebno(i_L) = EbNo(idx);
        fprintf('  L = %d: Target BLER (%.2f) reached at %.1f dB Eb/No\n', L_vec(i_L), target_bler, required_ebno(i_L));
    else
        % Never reached the target in our simulated Eb/No range
        required_ebno(i_L) = NaN; % Store 'Not a Number'
        fprintf('  L = %d: Never reached target BLER (%.2f) in simulated range.\n', L_vec(i_L), target_bler);
    end
end

% --- Plot the Task 2 Bar Chart ---
figure('Name', 'Task 2: Power vs. Complexity Trade-off (Cell-Edge)'); % Create Figure 2 (a NEW figure window)
% Create categorical labels for the X-axis
string_L_vec = string(L_vec); % Convert numbers to strings
% Create an ORDINAL categorical array to maintain the correct sort order
x_categories = categorical(string_L_vec, string_L_vec, 'Ordinal', true);

% Replace NaN values with 0 for plotting, so the bar is just 'missing'
plot_data = required_ebno;
plot_data(isnan(plot_data)) = 0;

b = bar(x_categories, plot_data);
grid on;
xlabel('Decoder List Size (L)');
ylabel(sprintf('Required Eb/No (dB) for %.2f BLER', target_bler));
title('Task 2: Power vs. Complexity Trade-off (Cell-Edge)');

% Add text labels on top of each bar
xtips = b.XEndPoints;
ytips = b.YEndPoints;
% Create labels, replacing 0 (from NaN) with 'N/A'
labels = strings(1, length(ytips));
for i = 1:length(ytips)
    if ytips(i) == 0
        labels(i) = 'N/A';
    else
        labels(i) = string(round(ytips(i), 2));
    end
end

text(xtips, ytips, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 10, 'FontWeight', 'bold');

ylim([min(EbNo)-2 max(EbNo)+2]); % Set Y-axis limits dynamically

fprintf('===== TASK 2 ANALYSIS COMPLETE. =====\n');
