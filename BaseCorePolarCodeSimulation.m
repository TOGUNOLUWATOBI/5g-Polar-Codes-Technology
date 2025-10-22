% A Comprehensive 5G NR Polar Code Simulation in MATLAB
% This script models a baseline PDCCH link over an AWGN channel
% to evaluate Block Error Rate (BLER) performance
%
% MODIFIED TO SWEEP MULTIPLE SCL LIST SIZES (L)

% Clear workspace and command window, seed the random number generator for repeatability
clear; clc;
rng('default');

% =========================================================================
% SECTION 1: MASTER SIMULATION PARAMETERS
% =========================================================================

% -- Code Configuration
K = 54;         % Message length in bits, including CRC (e.g., 30 info + 24 CRC)
E = 124;        % Rate-matched output length in bits

% -- List Sizes to Sweep
L_vec = [2, 4, 8, 16, 32]; % Vector of SCL decoder list sizes to test

% -- Link Configuration
linkDir = 'DL'; % 'DL' for Downlink (PDCCH) or 'UL' for Uplink (PUCCH)
modulation = 'QPSK'; % Modulation scheme for control channels
EbNo = -4:0.5:6;     % Eb/No range in dB to sweep for BLER curve

% -- Simulation Control
maxNumFrames = 10000;   % Maximum number of frames to simulate at each Eb/No point
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
% SECTION 3: SIMULATION PROCESSING LOOP
% =========================================================================

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
            % 1. Generate random message bits
            msg = randi([0 1], K-crcLen, 1);
            
            % 2. Attach CRC
            msg_crc = nrCRCEncode(msg, crcPoly);
            
            % 3. Polar Encode and Rate Match
            rateMatchedBits = nrPolarEncode(msg_crc, E, nMax, iIL);
            
            % 4. Modulate
            modulatedSymbols = nrSymbolModulate(rateMatchedBits, modulation);
            
            % -----------------------------------------------------------------
            % CHANNEL
            % -----------------------------------------------------------------
            % 6. Pass through AWGN channel
            receivedSymbols = channel(modulatedSymbols);
            
            % -----------------------------------------------------------------
            % RECEIVER
            % -----------------------------------------------------------------
            % 7. Soft Demodulate to get LLRs
            rxLLRs = nrSymbolDemodulate(receivedSymbols, modulation, noiseVar);
            
            % 8. Polar Decode (CA-SCL decoder, includes rate recovery)
            % *** MODIFICATION: Use the current loop variable 'L' ***
            decodedBits = nrPolarDecode(rxLLRs, K, E, L, nMax, iIL, crcLen);
            
            % -----------------------------------------------------------------
            % ERROR CALCULATION
            % -----------------------------------------------------------------
            % 10. Check for a block error
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
% SECTION 4: PLOT RESULTS
% =========================================================================
figure;
hold on; % Allow multiple lines to be plotted on the same axes

% Define markers for the different plots
markers = {'-o', '-s', '-^', '-d', '-v', '-x'};

% Loop through each list size result and plot it
for i_L = 1:length(L_vec)
    % Select marker, wrap around if we have more lines than markers
    marker_style = markers{mod(i_L-1, length(markers)) + 1};
    
    semilogy(EbNo, bler_results(i_L, :), marker_style, 'LineWidth', 1.5, 'MarkerSize', 6);
end

hold off;
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate (BLER)');

% Add the dynamic legend
legend(legend_entries, 'Location', 'SouthWest');

% Add a title to describe the simulation parameters
title(sprintf('Polar Code BLER vs. Eb/No (K=%d, E=%d, %s, %s)', K, E, linkDir, modulation));

ylim([1e-5 1]); % Set Y-axis limits