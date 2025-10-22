% A Comprehensive 5G NR Polar Code Simulation in MATLAB
% This script models a baseline PDCCH link over an AWGN channel
% to evaluate Block Error Rate (BLER) performance.

% Clear workspace and command window, seed the random number generator for repeatability
clear; clc;
rng('default');

% =========================================================================
% SECTION 1: MASTER SIMULATION PARAMETERS
% =========================================================================

% -- Code Configuration
K = 54;         % Message length in bits, including CRC (e.g., 30 info + 24 CRC) [6, 8]
E = 124;        % Rate-matched output length in bits [6, 8]
L = 8;          % SCL decoder list size, a power of two [1, 2, 3, 4][6, 8]

% -- Link Configuration
linkDir = 'DL'; % 'DL' for Downlink (PDCCH) or 'UL' for Uplink (PUCCH)
modulation = 'QPSK'; % Modulation scheme for control channels [6]
EbNo = -4:0.5:6;     % Eb/No range in dB to sweep for BLER curve

% -- Simulation Control
maxNumFrames = 10000;   % Maximum number of frames to simulate at each Eb/No point
minFrameErrors = 100;   % Minimum number of frame errors to collect for statistical significance [6]

% -- Initialize Results Storage
bler = zeros(size(EbNo)); % Array to store Block Error Rate results

% =========================================================================
% SECTION 2: DERIVED PARAMETERS (based on Link Direction)
% =========================================================================
% These parameters are set according to the 3GPP standard for Polar codes [7]

if strcmpi(linkDir,'DL')
    % Downlink (PDCCH) specific parameters
    crcPoly = '24C';      % CRC polynomial for DCI [6]
    crcLen = 24;          % CRC length for DCI [6]
    nMax = 9;             % Max value of n for mother code length 2^n [7]
    iIL = true;           % Input interleaving is enabled for downlink [6, 7]
    iBIL = false;         % Coded-bit interleaving is disabled for downlink [8]
else % Uplink
    % Uplink (PUCCH) specific parameters
    crcPoly = '11';       % CRC polynomial for UCI (for K > 30) [8]
    crcLen = 11;          % CRC length for UCI [8]
    nMax = 10;            % Max value of n for mother code length 2^n [7]
    iIL = false;          % Input interleaving is disabled for uplink [7, 8]
    iBIL = true;          % Coded-bit interleaving is enabled for uplink [8]
end

% =========================================================================
% SECTION 3: SIMULATION PROCESSING LOOP
% =========================================================================

for i = 1:length(EbNo)
    
    % -- Calculate Noise Variance from Eb/No
    codeRate = K/E;
    bitsPerSymbol = 2; % For QPSK
    EsNo = EbNo(i) + 10*log10(bitsPerSymbol);
    snr_dB = EsNo + 10*log10(codeRate);
    noiseVar = 1./(10.^(snr_dB/10));
    
    % -- Instantiate AWGN Channel with the calculated noise variance
    channel = comm.AWGNChannel('NoiseMethod','Variance','Variance',noiseVar);
    
    % -- Loop counters for each Eb/No point
    frameCount = 0;
    errorCount = 0;
    
    fprintf('Simulating Eb/No = %.1f dB... \n', EbNo(i));

    while frameCount < maxNumFrames && errorCount < minFrameErrors
        
        % -----------------------------------------------------------------
        % TRANSMITTER
        % -----------------------------------------------------------------
        % 1. Generate random message bits
        msg = randi([0 1], K-crcLen, 1);
        
        % 2. Attach CRC
        msg_crc = nrCRCEncode(msg, crcPoly);
        
        % 3. Polar Encode and Rate Match
        % The nrPolarEncode function handles both encoding and rate matching to length E.
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
        % Provide perfect noise variance knowledge to the demodulator
        rxLLRs = nrSymbolDemodulate(receivedSymbols, modulation, noiseVar);
        
        % 8. Polar Decode (CA-SCL decoder, includes rate recovery)
        % Pass the LLRs directly from the demodulator to the decoder.
        % The nrPolarDecode function handles rate recovery internally.
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
    
    % -- Calculate and store BLER for this Eb/No point
    bler(i) = errorCount / frameCount;
    % disp();
end

% =========================================================================
% SECTION 4: PLOT RESULTS
% =========================================================================
figure;
semilogy(EbNo, bler, 'o-');
grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate (BLER)');
legend(['K=' num2str(K) ', E=' num2str(E) ', L=' num2str(L)]);
ylim([1e-5 1]);
