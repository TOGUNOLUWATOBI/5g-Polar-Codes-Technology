% A basic MATLAB script to demonstrate a simple Polar Code simulation loop.

% --- 1. Simulation Parameters ---
K = 50;         % Message length in bits
E = 128;        % Rate-matched output length
EbNo_dB = 2;    % Eb/No in decibels - THIS IS THE SIGNAL QUALITY!
L = 8;          % Decoder list length - 

% Modulation scheme
modulation = 'QPSK';
info = nrDLSCHInfo(K, 1/3); % Using a nominal rate for DLSCH info

% --- Simulation Setup ---
num_block_errors = 0;
num_sim_blocks = 100; % Number of blocks to simulate for a quick test
                      % For real results, this should be much higher (e.g., 10000)

fprintf('Starting simulation for EbNo = %.1f dB\n', EbNo_dB);

% --- 2. Main Simulation Loop ---
for block_idx = 1:num_sim_blocks
    
    % --- Transmitter Side ---
    
    % Generate a random message (transport block)
    trBlk = randi([0 1], K, 1);
    
    % DLSCH Encoding (includes CRC attachment, segmentation, and Polar encoding)
    codedTrBlk = nrDLSCH(trBlk, info.Rate, E);

    % Modulate
    symbols = nrSymbolModulate(codedTrBlk, modulation);

    % --- Channel ---
    
    % Calculate noise variance from Eb/No
    % This conversion depends on the code rate and bits per symbol
    bitsPerSymbol = 2; % For QPSK
    codeRate = K/E;
    snr_dB = EbNo_dB + 10*log10(codeRate) + 10*log10(bitsPerSymbol);
    noiseVar = 10.^(-snr_dB/10);

    % Add Additive White Gaussian Noise (AWGN)
    % This is the simplest channel model.
    rx_symbols = awgn(symbols, snr_dB, 'measured');

    % --- Receiver Side ---
    
    % Demodulate - get Log-Likelihood Ratios (LLRs)
    llrs = nrSymbolDemodulate(rx_symbols, modulation, noiseVar);

    % DLSCH Decoding (the core part for you)
    % This function performs rate recovery and then Polar decoding.
    % The 'L' parameter is passed to the underlying Polar decoder.
    [decTrBlk, blk_err] = nrDLSCHDecode(llrs, info, E, L);

    % --- 3. Tally Errors ---
    if blk_err
        num_block_errors = num_block_errors + 1;
    end
    
end

% --- 4. Display Results ---
BLER = num_block_errors / num_sim_blocks;

fprintf('Simulation Finished!\n');
fprintf('---------------------\n');
fprintf('Signal Quality (Eb/No): %.1f dB\n', EbNo_dB);
fprintf('Decoder List Size (L): %d\n', L);
fprintf('Total Blocks Simulated: %d\n', num_sim_blocks);
fprintf('Blocks in Error: %d\n', num_block_errors);
fprintf('Block Error Rate (BLER): %.4f\n', BLER);
