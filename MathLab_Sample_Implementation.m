% A basic MATLAB script to demonstrate a simple Polar Code simulation loop.
% --- REVISED VERSION ---
% This version uses the core polar coding functions directly to avoid
% potential versioning issues and to focus on the key components.

% --- 1. Simulation Parameters ---
K = 50;         % Message length in bits (the actual information)
E = 128;        % Rate-matched output length (the length after encoding)
EbNo_dB = 2;    % Eb/No in decibels - THIS IS THE SIGNAL QUALITY!
L = 8;          % Decoder list length - A KEY PARAMETER FOR YOU!

% Modulation scheme
modulation = 'QPSK';
bitsPerSymbol = 2; % For QPSK

% --- Simulation Setup ---
num_block_errors = 0;
num_sim_blocks = 1000; % Let's run a few more blocks for a better result

fprintf('Starting simulation for EbNo = %.1f dB\n', EbNo_dB);

% --- 2. Main Simulation Loop ---
for block_idx = 1:num_sim_blocks
    
    % --- Transmitter Side ---
    
    % Generate a random message
    tx_bits = randi([0 1], K, 1);
    
    % --- Core Polar Encoding Step ---
    % This is the fundamental function for Polar encoding.
    encoded_bits = nrPolarEncode(tx_bits, E);

    % Modulate
    symbols = nrSymbolModulate(encoded_bits, modulation);

    % --- Channel ---
    
    % Calculate noise variance from Eb/No
    codeRate = K/E;
    snr_dB = EbNo_dB + 10*log10(codeRate) + 10*log10(bitsPerSymbol);
    noiseVar = 10.^(-snr_dB/10);

    % Add Additive White Gaussian Noise (AWGN)
    rx_symbols = awgn(symbols, snr_dB, 'measured');

    % --- Receiver Side ---
    
    % Demodulate - get Log-Likelihood Ratios (LLRs)
    llrs = nrSymbolDemodulate(rx_symbols, modulation, noiseVar);

    % --- Core Polar Decoding Step ---
    % This is the fundamental function for Polar decoding.
    % We tell it the expected message length (K), the input length (E),
    % and the list size (L) to use.
    rx_bits = nrPolarDecode(llrs, K, E, L);

    % --- 3. Tally Errors ---
    % Check if the decoded bits are different from the transmitted ones.
    if any(rx_bits ~= tx_bits)
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

