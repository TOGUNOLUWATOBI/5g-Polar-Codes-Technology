%% Evaluation Task 5: Impact of UE Mobility on PDCCH List Decoding with TDL-C
%
% Goal: Analyze how UE speed (Doppler) degrades the performance
%       gains from polar list decoding in a realistic fading channel.
%
% Channel: nrTDLChannel (3GPP Tapped Delay Line - TDL-C)
% Fixed:   Eb/No = 15 dB (higher SNR needed for fading channel)
% X-axis:  UE Speed (km/h)
% Y-axis:  BLER
% Lines:   Decoder List Size (L)
%

clear;
clc;
close all;

%% Simulation Parameters
% ---------------------
% Fixed Parameters
EbNo_dB = 12;          % Fixed Eb/No in dB (TDL-C needs ~6-8 dB more than AWGN due to fading)
carrierFreq = 4e9;      % Carrier frequency (e.g., 4 GHz, FR1)
K = 40;                 % DCI payload size in bits (e.g., Format 1_1)
RNTI = 1234;            % A test RNTI

% Channel Selection: 'AWGN' or 'TDL'
channelType = 'TDL';    % Use TDL-C fading channel

% Monte Carlo Parameters
simParams.minErrors = 100;  % Minimum errors to collect per point
simParams.maxFrames = 2000; % Maximum frames per point

% --- Loop Variables (as requested) ---
ueSpeeds_kmh = [3, 30, 60, 120]; % X-axis
listSizes = [2, 4, 8, 16, 32];  % Lines on graph

% Results storage
results.bler = zeros(length(ueSpeeds_kmh), length(listSizes));

%% 5G NR Carrier and PDCCH Configuration
% -------------------------------------
% Carrier configuration
carrier = nrCarrierConfig;
carrier.NSizeGrid = 52;     % Bandwidth in RBs
carrier.SubcarrierSpacing = 30; % 30 kHz SCS
carrier.CyclicPrefix = 'Normal';

% PDCCH configuration
pdcch = nrPDCCHConfig;
pdcch.NSizeBWP = 52; % BWP size
pdcch.RNTI = RNTI;

% Calculate total RBGs in the BWP
numRBGs_BWP = ceil(pdcch.NSizeBWP / 6); % = 9 in this case

% Create a frequency bitmap for the CORESET
numCORESETRBGs = 2; % Use 2 RBGs (12 PRBs)
freqResBitmap = zeros(1, numRBGs_BWP);
freqResBitmap(1:numCORESETRBGs) = 1;

% Apply CORESET settings
pdcch.CORESET.CORESETID = 0;
pdcch.CORESET.FrequencyResources = freqResBitmap;
pdcch.CORESET.Duration = 2; % 2 symbols
pdcch.CORESET.RBOffset = 0;
pdcch.CORESET.CCEREGMapping = 'interleaved';
pdcch.CORESET.REGBundleSize = 6;

% Configure a single candidate (AL=4) for this test
pdcch.SearchSpace.SearchSpaceID = 1;
pdcch.SearchSpace.CORESETID = 0;
pdcch.SearchSpace.NumCandidates = [0 0 1 0 0]; % One candidate at AL=4
pdcch.SearchSpace.StartSymbolWithinSlot = 0;
pdcch.SearchSpace.Duration = 1;
pdcch.AggregationLevel = 4;
pdcch.AllocatedCandidate = 1;

% Get resource indices and calculate E
[pdcchInd, ~] = nrPDCCHResources(carrier, pdcch);

% Check if any indices were returned
if isempty(pdcchInd)
     error(['No valid PDCCH candidate found. ' ...
           'Check SearchSpace config (AL=4, Cand=1).']);
end

% E is the number of REs times bits/RE (2 for QPSK)
E = numel(pdcchInd) * 2;

% Check if payload fits
if E < K
    error('Payload size K=%d is too large for coded bits E=%d. Increase AL or CORESET resources.', K, E);
end

% Calculate fixed code rate and SNR
R = K / E;
EbNo_lin = 10^(EbNo_dB / 10);
snr_lin = EbNo_lin * R * 2; % 2 for QPSK (log2(M))
noiseVar = 1.0 / snr_lin;

fprintf('Simulation Parameters:\n');
fprintf('  K = %d bits, E = %d bits, Code Rate R = %.4f\n', K, E, R);
fprintf('  Eb/No = %d dB, SNR = %.2f dB, Noise Variance = %.4f\n', EbNo_dB, 10*log10(snr_lin), noiseVar);
fprintf('  Channel: %s\n\n', channelType);

%% 3GPP TDL Channel Configuration
% ------------------------------
channel = nrTDLChannel;
channel.DelayProfile = 'TDL-C'; % Common non-line-of-sight profile
channel.DelaySpread = 300e-9;   % 300 ns delay spread
channel.NumReceiveAntennas = 1;
channel.NumTransmitAntennas = 1;

% Get channel sample rate from OFDM info
ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;

%% Main Simulation Loop
% --------------------
hWait = waitbar(0, 'Starting TDL-C simulation...');
totalSims = length(ueSpeeds_kmh) * length(listSizes);
simCount = 0;

for iSpeed = 1:length(ueSpeeds_kmh)
    speed_kmh = ueSpeeds_kmh(iSpeed);
    
    % Convert speed (km/h) to max Doppler shift (Hz)
    speed_ms = speed_kmh * 1000 / 3600;
    dopplerShift = (speed_ms / 3e8) * carrierFreq;
    
    % Release channel object before changing non-tunable property
    release(channel);
    channel.MaximumDopplerShift = dopplerShift;
    
    fprintf('Simulating Speed: %d km/h (Doppler: %.1f Hz)\n', speed_kmh, dopplerShift);

    for iList = 1:length(listSizes)
        L = listSizes(iList);
        
        numErrors = 0;
        numFrames = 0;
        
        % Reset channel for each new simulation point
        reset(channel);

        while numErrors < simParams.minErrors && numFrames < simParams.maxFrames
            
            % 1. Generate DCI payload
            dciBits = randi([0 1], K, 1);
            
            % 2. DCI Encode (Polar)
            codedBits = nrDCIEncode(dciBits, RNTI, E);
            
            % 3. Modulate and scramble
            pdcchSym = nrPDCCH(codedBits, pdcch.RNTI, pdcch.RNTI);
            
            % 4. Map to grid
            txGrid = nrResourceGrid(carrier);
            txGrid(pdcchInd) = pdcchSym;

            % 5. Channel processing
            if strcmp(channelType, 'AWGN')
                % AWGN Channel - direct symbol-level processing
                noise = sqrt(noiseVar/2) * (randn(size(pdcchSym)) + 1i*randn(size(pdcchSym)));
                rxSym = pdcchSym + noise;
                eqSym = rxSym;
                
            else
                % TDL Channel - full OFDM processing
                txWaveform = nrOFDMModulate(carrier, txGrid);
                
                % Pass through TDL channel
                [rxWaveform, pathGains, sampleTimes] = channel(txWaveform);
                
                % Perfect Channel Estimation (with timing offset)
                offset = channel.info.ChannelFilterDelay;
                pathFilters = getPathFilters(channel);
                hest = nrPerfectChannelEstimate(carrier, pathGains, pathFilters, offset);
                
                % Timing Synchronization and OFDM Demodulation
                % Trim the received waveform by the channel filter delay
                rxWaveform_trim = rxWaveform(1+offset:end, :);
                rxGrid = nrOFDMDemodulate(carrier, rxWaveform_trim);
                
                % Ensure dimensions match between rxGrid and hest
                if size(hest, 2) > size(rxGrid, 2)
                    % Trim hest to match rxGrid
                    hest = hest(:, 1:size(rxGrid, 2));
                end
                
                % Add AWGN noise AFTER demodulation (to match noiseVar scale)
                noise = sqrt(noiseVar/2) * (randn(size(rxGrid)) + 1i*randn(size(rxGrid)));
                rxGrid = rxGrid + noise;
                
                % Extract PDCCH symbols and channel estimates
                rxSym = nrExtractResources(pdcchInd, rxGrid);
                hestSym = nrExtractResources(pdcchInd, hest);
                
                % Equalize
                eqSym = rxSym ./ hestSym;
            end
            
            % 6. Demodulate QPSK to get LLRs
            scrambledLLRs = nrSymbolDemodulate(eqSym, 'QPSK', noiseVar);
            
            % 7. Descramble the LLRs (at bit level)
            [scramSeqBin, ~] = nrPDCCHPRBS(pdcch.RNTI, pdcch.RNTI, E);
            rxLLRs = scrambledLLRs;
            rxLLRs(scramSeqBin == 1) = -rxLLRs(scramSeqBin == 1);
            
            % 8. Polar Decode (using specified List Size)
            [decodedBits, crcErr] = nrDCIDecode(rxLLRs, K, L, RNTI);
            
            % 9. Count errors
            if crcErr ~= 0
                numErrors = numErrors + 1;
            end
            numFrames = numFrames + 1;
        end
        
        % Store result
        results.bler(iSpeed, iList) = numErrors / numFrames;
        
        % Update waitbar
        simCount = simCount + 1;
        waitbar(simCount / totalSims, hWait, ...
            sprintf('Speed: %d km/h, L=%d. BLER=%.2e', ...
            speed_kmh, L, results.bler(iSpeed, iList)));
    end
end
close(hWait);
disp('Simulation finished.');

%% Plotting Results
% ----------------
% Create legend entries
legendStr = "L = " + string(listSizes);

figure;
% Plot BLER vs. Speed, with one line per List Size
hPlot = semilogy(ueSpeeds_kmh, results.bler, 'o-');
set(hPlot, 'LineWidth', 1.5);

grid on;
xlabel('UE Speed (km/h)');
ylabel('PDCCH BLER');
title(sprintf('PDCCH BLER vs. UE Speed (Eb/No = %d dB, %s)', EbNo_dB, channelType));
legend(legendStr, 'Location', 'northwest');
ylim([simParams.minErrors/simParams.maxFrames*0.5 1]); % Set logical Y-axis limits

% Save results to file
save('Task5_UEMobility_TDL_Results.mat', 'results', 'ueSpeeds_kmh', 'listSizes', 'EbNo_dB', 'channelType');
fprintf('Results saved to Task5_UEMobility_TDL_Results.mat\n');

% Save figure
saveas(gcf, 'Task5_UEMobility_TDL_BLER_vs_Speed.png');
fprintf('Figure saved to Task5_UEMobility_TDL_BLER_vs_Speed.png\n');

% Display results table
fprintf('\n=== Simulation Results ===\n');
fprintf('Eb/No = %d dB, Channel: %s\n\n', EbNo_dB, channelType);
fprintf('Speed (km/h) | ');
fprintf('L=%d\t\t', listSizes);
fprintf('\n');
fprintf(repmat('-', 1, 60));
fprintf('\n');
for iSpeed = 1:length(ueSpeeds_kmh)
    fprintf('%12d | ', ueSpeeds_kmh(iSpeed));
    for iList = 1:length(listSizes)
        fprintf('%.2e\t', results.bler(iSpeed, iList));
    end
    fprintf('\n');
end

fprintf('\n=== Analysis ===\n');
fprintf('Expected behavior with TDL-C channel:\n');
fprintf('- BLER should increase with UE speed (higher Doppler)\n');
fprintf('- Higher list sizes (L) should provide better performance\n');
fprintf('- At high speeds, the benefit of large L may diminish\n');
