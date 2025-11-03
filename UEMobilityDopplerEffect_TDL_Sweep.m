%% Evaluation Task 5: Comprehensive Eb/No Sweep for UE Mobility Impact
%
% Goal: Analyze mobility impact across multiple Eb/No levels
%       This creates multiple figures showing BLER vs Speed at different SNR levels
%
% Channel: nrTDLChannel (3GPP Tapped Delay Line - TDL-C)
% X-axis:  UE Speed (km/h)
% Y-axis:  BLER
% Lines:   Decoder List Size (L)
% Sweep:   Multiple Eb/No values
%

clear;
clc;
close all;

%% Simulation Parameters
% ---------------------
% Sweep Parameters
EbNo_dB_sweep = [8, 10, 12, 14];  % Multiple Eb/No values to sweep
carrierFreq = 4e9;      % Carrier frequency (e.g., 4 GHz, FR1)
K = 40;                 % DCI payload size in bits (e.g., Format 1_1)
RNTI = 1234;            % A test RNTI

% Channel Selection
channelType = 'TDL';    % Use TDL-C fading channel

% Monte Carlo Parameters
simParams.minErrors = 100;  % Minimum errors to collect per point
simParams.maxFrames = 2000; % Maximum frames per point

% Loop Variables
ueSpeeds_kmh = [3, 30, 60, 120]; % X-axis
listSizes = [2, 4, 8, 16, 32];  % Lines on graph

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

% Calculate fixed code rate
R = K / E;
fprintf('Configuration: K = %d bits, E = %d bits, Code Rate R = %.4f\n\n', K, E, R);

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

%% Main Simulation Loop - Sweep Across Eb/No Values
% --------------------------------------------------
fprintf('=== Starting Comprehensive Eb/No Sweep ===\n\n');
fprintf('Sweeping Eb/No: %s dB\n', mat2str(EbNo_dB_sweep));
fprintf('UE Speeds: %s km/h\n', mat2str(ueSpeeds_kmh));
fprintf('List Sizes: %s\n\n', mat2str(listSizes));

% Storage for all results
allResults = struct();

for iEbNo = 1:length(EbNo_dB_sweep)
    EbNo_dB = EbNo_dB_sweep(iEbNo);
    
    % Calculate SNR and noise variance for this Eb/No
    EbNo_lin = 10^(EbNo_dB / 10);
    snr_lin = EbNo_lin * R * 2; % 2 for QPSK (log2(M))
    noiseVar = 1.0 / snr_lin;
    
    fprintf('========================================\n');
    fprintf('Eb/No = %d dB, SNR = %.2f dB\n', EbNo_dB, 10*log10(snr_lin));
    fprintf('========================================\n\n');
    
    % Results storage for this Eb/No
    results = struct();
    results.EbNo_dB = EbNo_dB;
    results.bler = zeros(length(ueSpeeds_kmh), length(listSizes));
    results.ueSpeeds_kmh = ueSpeeds_kmh;
    results.listSizes = listSizes;
    
    % Progress tracking
    totalSims = length(ueSpeeds_kmh) * length(listSizes);
    simCount = 0;
    hWait = waitbar(0, sprintf('Eb/No = %d dB: Starting...', EbNo_dB));
    
    for iSpeed = 1:length(ueSpeeds_kmh)
        speed_kmh = ueSpeeds_kmh(iSpeed);
        
        % Convert speed (km/h) to max Doppler shift (Hz)
        speed_ms = speed_kmh * 1000 / 3600;
        dopplerShift = (speed_ms / 3e8) * carrierFreq;
        
        % Release channel object before changing non-tunable property
        release(channel);
        channel.MaximumDopplerShift = dopplerShift;
        
        fprintf('  Speed: %d km/h (Doppler: %.1f Hz)\n', speed_kmh, dopplerShift);
        
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
                
                % 5. TDL Channel - full OFDM processing
                txWaveform = nrOFDMModulate(carrier, txGrid);
                
                % Pass through TDL channel
                [rxWaveform, pathGains, sampleTimes] = channel(txWaveform);
                
                % Perfect Channel Estimation (with timing offset)
                offset = channel.info.ChannelFilterDelay;
                pathFilters = getPathFilters(channel);
                hest = nrPerfectChannelEstimate(carrier, pathGains, pathFilters, offset);
                
                % Timing Synchronization and OFDM Demodulation
                rxWaveform_trim = rxWaveform(1+offset:end, :);
                rxGrid = nrOFDMDemodulate(carrier, rxWaveform_trim);
                
                % Ensure dimensions match between rxGrid and hest
                if size(hest, 2) > size(rxGrid, 2)
                    hest = hest(:, 1:size(rxGrid, 2));
                end
                
                % Add AWGN noise AFTER demodulation
                noise = sqrt(noiseVar/2) * (randn(size(rxGrid)) + 1i*randn(size(rxGrid)));
                rxGrid = rxGrid + noise;
                
                % Extract PDCCH symbols and channel estimates
                rxSym = nrExtractResources(pdcchInd, rxGrid);
                hestSym = nrExtractResources(pdcchInd, hest);
                
                % Equalize
                eqSym = rxSym ./ hestSym;
                
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
                sprintf('Eb/No=%d dB | Speed=%d km/h, L=%d | BLER=%.2e', ...
                EbNo_dB, speed_kmh, L, results.bler(iSpeed, iList)));
        end
    end
    close(hWait);
    
    % Store results for this Eb/No
    allResults.(sprintf('EbNo_%ddB', EbNo_dB)) = results;
    
    % Print results table for this Eb/No
    fprintf('\n--- Results for Eb/No = %d dB ---\n', EbNo_dB);
    fprintf('Speed (km/h) |');
    for L = listSizes
        fprintf(' L=%-6d', L);
    end
    fprintf('\n');
    fprintf('-------------|');
    for L = listSizes
        fprintf('--------');
    end
    fprintf('\n');
    for iSpeed = 1:length(ueSpeeds_kmh)
        fprintf('%12d |', ueSpeeds_kmh(iSpeed));
        for iList = 1:length(listSizes)
            fprintf(' %.2e\t', results.bler(iSpeed, iList));
        end
        fprintf('\n');
    end
    fprintf('\n');
    
    % Create and save individual figure for this Eb/No
    figure('Position', [100, 100, 800, 600]);
    legendStr = "L = " + string(listSizes);
    
    hPlot = semilogy(ueSpeeds_kmh, results.bler, 'o-');
    set(hPlot, 'LineWidth', 2, 'MarkerSize', 8);
    
    grid on;
    xlabel('UE Speed (km/h)', 'FontSize', 12);
    ylabel('PDCCH BLER', 'FontSize', 12);
    title(sprintf('PDCCH BLER vs. UE Speed\n(Eb/No = %d dB, TDL-C Channel)', EbNo_dB), 'FontSize', 14);
    legend(legendStr, 'Location', 'best', 'FontSize', 10);
    ylim([simParams.minErrors/simParams.maxFrames*0.5 1]);
    set(gca, 'FontSize', 11);
    
    % Save figure
    figFilename = sprintf('Task5_Sweep_EbNo_%ddB.png', EbNo_dB);
    saveas(gcf, figFilename);
    fprintf('Figure saved: %s\n\n', figFilename);
end

%% Create Combined Comparison Figure
% ----------------------------------
fprintf('Creating combined comparison figure...\n');

figure('Position', [100, 100, 1200, 800]);

% Create subplots for each Eb/No
numEbNo = length(EbNo_dB_sweep);
for iEbNo = 1:numEbNo
    EbNo_dB = EbNo_dB_sweep(iEbNo);
    results = allResults.(sprintf('EbNo_%ddB', EbNo_dB));
    
    subplot(2, 2, iEbNo);
    legendStr = "L = " + string(listSizes);
    
    hPlot = semilogy(ueSpeeds_kmh, results.bler, 'o-');
    set(hPlot, 'LineWidth', 1.5, 'MarkerSize', 6);
    
    grid on;
    xlabel('UE Speed (km/h)', 'FontSize', 11);
    ylabel('PDCCH BLER', 'FontSize', 11);
    title(sprintf('Eb/No = %d dB', EbNo_dB), 'FontSize', 12, 'FontWeight', 'bold');
    legend(legendStr, 'Location', 'best', 'FontSize', 9);
    ylim([simParams.minErrors/simParams.maxFrames*0.5 1]);
    set(gca, 'FontSize', 10);
end

sgtitle('PDCCH BLER vs. UE Speed Across Multiple Eb/No Levels (TDL-C Channel)', ...
    'FontSize', 16, 'FontWeight', 'bold');

% Save combined figure
saveas(gcf, 'Task5_Sweep_Combined.png');
fprintf('Combined figure saved: Task5_Sweep_Combined.png\n');

%% Save All Results
% ----------------
save('Task5_Sweep_AllResults.mat', 'allResults', 'EbNo_dB_sweep', ...
    'ueSpeeds_kmh', 'listSizes', 'simParams', 'channelType');
fprintf('All results saved: Task5_Sweep_AllResults.mat\n');

%% Final Summary
% -------------
fprintf('\n========================================\n');
fprintf('=== COMPREHENSIVE SWEEP COMPLETE ===\n');
fprintf('========================================\n\n');

fprintf('Generated Files:\n');
for iEbNo = 1:length(EbNo_dB_sweep)
    fprintf('  - Task5_Sweep_EbNo_%ddB.png\n', EbNo_dB_sweep(iEbNo));
end
fprintf('  - Task5_Sweep_Combined.png\n');
fprintf('  - Task5_Sweep_AllResults.mat\n\n');

fprintf('Simulation Parameters Summary:\n');
fprintf('  Eb/No Range: %d to %d dB\n', min(EbNo_dB_sweep), max(EbNo_dB_sweep));
fprintf('  Speed Range: %d to %d km/h\n', min(ueSpeeds_kmh), max(ueSpeeds_kmh));
fprintf('  List Sizes: %s\n', mat2str(listSizes));
fprintf('  Channel: TDL-C (300ns delay spread)\n');
fprintf('  Frames per point: up to %d\n', simParams.maxFrames);
fprintf('  Errors per point: at least %d\n\n', simParams.minErrors);

disp('All simulations completed successfully!');
