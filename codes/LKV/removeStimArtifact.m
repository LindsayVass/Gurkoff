function [newWaveform, badEpochs, stimEpochs] = removeStimArtifact(fileName, waveform, stimFreq, samplingRate, manualOrAuto, preSpike, postSpike, autocorrBuffer, stimBuffer)
% removeStimArtifact - remove stimulation artifact using the method in
%   Wichmann (2000) J Neurosci Methods
% >> newWaveform = removeStimArtifact(waveform, stim_freq, samplingRate)
%
% Inputs:
%   fileName: short name of the recording being analyzed; this is only used
%       to label figures
%   waveform: raw signal containing the data from which the artifact is to
%       be removed
%   stim_freq: frequency in Hz of the stimulation
%   samplingRate: frequency in Hz of the recording sampling rate
%   manualOrAuto: 'manual' = manually reject questionable stimulation peaks
%       or 'auto' = automatically reject ALL questionable stimulation peaks
%       (more conservative, but less time consuming)
%   preSpike: number of samples before the identified stimulation peak to
%       include in interpolation (suggested = 3)
%   postSpike: number of samples after the identified stimulation peak to
%       include in interpolation (suggested = 12)
%
% Optional Input:
%   autocorrBuffer: (default = 50) number of lags to test above and below the
%       expected lag between stimulation peaks; for example, for stimFreq =
%       7.7 Hz and samplingRate = 1000 Hz, the expected lag is 130 samples
%       between stimulation peaks. A buffer of 50 will test lags of 80 -
%       180 samples.
%   stimBuffer: (default = 100) number of samples before and after
%       stimulation peak to include in calculation to remove stimulation
%       artifact
%
% Output:
%   newWaveform: waveform with stimulation artifact removed
%   badEpochs: two column vector containing intervals where stimulation was
%       not found at the expected frequency; these intervals typically
%       correspond to artifacts in the data (e.g., saturation); first
%       column is start index in waveform and second column is end index in
%       waveform; you will want to exclude these time points from future
%       analyses
%   stimEpochs: two column vector containing intervals where stimulation
%       was identified and removed; this could be used to exclude these
%       timepoints from later analyses
%
% Lindsay Vass
% 25 August 2015

if ~exist('autocorrBuffer', 'var')
    autocorrBuffer = 50;
end

if ~exist('stimBuffer', 'var')
    stimBuffer = 100;
end

% expected lag in samples
expAutocorr = round(samplingRate / stimFreq); 


%% identify max autocorrelation in signal, which should match expAutocorr

% this is basically a sanity check to make sure that we find stimulation at
% the expected frequency
[autocorr, lags] = autocorrel(waveform, expAutocorr - autocorrBuffer, expAutocorr + autocorrBuffer);
maxAutocorr = lags(autocorr == max(autocorr));

if maxAutocorr ~= expAutocorr
    warning(['Based on stimulation frequency, expected a lag of ' num2str(expAutocorr) ' samples between stimulations, but identified a lag of ' num2str(maxAutocorr) '. Are you sure there''s stimulation in this dataset?'])
end

%% identify time points of stimulations

% peaks should be separated by maxAutocorr
[~, locs] = findpeaks(waveform, 'minpeakdistance', maxAutocorr - 1);
locdiff   = diff(locs);
peakInds = find(locdiff == maxAutocorr);

% ask user to confirm first and last peak
peakLocs = locs(peakInds);
newPeakLocs = peakFinder(waveform, maxAutocorr, peakLocs, 'first');
newPeakLocs = peakFinder(waveform, maxAutocorr, newPeakLocs, 'last');

% handle peaks separated by more than expAutocorr
[finalPeakLocs, badEpochs] = interpolateMissingPeaks(waveform, maxAutocorr, newPeakLocs, manualOrAuto);

% report how much data is excluded
badTime = sum(badEpochs(:, 2) - badEpochs(:, 1)) / samplingRate;
totTime = length(waveform) / samplingRate;
badPct  = (badTime / totTime) * 100;

warning(['Excluding ' num2str(badTime) ' seconds of data, ' num2str(badPct) '% of total data.'])

%% remove artifact

% select intervals around stimulations
stimInds  = repmat(finalPeakLocs, [1, stimBuffer * 2 + 1]);
sampShift = repmat([-stimBuffer:stimBuffer], [size(stimInds, 1), 1]);
stimInds  = stimInds + sampShift;
stimInds(stimInds > length(waveform)) = length(waveform);
stimData  = waveform(stimInds);

startInterp = stimBuffer - preSpike;
stopInterp  = stimBuffer + postSpike;

deSpikedData = zeros(size(stimData));
for thisSpike = 1:size(deSpikedData, 1)
    thisData = stimData(thisSpike, :);
    x  = [1:startInterp stopInterp:length(thisData)];
    y  = thisData(x);
    xq = [startInterp + 1:stopInterp - 1];
    yq = interp1(x, y, xq, 'pchip');
    deSpikedData(thisSpike, :) = [thisData(1:startInterp) yq thisData(stopInterp:end)];    
end

% plot a few example spikes
spikeList = randsample(size(stimData, 1), 9);
figure;
for i = 1:length(spikeList)
    subplot(3, 3, i);
    plot(stimData(i, :));
    hold on;
    plot(deSpikedData(i, :), 'r');
    set(gca, 'xlim', [1 stimBuffer*2]);
end
suptitle({['Example artifact removals for ' fileName ], 'Blue = Before    Red = After'});

newWaveform = waveform;
newWaveform(stimInds) = deSpikedData;

stimEpochs = [finalPeakLocs - preSpike, finalPeakLocs + postSpike];