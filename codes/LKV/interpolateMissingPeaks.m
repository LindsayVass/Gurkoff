function [finalPeakLocs, badEpochs] = interpolateMissingPeaks(maxAutocorr, newPeakLocs, tol)
% interpolateMissingPeaks - add peaks where they failed to be identified
% >> finalPeakLocs = interpolateMissingPeaks(waveform, maxAutocorr, newPeakLocs, badPeaks)
%
% Inputs:
%   waveform: raw signal containing the data from which the artifact is to
%       be removed
%   maxAutocorr: expected lag in EEG samples between peaks
%   newPeakLocs: output of peakFinder containing the indices of waveform
%       where stimulation peaks are located
%
% Optional Input:
%   tol: (default = 0.01) tolerance for identifying peaks at the expected
%       intervals; higher tolerance means the algorithm will accept peaks
%       that are further apart than expected. For example, for a 7.7 Hz
%       oscillation at 1000 Hz sampling rate, 0.01 allows peaks separated
%       by 131 samples rather than 130.
%
% Output:
%   finalPeakLocs: updated version of newPeakLocs containing indices of all
%        stimulation peaks in waveform
%   badEpochs: two column vector containing intervals where stimulation was
%       not found at the expected frequency; these intervals typically
%       correspond to artifacts in the data (e.g., saturation); first
%       column is start index in waveform and second column is end index in
%       waveform
%
% Lindsay Vass
% 26 August 2015

if ~exist('tol', 'var')
    tol = 0.01;
end

peakDiff = diff(newPeakLocs);
badPeaks = find(peakDiff ~= maxAutocorr);
numMissingPeaks = peakDiff(badPeaks) / maxAutocorr;

finalPeakLocs = newPeakLocs;

badEpochs = [];
for thisPeak = 1:length(badPeaks)
    
    firstPeak = newPeakLocs(badPeaks(thisPeak));
    lastPeak  = newPeakLocs(badPeaks(thisPeak) + 1);
    
    firstPeakInd = find(firstPeak == finalPeakLocs);
    lastPeakInd  = find(lastPeak == finalPeakLocs);
    
    numMissingInt = abs(round(numMissingPeaks(thisPeak)) - numMissingPeaks(thisPeak));
    
    if numMissingInt < tol
        newPeaks = [firstPeak:maxAutocorr:lastPeak]';
        if lastPeak - newPeaks(end) <= tol * maxAutocorr
            newPeaks(end) = [];
        end
        finalPeakLocs = [finalPeakLocs(1:firstPeakInd); newPeaks; finalPeakLocs(lastPeakInd:end)];
        finalPeakLocs = unique(finalPeakLocs);
    else
        badEpochs = cat(1, badEpochs, [firstPeak lastPeak]);
    end
end