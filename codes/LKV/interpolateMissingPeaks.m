function [finalPeakLocs, badEpochs] = interpolateMissingPeaks(waveform, maxAutocorr, newPeakLocs, manualOrAuto, tol)
% interpolateMissingPeaks - add peaks where they failed to be identified
% >> finalPeakLocs = interpolateMissingPeaks(waveform, maxAutocorr, newPeakLocs, badPeaks)
%
% Inputs:
%   waveform: raw signal containing the data from which the artifact is to
%       be removed
%   maxAutocorr: expected lag in EEG samples between peaks
%   newPeakLocs: output of peakFinder containing the indices of waveform
%       where stimulation peaks are located
%   manualOrAuto: 'manual' = manually reject questionable stimulation peaks
%       or 'auto' = automatically reject ALL questionable stimulation peaks
%       (more conservative, but less time consuming)
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
    
    if strcmpi(manualOrAuto, 'auto') == 1
        badEpochs = cat(1, badEpochs, [firstPeak lastPeak]);
    else
        
        if numMissingInt < tol
            newPeaks = [firstPeak:maxAutocorr:lastPeak]';
            
            if lastPeak - newPeaks(end) <= tol * maxAutocorr
                newPeaks(end) = [];
            end
            
            h = plotAllPeaks(waveform, newPeaks, maxAutocorr);
            prompt = 'Keep or reject peaks?';
            answer = questdlg(prompt, 'Keep peaks?', 'Keep All', 'Reject All', 'Keep All');
            
            if strcmpi(answer, 'Keep All') == 1
                finalPeakLocs = [finalPeakLocs(1:firstPeakInd); newPeaks; finalPeakLocs(lastPeakInd:end)];
                finalPeakLocs = unique(finalPeakLocs);
            else
                badEpochs = cat(1, badEpochs, [firstPeak lastPeak]);
            end
            
            close(h);
        else
            badEpochs = cat(1, badEpochs, [firstPeak lastPeak]);
        end
    end
end

function h = plotAllPeaks(waveform, thePeaks, maxAutocorr, plotXBuffer, plotYBuffer)

if ~exist('plotXBuffer', 'var')
    plotXBuffer = 10;
end

if ~exist('plotYBuffer', 'var')
    plotYBuffer = 2;
end

h = figure;
plot(waveform);
hold on;

for i = 1:length(thePeaks)
    plot(thePeaks, waveform(thePeaks), 'k^', 'markerfacecolor', [1 0 0]);
end
hold off;
axis([thePeaks(1) - plotXBuffer * maxAutocorr thePeaks(end) + plotXBuffer * maxAutocorr mean(waveform) - plotYBuffer*std(waveform) inf])
