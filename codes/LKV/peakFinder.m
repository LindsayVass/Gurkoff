function peakLocs = peakFinder(waveform, maxAutocorr, peakLocs, startingPeak)
% peakFinder - ask user to identify first and last stimulation pulse using
% interactive interface
% >> peakLoc = peakFinder(waveform, maxAutocorr, peakLocs, startingPeak)
%
% Inputs:
%   waveform: raw signal containing the data from which the artifact is to
%       be removed
%   maxAutocorr: expected lag in EEG samples between peaks
%   peakLocs: indices of the peaks identified by findpeaks
%   startingPeak: either 'first' if identifying the first peak of the
%       dataset or 'last' if identifying the last peak of the dataset
%
% Output:
%   peakLoc: index of the starting peak identified by the user; if the
%       algorithm correctly identified the starting peak, this value will
%       be the same as "startingPeak"
%
% Lindsay Vass
% 26 August 2015

% Plot the starting peak
if strcmpi(startingPeak, 'first') == 1
    startPeakInd = 1;
elseif strcmpi(startingPeak, 'last') == 1
    startPeakInd = length(peakLocs);
else
    error('startingPeak must be set to either ''first'' or ''last''');
end

startPeak = peakLocs(startPeakInd);

% Ask user to identify if starting peak is correctly labeled
h = plotPeak(waveform, startPeak, maxAutocorr);

while 1
    
    % Ask user if peak is correctly identified
    if strcmpi(startingPeak, 'first') == 1
        prompt = 'Is the first stimulation peak correctly labeled?';
    else
        prompt = 'Is the last stimulation peak correctly labeled?';
    end
    
    answer = questdlg(prompt);
    
    % Ask user which direction to shift peak
    if strcmpi(answer, 'Cancel') == 1
        error('User ended peak finding operation.');
    elseif strcmpi(answer, 'Yes') == 1
        if strcmpi(startingPeak, 'first') == 1
            if startPeakInd > 1
                peakLocs(1:startPeakInd - 1) = [];
            end
        else
            if startPeakInd < length(peakLocs)
                peakLocs(startPeakInd + 1:end) = [];
            end
        end
            
        close(h);
        return
    else
        prompt  = 'Which direction should we shift the peak?';
        peakDir = questdlg(prompt, 'Move Peak', 'Forward', 'Backward', 'Forward');
    end
    
    % Shift peak
    close(h);
    warnText = 'There are no pre-identified peaks remaining. Estimating peak location based on expected lag between peaks.';
    if strcmpi(peakDir, 'Forward') == 1
        
        startPeakInd = startPeakInd + 1;
        
        if startPeakInd > length(peakLocs)
            warning(warnText);
            startPeak = startPeak + maxAutocorr;
            peakLocs(end + 1) = startPeak;
            startPeakInd = length(peakLocs);
        elseif startPeakInd < 1
            warning(warnText);
            startPeak = startPeak - maxAutocorr;
            peakLocs = cat(1, startPeak, peakLocs);
            startPeakInd = 1;
        else
            startPeak = peakLocs(startPeakInd);
        end
        
        h = plotPeak(waveform, startPeak, maxAutocorr);
    else
        startPeakInd = startPeakInd - 1;
        
        if startPeakInd < 1
            warning(warnText);
            startPeak = startPeak - maxAutocorr;
            peakLocs = cat(1, startPeak, peakLocs);
            startPeakInd = 1;
        elseif startPeakInd > length(peakLocs)
             warning(warnText);
            startPeak = startPeak + maxAutocorr;
            peakLocs(end + 1) = startPeak;
            startPeakInd = length(peakLocs);
        else
            startPeak = peakLocs(startPeakInd);
        end
        
        h = plotPeak(waveform, startPeak, maxAutocorr);
    end
    
end

function h = plotPeak(waveform, startPeak, maxAutocorr, plotXBuffer, plotYBuffer)

if ~exist('plotXBuffer', 'var')
    plotXBuffer = 30;
end

if ~exist('plotYBuffer', 'var')
    plotYBuffer = 2;
end

h = figure;
plot(waveform);
hold on;
plot(startPeak, waveform(startPeak), 'k^', 'markerfacecolor', [1 0 0]);
hold off;
axis([startPeak - plotXBuffer * maxAutocorr startPeak + plotXBuffer * maxAutocorr mean(waveform) - plotYBuffer*std(waveform) inf])

