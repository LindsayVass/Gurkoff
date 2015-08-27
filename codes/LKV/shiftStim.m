function [shiftedStimData, timeShifts, meanStim] = shiftStim(stimData)
% shiftStim - shift each stimulation in time so that it has maximal
% correlation with meanStim
% >> [shiftedStimData, timeShifts] = shiftStim(stimData)
%
% Inputs:
%   stimData: matrix in which each row contains a timeseries of the raw
%       data surrounding the previously identified stimulation peak
%
% Optional Input:
%   maxShift: (default = 3) max number of samples to shift the data by
%
% Output:
%   shiftedStimData: zero-padded matrix containing all the aligned
%       stimulation data
%   timeShifts: column vector indicating the number of samples each row of
%       the data in stimData was shifted
%   meanStim: average stimulation peak (before time shifting)
%
% Lindsay Vass
% 27 August 2015

if ~exist('maxShift', 'var')
    maxShift = 3;
end

meanStim  = mean(stimData, 1);

tempShift  = cell(size(stimData, 1), 1);
timeShifts = zeros(size(stimData, 1), 1);

maxL = 0;
maxR = 0;

for i = 1:size(stimData, 1)
    % identify best shift
    thisStim  = stimData(i, :);
    [r, lags] = xcorr(meanStim, thisStim, maxShift);
    thisShift = lags(r == max(r));
    timeShifts(i) = thisShift;
    
    % update max shift
    if thisShift < 0 & abs(thisShift) > maxL
        maxL = abs(thisShift);
    end
    
    if thisShift > 0 & thisShift > maxR
        maxR = thisShift;
    end
    
    % shift the data
    if thisShift == 0
        shiftedStim = thisStim;
    elseif thisShift > 0
        shiftedStim = [zeros(1, thisShift), thisStim];
    else
        shiftedStim = [thisStim, zeros(1, abs(thisShift))];
    end
    
    % insert into the cell array
    tempShift{i} = shiftedStim;
end

shiftedStimData = zeros(size(stimData, 1), size(meanStim, 2) + maxL + maxR);

for i = 1:size(shiftedStimData, 1)
    
    thisShift = timeShifts(i);
    thisData  = tempShift{i};
    
    if thisShift == 0
        shiftedStimData(i, :) = [zeros([1 maxL]) thisData zeros([1 maxR])];
    elseif thisShift > 0
        shiftedStimData(i, :) = [zeros([1 maxL]) thisData zeros([1 maxR - thisShift])];
    else
        shiftedStimData(i, :) = [zeros([1 maxL + thisShift]) thisData zeros([1 maxR])];
    end
    
end

meanStim = [zeros([1 maxL]) meanStim zeros([1 maxR])];