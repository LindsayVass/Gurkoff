function [corr, lags] = autocorrel(y,firstlag,lastlag)
% autocorrel - Compute Autocorrelations Through Lags
% >> [corr, lags] = autocorrel(y, firstlag, lastlag)
%
% Inputs:
% y - series to compute autocorrelation for, nx1 column vector
% firstlag - starting value for lag, 1x1 integer
% lastlag - last value for lag, 1x1 integer
%
% Output:
% corr - vector containing autocorrelations
% lags - vector containing lags at which autocorrelation was calculated
%
%
% Example:
% >> acf(randn(100,1), 5, 10)
%
% Original code from: http://www.mathworks.com/matlabcentral/fileexchange/30540-autocorrelation-function--acf-
%
% Modified by Lindsay Vass (25 August 2015)


% --------------------------
% USER INPUT CHECKS
% --------------------------

[n1, n2] = size(y) ;
if n2 ~=1
    error('Input series y must be an nx1 column vector')
end

[a1, a2] = size(lastlag) ;
if ~((a1==1 & a2==1) & (lastlag<n1))
    error('Input number of lastlag must be a 1x1 scalar, and must be less than length of series y')
end

[b1, b2] = size(firstlag);
if ~((b1 == 1 & b2 == 1) & firstlag < lastlag)
    error('Input number of firstlag must be a 1x1 scalar, and must be less than the lastlag')
end



% -------------
% BEGIN CODE
% -------------
lags = [firstlag:1:lastlag];
corr = zeros(length(lags),1) ;
N = max(size(y)) ;
ybar = mean(y); 

% Collect ACFs at each lag i
for i = firstlag:lastlag 
cross_sum=(y(i+1:N)-ybar)'*(y(1:N-i)-ybar); 
yvar = (y-ybar)'*(y-ybar) ; 
corr(i - firstlag + 1) = (cross_sum / yvar)*(N/(N-i)) ; 
end