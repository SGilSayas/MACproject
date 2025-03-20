function outputsignal = ampreduce(varargin)
%ampreduce   Reduces the amplitude of a signal
%
%   For vectors, AMPREDUCE(X) is a vector with the same mean but 
%   smaller range than the original vector. The range is reduced by a factor
%   defined by the user. If the user does not define a factor, the range is reduced by a factor of 2
%
%   Example: If X = [0 2 4]
%
%   then ampreduce(X) is [1 2 3]
%
%   See attached license for terms of use.


%Initial checks
if nargin == 1
    signal = varargin{1};
    reduction = 2;
elseif nargin == 2
    signal = varargin{1};
    reduction = varargin{2};
else
    error('Wrong number of inputs');
end

%Processing...
temp_signal = signal / reduction;
outputsignal = (mean(signal)) - (mean(temp_signal))  + temp_signal;

end