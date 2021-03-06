%[rasterdata, truefreq, points_per_line] = raster(freq, Ts, Nperiod, varargin)
%    
% 
%   Constructs a triangle wave of freqency freq, using sample period Ts for
%   duration [0, tend]. If 'coerce', 1 is passed, then the frequency is
%   coerced such that the number of points in a quarter period will be an integer 
%   multiple of the sampling frequency
% 
%  Inputs
%  ------
%   freq - desired triangle wave freqency
%   Ts     sample period
%   Nperiod total number of periods. Default is 1.
%   coerce 0, 1 (optional)
%   shift: a time to shift the raster pattern by. if shift= -1/4*period,
%           then the pattern will start at the negative peak value.
% 
%  Outputs
%  ------
%   rasterdata  a timeseries object containing the triangle wave time, and
%   y data
%   truefreq   if coerce = 1, then this contains the frequency that freq
%   was coerced into.


function [rasterdata, truefreq, points_per_line]...
                       = raster(freq, Ts, Nperiod, varargin)

if ~exist('Nperiod', 'var') || isempty(Nperiod)
    Nperiod = 1;
end
defcoerce = 0;
defshift = 0;
p = inputParser;
addParameter(p, 'coerce', defcoerce);
addParameter(p, 'shift', defshift);

parse(p, varargin{:});
opts = p.Results;

period = 1/freq;
if opts.coerce
    % We want to get points at the vertices, so coerce the triangle frequency
    % to so that the number of points in a quarter period will be an integer 
    % multiple of the sampling frequency
    quarterperiod = period/4;
    pointsperquarterperiod = round(quarterperiod/Ts);
    period = 4*(pointsperquarterperiod*Ts);
else
    quarterperiod = period/4;
    pointsperquarterperiod = (quarterperiod/Ts);
end


points_per_line = 2*pointsperquarterperiod;
numpoints = 4*pointsperquarterperiod;
tvec = [0:1:numpoints-1]'*Ts;

% From https://en.wikipedia.org/wiki/Triangle_wave
% (2/a)*|     |  t    1  |    |_ t/a + 1/2 _|
%       |t - a|  -- + -- |(-1)
%       |     |_ a    2 _|
a = period/2;
taplushalf = floor((tvec + opts.shift)/a + 1/2);
y = (2/a)*((tvec + opts.shift)- a*taplushalf);
y = y.*(-1).^taplushalf;

Y_all = [];
for k=1:Nperiod
    Y_all = [Y_all; y];
end
tvec_all = (0:length(Y_all)-1)'*Ts;
rasterdata = timeseries(Y_all, tvec_all);
truefreq = 1/period;
end




