function [transient, prompt, time, xincr, start, fit_start] = ...
    createDecay(tau, nrPhots, noise, promptSigma, nrBins, range)
% [transient, prompt, xincr, time, start, fit_start] = ...
%      createDecay(tau, nrPhots, noise, promptSigma, nrBins, range, offset)
%   Creates single- or multi-exponential decays
%
%       tau         a scalar or vector of lifetime components
%                   [default = 2 ns]
%       nrPhots     number of photons of each lifetime component
%                   [default = 10000]
%       noise       number of background white noise photons
%                   [default = 300]
%       promptSigma sigma of the normal distribution of the prompt
%                   [default = 0.05 ns]
%       nrBins      microtime bin count [default = 256]
%       range       maximum microtime range [default = 10 ns]
%       offset      microtime of transient start [default = 1 ns]
%
%
%       transient   fluorescence lifetime transient
%       prompt      simulated instrument response function
%       time        time histogram bins
%       xincr       time duration of the bin
%       start       bin number of transient rise start
%       fit_start   bin number of the start of fitting of the transient
%
%   Run without any arguments to get the default prompt.
%   Supply one or more arguments to achieve the transient of the required
%   parameters
%
% GNU GPL license 3.0
% copyright 2013 Jakub Nedbal

% parse the inputs and substitute defaults if necessary
if ~exist('tau', 'var')
    tau = 2;
end

if ~exist('nrPhots', 'var')
    nrPhots = 10000;
end

if ~exist('noise', 'var')
    noise = 300;
end

if ~exist('promptSigma', 'var')
    promptSigma = 0.05;
end

if ~exist('nrBins', 'var')
    nrBins = 256;
end

if ~exist('range', 'var')
    range = 10;
end

if ~exist('offset', 'var')
    offset = 1;
end

% check for sensible inputs
if numel(tau) ~= numel(nrPhots)
    warning(['The number of elements in "tau" must match the number ', ...
             'of elements in "nrPhots".\nYou gave %d "tau" but %d ', ...
             '"nrPhots".\n'], numel(tau), numel(nrPhots));
    return
end

if any(tau <= 0)
    warning('"tau" must be positive.\nYou gave %g.\n', min(tau)); %#ok<WNTAG>
    return
end

if any(nrPhots < 0)
    warning('"nrPromts" must be non-negative.\nYou gave %g.\n', ...
            min(nrPrompts)); %#ok<WNTAG>
    return
end

if numel(promptSigma) ~= 1
    warning(['"promptSigma" must be a scalar.\nYour''s contains %d ', ...
             'elements.\n'], numel(promptSigma));
	return
end

if promptSigma < 0
    warning('"promptSigma" must be non-negative.\nYou gave %g.\n', ...
            min(promptSigma)); %#ok<WNTAG>
    return
end

if numel(range) ~= 1
    warning(['"range" must be a scalar.\nYour''s contains %d ', ...
             'elements.\n'], numel(range));
    return
end

if range <= 0
    warning('"range" must be a positive number.\nYou gave %g.\n', range); %#ok<WNTAG>
    return
end

if numel(offset) ~= 1
    warning(['"offset" must be a scalar.\nYour''s contains ', ...
             '%d elements.\n'], numel(offset));
    return
end

if any(offset <= 0 || offset >= range)
    warning(['"offset" must be positive and less than "range".\n', ...
             'You gave %g.\n'], max(offset));
    return
end

    


xincr = range / nrBins;             % nanoseconds per bin
time = linspace(0, range, nrBins);  % histogram bins
% create an single or multi exponential decay
y = zeros(nrBins, 1);
for i = 1 : numel(tau)
    y = y + histc(offset - log(rand(nrPhots(i), 1)) * tau(i), time);
end
% normalize the prompt
prompt = histc(range / 2 + promptSigma * randn(mean(nrPhots), 1), time);
prompt = prompt / sum(prompt);
% transient convolved with the prompt
transient = abs(ifftshift(ifft(fft(y) .* fft(prompt))));
% get rid of extra zeros in the prompt
prompt = prompt(prompt > max(prompt) / 100);
% add white noise
transient = transient + histc(rand(noise, 1) * range, time);
% start of the transient at the point of highest gradient
[dummy, dTrans] = max(gradient(transient));
% minus the point of highest gradient in the prompt
[dummy, dPrompt] = max(gradient(prompt));

% start of the transient at the point of rise of the transient
start = dTrans - dPrompt+1;
% start of the fit is half the prompt width beyond the maximum of the 
% transient
[dummy, fit_start] = max(transient);
fit_start = fit_start + floor(numel(prompt) / 2) - start;
