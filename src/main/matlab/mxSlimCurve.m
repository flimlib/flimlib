%MXSLIMCURVE Perform lifetime fit on TCSPC data
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START) Fits a
%   a TCSPC decay with a single-exponential model. It can fit a single
%   TCSPC decay curve or number of them at the same time. TRANSIENT 
%   contains a single or multiple TCSPC decay curves from the start of the 
%   transient to its end. For a single decay TRANSIENT must be in a column
%   vector Mx1, for N different decays it must be an array shape of MxN
%   elements. M are the number of micro time bins in each transient. PROMPT
%   is a normalized instrument response function, i.e. sum(PROMPT) == 1;
%   There may be a single prompt for all decays in a Px1 column vector or
%   each decay may have its own prompt in an PxN array. X_INC is the micro
%   time bin size of the decay. There may be a single X_INC for all decays
%   or a vector of N elements with potentially a different X_INC value for
%   each transient, FIT_START is the index of the TRANSIENT, where fitting
%   should start. Normally, this is the index of the transient's maximum
%   plus half the number of elements in the prompt. LMA_PARAM is the result
%   of the Marquardt-Levenberg Algorithm fitting a single-exponential decay
%   curve to the data in TRANSIENT. It is a 4xN array array taking the form
%   of [Z; A; tau; chi_sq], where:
%       Z:      Noise baseline
%       A:      Decay curve maximum
%       tau:    Lifetime
%       chi_sq: Chi^2 goodness of the fit
%   
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START, FIT_TYPE)
%   By default a single-exponential decay model is assumed. Choose one of
%   the following FIT_TYPE values for different fitting models:
%       1: Single exponential
%       2: Double exponential
%       3: Triple exponential
%       4: Stretched exponential (fitting function GCI_stretchedexp)
%   With different FIT_TYPE the LMA_PARAM gets a different form:
%       FIT_TYPE        LMA_PARAM form
%           1           [Z; A; tau; chi_sq]
%           2           [Z; A1; tau1; A2; tau2; chi_sq]
%           3           [Z; A1; tau1; A2; tau2; A3; tau3; chi_sq]
%           4           [Z; A; tau; H; chi_sq]
%
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START, ...
%                           FIT_TYPE, NOISE_MODEL) By default a Poisson
%   fit-based model for noise is being used. Choose on of the following
%   NOISE_MODEL values for different noise models:
%       0: NOISE_CONST
%       1: NOISE_GIVEN
%       2: NOISE_POISSON_DATA
%       3: NOISE_POISSON_FIT
%       4: NOISE_GAUSSIAN_FIT
%       5: NOISE_MLE
%
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START, ...
%                           FIT_TYPE, NOISE_MODEL, CHI_SQ_TARGET) By 
%   default the fit stops once chi^2 reaches a value of 1.1 or less. Choose
%   CHI_SQ_TARGET greater than one to modify this limit. The resulting fits
%   are likely to be more accurate with smaller values at the expense of
%   longer fitting time.
%
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START, ...
%                           FIT_TYPE, NOISE_MODEL, CHI_SQ_TARGET, ...
%                           CHI_SQ_DELTA) By default the fitting terminates
%   once the chi^2 changes by less than 0.001 after each iterration. Choose
%   a different CHI_SQ_DELTA to modify this cut-off. The resulting fits
%   are likely to be more accurate with smaller values at the expense of
%   longer fitting time.
%
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START, ...
%                           FIT_TYPE, NOISE_MODEL, CHI_SQ_TARGET, ...
%                           CHI_SQ_DELTA, FIT_END) Normaly the entire
%   length of the decay curve from the FIT_START to the end are esed for
%   the fitting. Sometimes the very end is corrupted and it might be
%   desirable to omit the values from the end of the decay curve from the
%   fitting. Choose FIT_END less than the number of elements in TRANSIENT
%   to truncate the section of the decay curve used for the fitting.
%
%   LMA_PARAM = MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START, ...
%                           FIT_TYPE, NOISE_MODEL, CHI_SQ_TARGET, ...
%                           CHI_SQ_DELTA, FIT_END, SIGMA_VALUES) With
%   NOISE_TYPE == 1 (GIVEN) a custom set of standard deviations supplies
%   the individual noise characteristics to each data point. There may be a
%   single sigma array for all decays in a Px1 column vector or each decay
%   may have have its own sigma array in an PxN array.
%
%   [LMA_PARAM, RLD_PARAM, LMA_FIT, RLD_FIT] = ...
%       MXSLIMCURVE(TRANSIENT, PROMPT, X_INC, FIT_START) RLD_param provides
%   the fitted parameters based on Rapid Lifetime Determination in the form
%   of [Z; A; tau; chi_sq]. LMA_FIT provides the fitted decay curve using
%   the Marquardt-Levenberg Algorithm. LMA_FIT provides the fitted decay
%   curve using the Rapid Lifetime Determination Algorithm. 
%
%   EXAMPLE:
%       Create a single exponential decay with 10 nanoseconds duration
%       divided into 256 micro-time points. Create a Gaussian prompt of 
%       50 ps standard deviation. Calculate the start of the transient's
%       rise. Truncate the transient at the beginning. Calculate the start
%       of the fitting and perform the fit.
%
%   % *** Supply decay parameters ***
%   tau = 2;                            % simulated lifetime 2 nanoseconds
%   nrPhots = 1e4;                      % 10000 photons in the decay
%   nrNoise = 500;                      % 500 dark-count background photons
%   range = 10;                         % micro-time range 10 nanoseconds
%   nrBins = 256;                       % number of micro-time bins
%   offset = 1;                         % transient rise time 1 nanosecond
%   pSigma = 0.05                       % prompt sigma is 50 picoseconds
%   time = linspace(0, range, nrBins);  % histogram bins
%   xincr = range / nrBins;             % nanoseconds per bin
%
%   % *** Create transient and prompt (instrument response function) ***
%   % Create a prompt
%   prompt = histc(range / 2 + pSigma * randn(mean(nrPhots), 1), time);
%   prompt = prompt / sum(prompt);
%   % Create a transient
%   transient = histc(offset - log(rand(nrPhots, 1)) * tau, time);
%   % Convolve the transient with the prompt
%   transient = abs(ifftshift(ifft(fft(transient) .* fft(prompt))));
%   % Get rid of extra zeros in the prompt
%   prompt = prompt(prompt > max(prompt) / 1000);
%   % add white noise
%   transient = transient + histc(rand(nrNoise, 1) * range, time);
%
%   % *** Work out start of transient and start of fitting ***
%   % start of the transient at the point of highest gradient
%   [dummy, dTrans] = max(gradient(transient));
%   % minus the point of highest gradient in the prompt
%   [dummy, dPrompt] = max(gradient(prompt));
%   % start of the transient at the point of rise of the transient
%   start = dTrans - dPrompt+1;
%   % start of the fit is half the prompt width beyond the maximum of the 
%   % transient
%   [dummy, fit_start] = max(transient);
%   fit_start = fit_start + floor(numel(prompt) / 2) - start;
%   % truncate the transient before the start of the decay rise
%   transient = transient(start : end);
%
%   % *** Perform the fit ***
%   paramsLMA = mxSlimCurve(transient, prompt, xincr, fit_start);
%   
%   fprintf('A:\t%g\ntau:\t%g\nZ:\t%g\nchi2:\t%g\n', paramsLMA([2 3 1 4]));
%
% GNU GPL license 3.0
% copyright 2013 Jakub Nedbal