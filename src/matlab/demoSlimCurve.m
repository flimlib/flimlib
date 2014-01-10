function demoSlimCurve
% demoSlimCurve     This script runs SlimCurve to demonstrate its
%                   functions.
%
%   demoSlimCurve creates a figure with four panels:
%
%       The first panel shows default single-exponential fit
%
%       The second panel shows a double-exponential transient fitted with a
%       single-exponential model.
%
%       The third panel shows a double-exponential transient fitted with a
%       double-exponential model.
%
%       The fourth-panel shows three histograms of 10000 fits each with
%       varying number of photons in each transient. They include 500,
%       2000 or 10000 photons. The more photons are in a transient, the
%       the accurate the fit becomes.
%
% GNU GPL license 3.0
% copyright 2013 Jakub Nedbal

%% First make sure that the mxSlimCurve file exists.
if ~exist(['mxSlimCurve.' mexext], 'file')
    try
        compileSlimCurve;
    catch %#ok<CTCH>
        error('Could not compile SlimCurve.');
    end
end
    
%% Create a figure
close all
figure('Position', [0 0 640, 800])

%% Create, plot and annotate a single exponential decay
[transient, prompt, time, xincr, start, fit_start] = createDecay;
fit_end = numel(transient) - start - 1;
[paramsLMA, paramsRLD, fittedLMA] = ...
    mxSlimCurve(transient(start : end), prompt, xincr, fit_start, 1, 3);

% plot the curve and its description
createPlot(1, time, start, fit_start, fit_end, transient, paramsLMA, ...
           fittedLMA, 'Basic mono-exponential fit');
% annotate it
text(0.40, 0.8, '# h\nu:', 'FontSize', 16);
text(0.40, 0.6, '\tau:', 'FontSize', 16);
text(0.6, 0.8, '10 000', 'FontSize', 16);
text(0.6, 0.6, '2.00 ns', 'FontSize', 16);

%% Create, plot and annotate a double-exponential decay with
%% single-exponential fit
[transient, prompt, time, xincr, start, fit_start] = ...
    createDecay([3, 1.5], [2000000, 1000000], 20000, 0.05, 2048);
fit_end = numel(transient) - start - 1;

[paramsLMA, paramsRLD, fittedLMA] = ...
    mxSlimCurve(transient(start : end), prompt, xincr, fit_start, 1, 3);

% plot the curve and its description
createPlot(2, time, start, fit_start, fit_end, transient, paramsLMA, ...
           fittedLMA, '2-exponential decay, 1-exponential fit');
% annotate it
text(0.40, 0.72, '# h\nu_1:', 'FontSize', 16);
text(0.40, 0.52, '\tau_1:', 'FontSize', 16);
text(0.40, 0.32, '# h\nu_2:', 'FontSize', 16);
text(0.40, 0.12, '\tau_2:', 'FontSize', 16);
text(0.65, 0.8, '200 000', 'FontSize', 16);
text(0.65, 0.6, '3.00 ns', 'FontSize', 16);
text(0.65, 0.4, '100 000', 'FontSize', 16);
text(0.65, 0.2, '1.50 ns', 'FontSize', 16);

%% Create, plot and annotate a double-exponential decay with
%% double-exponential fit
[paramsLMA, paramsRLD, fittedLMA] = ...
    mxSlimCurve(transient(start : end), prompt, xincr, fit_start, 2);

% plot the curve and its description
createPlot(3, time, start, fit_start, fit_end, transient, paramsLMA, ...
           fittedLMA, '2-exponential decay, 2-exponential fit');
% annotate it
text(0.40, 0.72, '# h\nu_1:', 'FontSize', 16);
text(0.40, 0.52, '\tau_1:', 'FontSize', 16);
text(0.40, 0.32, '# h\nu_2:', 'FontSize', 16);
text(0.40, 0.12, '\tau_2:', 'FontSize', 16);
text(0.65, 0.8, '200 000', 'FontSize', 16);
text(0.65, 0.6, '3.00 ns', 'FontSize', 16);
text(0.65, 0.4, '100 000', 'FontSize', 16);
text(0.65, 0.2, '1.50 ns', 'FontSize', 16);

%% Run single exponential fit 10000x with three different photon counts in
%  each fit.
axes('Position', [0.08, 1.05 - 4 * 0.25, 0.5, 0.15]);
hbins = 1.5 : 0.025 : 3;
hists = zeros(4, numel(hbins));
nrPhots = [500, 2000, 10000];
nrDecays = 10000;
nrNoise = nrPhots / 20;
tau = 2.25;
fits = zeros(2, numel(nrPhots));

for i = 1 : numel(nrPhots)
    [transient, prompt, time, xincr, start, fit_start] = ...
        createManyDecays(nrDecays, tau, nrPhots(i), nrNoise(i));
    [paramsLMA] = ...
        mxSlimCurve(transient(start : end, :), prompt, xincr, fit_start);
    hists(i, :) = histc(paramsLMA(3, :), hbins);
    [fits(1, i), fits(2, i)] = ...
        roundtoerror(mean(paramsLMA(3, :)), std(paramsLMA(3, :)));
end
plot(hbins, hists, 'LineWidth', 2)
hold on
set(gca, 'FontSize', 16, 'YTick', [])
title('Lifetime Distribution vs. Photon Count')
legend({'500', '2000', '10000'});
plot([tau, tau], get(gca, 'YLim'), 'k:')
axes('Position', [0.60, 1.05 - 4 * 0.25, 0.35, 0.15], 'Visible', 'off');
hold on

text(0.02, 0.8, '\tau:', 'FontSize', 16);
text(0.27, 0.8, num2str(tau), 'FontSize', 16);
for i = 1 : numel(nrPhots)
    text(0.02, 0.72 - 0.2 * i, sprintf('\\tau_{%d}:', nrPhots(i)), ...
        'FontSize', 16);
    text(0.27, 0.8 - 0.2 * i, sprintf('%g\\pm%g', fits(:, i)), ...
        'FontSize', 16);
end

%% Export the figure
set(gcf, 'PaperUnits', 'centimeters', 'PaperSize', [21, 29.7], ...
    'PaperPosition', [0, 0, 21, 29.7], 'Units', 'centimeters')
saveas(gcf, 'demoSlimCurve', 'pdf')

function createPlot(in, time, start, fit_start, fit_end, transient, ...
                    paramsLMA, fittedLMA, tit)
% Create a plot with the simulated and fitted data

axes('Position', [0.08, 1.05 - in * 0.25, 0.5, 0.15]);
hold on
yl = [0 1.2 * max(transient)];
fill(time([1 start start 1 1]), [0 0 yl(2) yl(2) 0], [0.9 0.9 0.9], ...
     'EdgeColor', 'none')
plot(time(1:start), transient(1:start), 'Color', [0.5 0.5 0.5], ...
     'LineWidth', 2)
plot(time(start:start+fit_end), transient(start:start+fit_end), ...
     'Color', 'b', 'LineWidth', 2)
plot(time(start+fit_end:end), transient(start+fit_end:end), ...
     'Color', [0.5 0.5 0.5], 'LineWidth', 2)
plot(time(start * [1 1]), yl, 'k:', 'LineWidth', 2)
text(time(start), yl(2) / 3, 'start', 'FontSize', 16, 'Rotation', 90, ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
plot(time((start + fit_start) * [1 1]), yl, 'k:', 'LineWidth', 2)
text(time(start + fit_start), yl(2) / 3, 'fit\_start', 'FontSize', 16, ...
     'Rotation', 90, 'VerticalAlignment', 'cap', ...
     'HorizontalAlignment', 'center')
plot(time((start + fit_end) * [1 1]), yl, 'k:', 'LineWidth', 2)
text(time(start + fit_end), yl(2) / 3, 'fit\_end', 'FontSize', 16, ...
     'Rotation', 90, 'VerticalAlignment', 'bottom', ...
     'HorizontalAlignment', 'center')
set(gca, 'YLim', yl, 'YTick', [], 'FontSize', 16)
%xlabel('Microtime [ns]')
%ylabel('Photon count')
title(tit)

plot(time(start : start + fit_end), fittedLMA(1 : fit_end + 1, 1)', 'r-')

axes('Position', [0.60, 1.05 - in * 0.25, 0.35, 0.15], 'Visible', 'off');
hold on
if size(paramsLMA, 1) == 4
    text(0.02, 0.8, 'A:', 'FontSize', 16);
    text(0.02, 0.6, '\tau:', 'FontSize', 16);
    text(0.02, 0.4, 'Z:', 'FontSize', 16);
    text(0.02, 0.2, '\chi^2:', 'FontSize', 16);
    text(0.14, 0.8, num2str(round(paramsLMA(2))), 'FontSize', 16);
    text(0.14, 0.6, num2str(round(100 * paramsLMA(3)) / 100), ...
        'FontSize', 16);
    text(0.14, 0.4, num2str(round(100 * paramsLMA(1)) / 100), ...
        'FontSize', 16);
    text(0.14, 0.2, num2str(0 + round(100 * paramsLMA(4)) / 100), ...
        'FontSize', 16);
elseif size(paramsLMA, 1) == 6
    text(0.02, 0.88, 'A_1:', 'FontSize', 14);
    text(0.02, 0.72, '\tau_1:', 'FontSize', 14);
    text(0.02, 0.56, 'A_2:', 'FontSize', 14);
    text(0.02, 0.40, '\tau_2:', 'FontSize', 14);
    text(0.02, 0.25, 'Z:', 'FontSize', 14);
    text(0.02, 0.08, '\chi^2:', 'FontSize', 14);
    text(0.14, 0.93, num2str(round(paramsLMA(2))), 'FontSize', 14);
    text(0.14, 0.77, num2str(round(100 * paramsLMA(3)) / 100), ...
        'FontSize', 14);
    text(0.14, 0.61, num2str(round(paramsLMA(4))), 'FontSize', 14);
    text(0.14, 0.45, num2str(round(100 * paramsLMA(5)) / 100), ...
        'FontSize', 14);
    text(0.14, 0.25, num2str(round(100 * paramsLMA(1)) / 100), ...
        'FontSize', 14);
    text(0.14, 0.06, num2str(0 + round(100 * paramsLMA(6)) / 100), ...
        'FontSize', 14);
end

%% Print out the figure
print(gcf, 'demoSlimCurve.pdf', '-dpdf')

function [transient, prompt, time, xincr, start, fit_start] = ...
    createManyDecays(nrDecays, tau, nrPhots, nrNoise)
transient = zeros(256, nrDecays);
xincr = 10 / 256;     % nanoseconds per bin
time = linspace(0, 10, 256);
% create and normalize the prompt
prompt = histc(5 + 0.05 * randn(nrPhots, 1), time);
prompt = prompt / sum(prompt);
% create transients, convolve with prompt, add white noise
for i = 1 : nrDecays
    y = histc(1 - log(rand(nrPhots, 1)) * tau, time);
    y = abs(ifftshift(ifft(fft(y) .* fft(prompt))));
    transient(:, i) = y + histc(rand(nrNoise, 1) * 10, time);
end
% get rid of extra zeros in the prompt
prompt = prompt(prompt > max(prompt) / 100);
% start of the transient at the point of highest gradient
[dummy, dTrans] = max(gradient(mean(transient, 2)));
% minus the point of highest gradient in the prompt
[dummy, dPrompt] = max(gradient(prompt));

% start of the transient at the point of rise of the transient
start = dTrans - dPrompt+1;
% start of the fit is half the prompt width beyond the maximum of the 
% transient
[dummy, fit_start] = max(mean(transient, 2));
fit_start = fit_start + floor(numel(prompt) / 2) - start;


function [Vo, Eo] = roundtoerror(V, E)
%ROUNDTOERROR Rounds the values and their errors to a sensible number of 
%   significant digits. The error will be one significant digit or two
%   significant digits if between 1.0 and 1.9e-?. The value will be rounded
%   accordingly. V and E must be the same size.
%
%	ROUNDTOERROR(V, E) rounds the values V and E to the number of
%   significant digits matching the scale of the error.
%
%	Examples:
%		roundtoerror(234.3217, 0.5468) returns 234.3    and 0.5
%		roundtoerror(234.3217, 0.1345) returns 234.32   and 0.13
%		roundtoerror(234.3217, 0)      returns 234.3217 and 0
%

% check what decimal point is the first valid number
rounding = 10 .^ (floor(log10(abs(E))) - 1);
y = round(E ./ rounding);
% select those which start with 1
one = double(y > 9 & y < 20);
% round the result
rounding = 10 .^ (floor(log10(abs(E))) - one);
Eo = round(E ./ rounding) .* rounding;
Vo = round(V ./ rounding) .* rounding;
% make sure that undesirable NaNs are not introduced instead of zeros
inE = E == 0;
Eo(inE) = 0;
Vo(V == 0) = 0;
Vo(inE) = V(inE);   % keep the initial value if error is zero