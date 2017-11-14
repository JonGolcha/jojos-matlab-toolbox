% This script shows an example of phase randomization of a signal using
% Fourier transform with two shuffling methods (direct or amplitude adjusted).
%
% History:
% --------
% First version: 2014-07-13
% J.Chatel-Goldman @OIST, jonas.chatel.goldman(at)gmail.com


clear all
close all

% generation of a colored (band-pass) signals
x = randn(1,1000);
[Hbp, b, a] = designBandPassFilter_FIR(.15,.25,1);
x_colored_1 = filtfilt(b,a,x);
[Hbp, b, a] = designBandPassFilter_FIR(.25,.35,1);
x_colored_2 = filtfilt(b,a,x);
[Hbp, b, a] = designBandPassFilter_FIR(.35,.45,1);
x_colored_3 = filtfilt(b,a,x);
figure('name','Original signal','WindowStyle','docked'),
subplot(3,1,1),
plot(x_colored_1(1:500));
subplot(3,1,2),
pwelch(x_colored_1,128,96,256);
subplot(3,1,3),
hist(x_colored_1,100),


% generation of surrogate with phase randomization (same power spectrum)
x_shuffle = phaseSurrogate([x_colored_1',x_colored_2',x_colored_3'] ,100000); % this checks if doing alright on multiple signals
x_shuffle_1 = squeeze(x_shuffle(:,1,98400));    % arbitrary signal and shuffling indexes
figure('name','Shuffled signal','WindowStyle','docked'),
subplot(3,1,1),
plot(x_shuffle_1(1:500));
subplot(3,1,2),
pwelch(x_shuffle_1,128,96,256);
subplot(3,1,3),
hist(x_shuffle_1,100),

% generation of surrogate with phase randomization and adjusted amplitude (same power spectrum + same probability distribution)
x_shuffle_2 = phaseran2(x_colored_1,1);
figure('name','Shuffled signal (adjusted amplitude)','WindowStyle','docked'),
subplot(3,1,1),
plot(x_shuffle_2(1:500));
subplot(3,1,2),
pwelch(x_shuffle_2,128,96,256);
subplot(3,1,3),
hist(x_shuffle_2,100),



