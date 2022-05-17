%-------------------------------------------------------------------------%
%-- Important: Please ensure you run this code in MATLAB 2019b or higher
%-------------------------------------------------------------------------%
%-- Script to reproduce Figure 6 of the paper titled:
%-- "A Nonlinear Beamforming for Enhanced Spatiotemporal Sensitivity in 
%-- High Frame Rate Ultrasound Flow Imaging".
%-- Authors: Madhavanunni A N and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 05 - May - 2022
%-------------------------------------------------------------------------%
%-- Dependencies: 
%-- Beamformed data for B-mode and velocity estimation
%-- Time-Frequency Toolbox (https://tftb.nongnu.org/)
%-- ##Functions used from Time-Frequency Toolbox: tfrrpwv
%-------------------------------------------------------------------------%
%-- Acknowledgements:
%-- Franqois Auger and Patrick Flandrin. Improving the Readability of Time-Frequency 
%-- and Time-Scale Representations by the Reassignment Method. 
%-- IEEE Transactions on Signal Processing, 43(5), 1995. 
%-------------------------------------------------------------------------%
%%
clear;
addpath(genpath('lib\'));

%% -- Load DAS beamformed data
PRF = 16000;
font_size = 14;
bfMethod = 'DAS';
startFrame = 35000;
load('..\data\phantom_35_38_DAS_spectrogram_9deg_svd3.mat')

cLim = [-100,0];
qFrame = 1650;
%% -- Signal Plots -- Fig. 6(c) and (d)
hid = figure();
set(gcf,'Position',[50,100,1500,500])
tPlot = tiledlayout(2,4,'TileSpacing','compact');
tPlot.Padding = 'compact';

time_axis = ((qFrame+startFrame)/PRF):1/PRF:((qFrame+startFrame+511)/PRF);
figure(hid);nexttile(1)

plot_sig = real(env_samples_left(:,qFrame));
plot(time_axis,plot_sig/norm(plot_sig),'-b','DisplayName','Left sub-aperture');
title('Left sub-aperture signal')
xlabel('Time [s]')
ylabel('Normalised amplitude')
axis tight
ylim([-0.2,0.2]);
grid
set(gca,'FontSize',font_size)

nexttile(2)
plot_sig = real(env_samples_right(:,qFrame));
plot(time_axis,plot_sig/norm(plot_sig),'-r','DisplayName','Right sub-aperture');
xlabel('Time [s]')
title('Right sub-aperture signal')
axis tight
ylim([-0.45,0.45]);
grid
set(gca,'FontSize',font_size)

%% -- Time-Frequency plots -- Fig. 6(e) and (f)
nexttile(5)
plot_sig = real(env_samples_left(:,qFrame))/norm(real(env_samples_left(:,qFrame)))+...
    imag(env_samples_left(:,qFrame))/norm(imag(env_samples_left(:,qFrame)))*1i;
[s2,w2,t2] = tfrrpwv(plot_sig);
imagesc(time_axis,linspace(0,1,512),20*log10(abs(s2)/max(abs(s2(:)))));
colormap((jet))
caxis(cLim);
set(gca,'FontSize',font_size,'YDir','normal')
ylabel("Normalized Frequency" + newline + "(\times\pi radians/second)")
xlabel('Time [s]')

nexttile(6)
plot_sig = real(env_samples_right(:,qFrame))/norm(real(env_samples_right(:,qFrame)))+...
    imag(env_samples_right(:,qFrame))/norm(imag(env_samples_right(:,qFrame)))*1i;
[s2,w2,t2] = tfrrpwv(plot_sig);
imagesc(time_axis,linspace(0,1,512),20*log10(abs(s2)/max(abs(s2(:)))));
colormap((jet))
caxis(cLim);
set(gca,'FontSize',font_size,'YDir','normal')
xlabel('Time [s]')

%% Load Nonlinear beamformed data
bfMethod = 'NLHR';
startFrame = 35000;
PRF = 16000;

load('..\data\phantom_35_38_NLHR_spectrogram_9deg_svd5.mat')
cLim = [-100,0];
qFrame = 1650;
time_axis = ((qFrame+startFrame)/PRF):1/PRF:((qFrame+startFrame+511)/PRF);

%% -- Signal plots -- Fig. 6(g) and (h)
figure(hid);nexttile(3)
plot_sig = real(env_samples_left(:,qFrame));
plot(time_axis,plot_sig/norm(plot_sig),'-b','DisplayName','Left sub-aperture');
title('Left sub-aperture signal')
xlabel('Time [s]')
axis tight
ylim([-0.2,0.2]);
grid
set(gca,'FontSize',font_size)

nexttile(4)
plot_sig = real(env_samples_right(:,qFrame));
plot(time_axis,plot_sig/norm(plot_sig),'-r','DisplayName','Right sub-aperture');
xlabel('Time [s]')
title('Right sub-aperture signal')
axis tight
ylim([-0.45,0.45]);
grid
set(gca,'FontSize',font_size)

%% -- Time-frequency plots -- Fig. 6(i) and (j)

nexttile(7)
plot_sig = real(env_samples_left(:,qFrame))/norm(real(env_samples_left(:,qFrame)))+...
    imag(env_samples_left(:,qFrame))/norm(imag(env_samples_left(:,qFrame)))*1i;
[s2,w2,t2] = tfrrpwv(plot_sig);
imagesc(time_axis,linspace(0,1,512),20*log10(abs(s2)/max(abs(s2(:)))));
colormap((jet))
caxis(cLim);
set(gca,'FontSize',font_size,'YDir','normal')
xlabel('Time [s]')

nexttile(8)
plot_sig = real(env_samples_right(:,qFrame))/norm(real(env_samples_right(:,qFrame)))+...
    imag(env_samples_right(:,qFrame))/norm(imag(env_samples_right(:,qFrame)))*1i;
[s2,w2,t2] = tfrrpwv(plot_sig);
imagesc(time_axis,linspace(0,1,512),20*log10(abs(s2)/max(abs(s2(:)))));
colormap((jet))
caxis(cLim);
set(gca,'FontSize',font_size,'YDir','normal')
xlabel('Time [s]')
