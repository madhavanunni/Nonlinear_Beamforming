%-------------------------------------------------------------------------%
%-- Script to reproduce Figure 3 of the paper titled:
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
%-- MUST toolbox (https://www.biomecardio.com/MUST/):
%-- ##Functions used from MUST toolbox: vplot and smoothn
%-------------------------------------------------------------------------%
%-- Acknowledgements:
%-- Garcia D. Robust smoothing of gridded data in one and higher dimensions 
%-- with missing values. Computational Statistics & Data Analysis 2010; 54:1167-1178.
%-- Garcia D. A fast all-in-one method for automated post-processing of PIV data. 
%-- Exp Fluids 2011; 50:1247-1259.
%-------------------------------------------------------------------------%
%%
clear
close all;

cLim = [-1.5 1.5];
font_size = 14;

%% Load B-mode data
bmode_data = load('..\data\PWI_disk_Bmode.mat');
addpath('functions\')
addpath(genpath('lib\'))

%% Loading DAS beamformed data and its parameters
load('..\results\PWI_disk_RF\DAS\bfDataParams.mat')
load('..\results\PWI_disk_RF\DAS\bfData_frame_1.mat')
bfMethod = 'DAS';
    
%% Velocity estimation from the beamformed data
zl=length(gridParams.z_axis);
xl=length(gridParams.x_axis);
aCorrLen = 64;
nAlphas = bfParams.nAlphas;
nFrames = size(bfData.samples_left,3);
alpha_left = reshape(txRxAngles.alpha_left,xl,zl,nAlphas);
alpha_left = permute(alpha_left,[3,2,1]);
samples_left = reshape(bfData.samples_left,nAlphas,nFrames,xl,zl);
samples_right = reshape(bfData.samples_right,nAlphas,nFrames,xl,zl);
samples_left = permute(samples_left,[1,4,3,2]);
samples_right = permute(samples_right,[1,4,3,2]);

velocity = getVelocity(samples_left,samples_right,aCorrLen,bfParams.nAlphas,sysPara,alpha_left,bfParams);

[z_new,x_new]=meshgrid(double(gridParams.z_axis),double(linspace(min(gridParams.x_axis),max(gridParams.x_axis),zl)));
Vx = interp2(gridParams.Z,gridParams.X,double(reshape(mean(mean(velocity.Vz,1),3),zl,xl)).',z_new,x_new);
Vz = interp2(gridParams.Z,gridParams.X,double(reshape(mean(mean(velocity.Vx,1),3),zl,xl)).',z_new,x_new);

% Creating a binary mask to select the region of interest (ROI) 
ROImap = (hypot(x_new-(-7.7653e-04),z_new-0.02264)<0.01).'; 

%% Figure 3(a): Color Doppler map obtained with DAS without any smoothing
figure;
Vs = smoothn({-Vx.',-Vz.'},1,'robust');
Vs{1}(~ROImap) = 0;
Vs{2}(~ROImap) = 0;
imagesc((x_new(:,1))*100,(z_new(1,:))*100,-Vs{2});
colormap dopplermap
colorbar
caxis([-0.605 0.605])
set(gca,'FontSize',font_size,'Color','black')
title(strcat('Color Doppler map (',bfMethod,')'))
xlabel('Lateral position [cm]')
ylabel('Axial position [cm]')
axis equal ij tight
xlim([-1.25,1.25])
ylim([1,3.5]);

%% --- Figure 3(b): Vector flow image obtained with DAS 
Vs = smoothn({-Vx.',-Vz.'},5e4,'robust'); % Robust smoothing
Vs{1}(~ROImap) = 0;
Vs{2}(~ROImap) = 0;

figure()
ax1 = axes;
imagesc(bmode_data.x(1,:)*100,bmode_data.z(:,1)*100,squeeze(bmode_data.I(:,:,1)))
colormap(ax1,gray)
axis equal ij tight
set(gca,'FontSize',font_size)
set(ax1,'YLim',[1,3.5],'XLim',[-1.25,1.25])
colorbar('Color','none')
hold on

ax2 = axes;
vplot(x_new.'*100,z_new.'*100,Vs{1},Vs{2})
linkaxes([ax1,ax2]);
colormap(ax2,flipud(hot));%Dopplermap
hold off
title('Motion field (in pix) by speckle tracking')
axis equal tight off ij
set(ax2,'color','none','visible','off');
c = colorbar;
caxis([0 1])
set(ax2,'YLim',[1,3.5],'XLim',[-1.25,1.25])
set(gca,'FontSize',font_size)
title(strcat('Vector Flow Image (',bfMethod,')'))

%% Loading Nonlinear (NLHR) beamformed data and its parameters
load('..\results\PWI_disk_RF\NLHR\bfDataParams.mat')
load('..\results\PWI_disk_RF\NLHR\bfData_frame_1.mat')
bfMethod = 'NLHR';

%% Velocity estimation from the beamformed data
zl=length(gridParams.z_axis);
xl=length(gridParams.x_axis);
aCorrLen = 64;
nAlphas = bfParams.nAlphas;
nFrames = size(bfData.samples_left,3);
alpha_left = reshape(txRxAngles.alpha_left,xl,zl,nAlphas);
alpha_left = permute(alpha_left,[3,2,1]);
samples_left = reshape(bfData.samples_left,nAlphas,nFrames,xl,zl);
samples_right = reshape(bfData.samples_right,nAlphas,nFrames,xl,zl);
samples_left = permute(samples_left,[1,4,3,2]);
samples_right = permute(samples_right,[1,4,3,2]);

velocity = getVelocity(samples_left,samples_right,aCorrLen,bfParams.nAlphas,sysPara,alpha_left,bfParams);

[z_new,x_new]=meshgrid(double(gridParams.z_axis),double(linspace(min(gridParams.x_axis),max(gridParams.x_axis),zl)));
Vx = interp2(gridParams.Z,gridParams.X,double(reshape(mean(mean(velocity.Vz,1),3),zl,xl)).',z_new,x_new);
Vz = interp2(gridParams.Z,gridParams.X,double(reshape(mean(mean(velocity.Vx,1),3),zl,xl)).',z_new,x_new);
velMag = hypot(-Vx,-Vz);
velAng = wrapTo180(atan2d(-Vx,-Vz));

% Creating a binary mask to select the region of interest (ROI) 
ROImap = (hypot(x_new-(-7.7653e-04),z_new-0.02264)<0.01).'; 

%% Figure 3(c): Color Doppler map obtained with NLHR without any smoothing
figure;
Vs = smoothn({-Vx.',-Vz.'},1,'robust');
Vs{1}(~ROImap) = 0;
Vs{2}(~ROImap) = 0;
imagesc((x_new(:,1))*100,(z_new(1,:))*100,-Vs{2});
colormap dopplermap
colorbar
caxis([-0.605 0.605])
set(gca,'FontSize',font_size,'Color','black')
title(strcat('Color Doppler map (',bfMethod,')'))
xlabel('Lateral position [cm]')
ylabel('Axial position [cm]')
axis equal ij tight
xlim([-1.25,1.25])
ylim([1,3.5]);

%% --- Figure 3(d): Vector flow image obtained with NLHR 
%Robust smoothing
Vs = smoothn({-Vx.',-Vz.'},5e4,'robust');
Vs{1}(~ROImap) = 0;
Vs{2}(~ROImap) = 0;

figure()
ax1 = axes;
imagesc(bmode_data.x(1,:)*100,bmode_data.z(:,1)*100,squeeze(bmode_data.I(:,:,1)))
colormap(ax1,gray)
axis equal ij tight
set(gca,'FontSize',font_size)
set(ax1,'YLim',[1,3.5],'XLim',[-1.25,1.25])
colorbar('Color','none')
hold on

ax2 = axes;
vplot(x_new.'*100,z_new.'*100,Vs{1},Vs{2})
linkaxes([ax1,ax2]);
colormap(ax2,flipud(hot));%Dopplermap
hold off
title('Motion field (in pix) by speckle tracking')
axis equal tight off ij
set(ax2,'color','none','visible','off');
c = colorbar;
caxis([0 1])
set(ax2,'YLim',[1,3.5],'XLim',[-1.25,1.25])
set(gca,'FontSize',font_size)
title(strcat('Vector Flow Image (',bfMethod,')'))
