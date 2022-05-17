%-------------------------------------------------------------------------%
%-- Script for estimating vector velocity using triangulation technique
%-- Beamformers: Delay and Sum, Nonlinear beamforming
%-- "A Nonlinear Beamforming for Enhanced Spatiotemporal Sensitivity in
%-- High Frame Rate Ultrasound Flow Imaging".
%-- Authors: Madhavanunni A N and Mahesh Raveendranatha Panicker
%-- Affiliation: Indian Institute of Technology Palakkad, India
%-------------------------------------------------------------------------%
%-- Version: v1.0
%-- Last modified on 05 - May - 2022
%-------------------------------------------------------------------------%
%-- Dependencies:
%-- MUST Toolbox(https://www.biomecardio.com/MUST/index.html) for the 
%-- in-vitro rotating disk dataset, robust smoothing (smoothn function) and 
%-- vector flow visualisation (vplot function): 
%-------------------------------------------------------------------------%
%-- Acknowledgements for the rotating disk dataset, smoothn and vplot functions :
%-- 1. Damien Garcia. Make the most of MUST, an open-source Matlab UltraSound Toolbox
%-- IEEE International Ultrasonics Symposium, IUS, 2021. 
%-- 2. Craig Madiena, et al., Color and Vector Flow Imaging in Parallel Ultrasound with 
%-- Sub-Nyquist Sampling IEEE TUFFC, 65(5):795–802, 2018. 
%-- 3. Garcia D. Robust smoothing of gridded data in one and higher dimensions 
%-- with missing values. Computational Statistics & Data Analysis 2010; 54:1167-1178.
%-- 4. Garcia D. A fast all-in-one method for automated post-processing of PIV data
%-- Exp Fluids 2011; 50:1247-1259.
%-------------------------------------------------------------------------%
%-- Set bfParams.DMAS = 0 for DAS beamforming
%-- Set bfParams.beamApod = 1 and bfParams.DMAS = 1 to reproduce the
%-- results for nonlinear beamforming reported in our article.
%-------------------------------------------------------------------------%
clear 
%% Set all the required paths

addpath('functions');
addpath(genpath('lib\'));

saveEnable = 1; % Set saveEnable = 1 if the results have to be saved
%% Load RF data and get the required parameters
dataPath = "..\data\";
fileName = "dataset_frame_";
rfDataType = 'phantom';
datasetName = "PWI_disk_RF";

[sysPara,trueVel] = loadPWIdiscParams(dataPath,datasetName,rfDataType);

%% Select the beamformer parameters
% bfParams.DMAS = 0 for DAS beamforming
% bfParams.DMAS = 1 for Nonlinear beamforming
% bfParams.beamApod = 0 for single stage apodization
% bfParams.beamApod = 1 for double stage apodization
% bfParams.DMAStype = 0 for conventional architecture
% bfParams.DMAStype = 1 for simplified architecture
bfParams.beamApod = 1;
bfParams.DMAS = 1;
bfParams.DMAStype = 1; % Donot alter this line: "bfParams.DMAStype = 1"

%% Define the beamforming parameters
autoCorrLength = 64;                         % Length of an ensemble for autocorelation
nFramesToBF = sysPara.nSupFrames;
nEnsemble = 1;                               % Number of ensembles
bfParams.alphaSet = single(6:3:15);          % Number of Tx-Rx angles to beamform
bfParams.nAlphas = size(bfParams.alphaSet,2);
PRF = sysPara.PRF;

if bfParams.beamApod==0 && bfParams.DMAS==0
    bfParams.method = "DAS";
elseif bfParams.beamApod==1 && bfParams.DMAS==0
    bfParams.method = "DAS-DualApod";
elseif bfParams.beamApod==0 && bfParams.DMAS==1 && bfParams.DMAStype == 0
    bfParams.method = "DMAS";
elseif bfParams.beamApod==1 && bfParams.DMAS==1 && bfParams.DMAStype == 0
    bfParams.method = "DMAS-DualApod";
elseif bfParams.beamApod==0 && bfParams.DMAS==1 && bfParams.DMAStype == 1
    bfParams.method = "simple DMAS";
elseif bfParams.beamApod==1 && bfParams.DMAS==1 && bfParams.DMAStype == 1
    bfParams.method = "NLHR";
end

outDirName = strcat('..\results\',datasetName,filesep,bfParams.method);

if saveEnable && (~exist(outDirName,'dir'))
    mkdir(outDirName);
end

%% Define the grid parameters
xAxisOrg = single(sysPara.elem_pos(:,1));
timeVector = single(((0:(sysPara.nSamples-1))/sysPara.fs).'+sysPara.rxStartTime(1));
timeVectorOrg = timeVector;
zAxisOrg = single(0.5*sysPara.c*timeVector);
gridParams.theta = 0; % Select the grid rotation angle for directional beamforming

% Choose the ROI --- x_axis and z_axis to beamform
gridParams.x_axis = xAxisOrg(xAxisOrg>-0.012 & xAxisOrg<0.01).';
gridParams.x_axis = gridParams.x_axis(1:1:end);
sysPara.fs = 2*sysPara.fs;
timeVector = single(((0:(2*sysPara.nSamples-1))/(sysPara.fs)).'+sysPara.rxStartTime(1));
zAxisOrg = 0.5*sysPara.c*timeVector;
gridParams.z_axis(1,:) = zAxisOrg(zAxisOrg>0.012 & zAxisOrg<0.033).';
gridParams.zResolution = mean(gradient(gridParams.z_axis));
trueVel.Mag = NaN(size(gridParams.z_axis));
trueVel.Ang = ones(size(gridParams.z_axis))*trueVel.theta;


%% Load the estimation points grid create the grid (grid for conventional BF)
zl=length(gridParams.z_axis);
xl=length(gridParams.x_axis);

[~,X] = meshgrid(gridParams.z_axis,gridParams.x_axis);
Z = repmat(gridParams.z_axis,[xl 1])+ X.*tand(-gridParams.theta);
est_points = [X(:), zeros(size(X(:))), Z(:)];

% Generate the beamforming points
bf_points = single(est_points);
clear est_points

%% Transmit-Recieve elements and Alpha estimation for all the emissions

%-- Get the transmit elements' index for plane wave
[txElements, vectors] = getTxElements_PW(bf_points,sysPara,sysPara.txAngles);

%-- Alpha estimation
[alphas, vectors, rxElements, txRxAngles] = getAlphas(sysPara,bf_points,bfParams.alphaSet,txElements,vectors);

%% Focal delay
focal_delay = sysPara.focal_delay;
%% Load the frame data
for loopIdx = 1:sysPara.nFrames
    if saveEnable && exist(strcat(outDirName,'\','bfData_frame_',num2str(loopIdx),'.mat'),'file')
        continue
    else
        load(strcat(dataPath,datasetName));
        rawData = RF;
        %%Resample the RF data in frame direction to avoid aliasing
        clear rfDataset
        rfDataset = zeros(size(rawData,1),size(rawData,2),size(rawData,3)*2);
        for ch = 1:sysPara.nElements
            rfDataset(:,ch,:) = resample(squeeze(rawData(:,ch,:)).',2,1).';
        end
        bfParams.autoCorrLength = autoCorrLength;
        bfParams.nFrames = sysPara.nSupFrames*2;
        bfParams.nEnsemble = bfParams.nFrames-autoCorrLength+1; % Ensemble length for resampled data
        sysPara.PRF = 2*PRF;
        
        clear rawData
        %%-------- Beamform the ensemble ----------- %%
        it_frames=1;
        bf_etypes = sysPara.no_emission_types;
        nFramesToBF = bfParams.nFrames;
        nAlphas = bfParams.nAlphas;
        clear bfData
        bfData.samples_left = single(zeros(length(bf_etypes),size(txRxAngles.alpha_left,3),nFramesToBF,size(bf_points,1)));
        bfData.samples_right = single(zeros(length(bf_etypes),size(txRxAngles.alpha_left,3),nFramesToBF,size(bf_points,1)));
        tic
        % For each emission type
        for  i=1:length(bf_etypes)
            it_em_type = bf_etypes(i);
            em_idx=it_em_type;
            
            % For each emission in the ensemble
            wb = waitbar(0,"DAS beamforming for Tx No: "+it_em_type);
            
            for frame=1:nFramesToBF
                waitbar(frame/nFramesToBF,wb,sprintf('Beamforming the super frame %d/%d in frame %d',frame,nFramesToBF,loopIdx));
                %  Load the samples for the current frame
                samples = squeeze(rfDataset(:,:,frame));
                
                % Remove DC component (from hardware if present)
                rfData = (samples)-repmat(mean(samples,1),[size(samples,1) 1]);
                
                % Resample in depth dierction to avoid aliasing
                rfData = single(resample(rfData,2,1));
                
                
                % Time of flight calculation
                tx_d = repmat(vectors.tx_distance(it_em_type,:),[size(sysPara.elem_pos,1),1]);
                vectors.tof_total = ((vectors.rx_distance+tx_d)/sysPara.c);
                timeVector = single(((0:(size(rfData,1)-1))/sysPara.fs).'+sysPara.rxStartTime(1));
                
                % Interpolation
                tof_samples = zeros(sysPara.nElements,xl*zl);
                for nrx=1:sysPara.nElements
                    tof_samples(nrx,:) = interp1(timeVector+focal_delay(nrx),rfData(:,nrx),vectors.tof_total(nrx,:),'spline',0);
                end
                
                % Beam Apodization
                if bfParams.beamApod==1
                    if ~exist('beamApod','var')
%                         beamApod = getBeamApod(vectors.rx_distance,sysPara,'gaussian8');
                        load('..\data\beamApod_gauss8.mat')
                    end
                    %beamApod_samples = squeeze(sum((repmat(tof_samples,1,1,128).*beamApod),1)).';
                    beamApod_samples = zeros(size(tof_samples));
                    for chIdx = 1:size(beamApod,3)
                        beamApod_samples(chIdx,:) = squeeze(sum(tof_samples.*beamApod(:,:,chIdx),1)).';
                    end
                else
                    beamApod_samples = tof_samples;
                end
                
                % Beamformer
                [samples_left, samples_right,dualApodLeft, dualApodRight] = simpleMASBeamformer(beamApod_samples,it_em_type,sysPara,bf_points,nAlphas,txElements,rxElements,txRxAngles,xl,zl,timeVector,gridParams,bfParams.DMAS);

                bfData.samples_left(i,:,frame,:) = samples_left;
                bfData.samples_right(i,:,frame,:) = samples_right;
                
            end
            close(wb);
        end
        
        clc;
        disp(strcat('Beamforming completed for: ',num2str(loopIdx),'/ ',num2str(sysPara.nFrames)))
        toc
        clear rawData
        if saveEnable
            save(strcat(outDirName,'\','bfData_frame_',num2str(loopIdx),'.mat'),'bfData');
            disp(strcat('Beamformed data is saved.'))
        end
    end
end
%% saving parameters
sysPara = rmfield(sysPara,'Receive');
gridParams.zAxisOrg = zAxisOrg;
gridParams.xAxisOrg = xAxisOrg;
gridParams.X = X;
gridParams.Z = Z;
if saveEnable
    save(strcat(outDirName,'\','bfDataParams','.mat'),'dataPath','datasetName',...
        'gridParams','bfParams','sysPara','timeVector','txRxAngles','bf_points',...
        'trueVel','vectors','rxElements','txElements');
    disp(strcat('Parameter matrices saved.'))
else
    disp('Saving skipped.')
end

%% Show results
%%-- Load B-mode data
bmode_data = load('..\data\PWI_disk_Bmode.mat');

bfMethod = bfParams.method;
cLim = [-1.5 1.5];
font_size = 14;

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

%% --- Color Doppler map
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

%% --- Vector flow image 
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
