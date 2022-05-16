%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the RF data for the in-vitro rotating disk                        %
%                                                                        %
% Author: Madhavanunni A N                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sysPara,trueVel] = loadPWIdiscParams(dataPath,datasetName,rfDataType)

load(strcat(dataPath,datasetName));

%%set parametrs
sysPara.txType = 'plane_wave';
sysPara.flowType = strcat(rfDataType,'rotating-disc');
sysPara.txAngles = 0;
sysPara.nElements = param.Nelements;
sysPara.c = param.c;                         % Sound speed [m/s]
sysPara.fc = param.fc;       % Tx center frequency [Hz]
sysPara.lambda = sysPara.c/sysPara.fc;                      % Wavelength [m]
sysPara.elem_pos(:,1) = ((0:127)*param.pitch)-mean((0:127)*param.pitch);% Element positions  [m]
sysPara.elem_pos(:,2:3) = 0;
sysPara.fs = param.fs;                    % Sampling frequency [Hz]
sysPara.focal_delay = param.TXdelay;
sysPara.nSupFrames = 32;
sysPara.nFrames = 1;
sysPara.rxStartTime = param.t0*ones(sysPara.nFrames,1);
sysPara.PRF = param.PRF;
sysPara.no_emission_types = 1;
sysPara.flow_em_per_frame = 1;          % Number of flow emissions per frame
sysPara.txFocus = [0 0 6000];
sysPara.nSamples = size(RF,1);
sysPara.Receive = 'NA';

%-- Load true/theoretical velocity
trueVel.flowRegion = NaN;% linspace(0.028,0.042,1001);
trueVel.theta = 0;
trueVel.v0 = NaN;
trueVel.velVector = NaN;

end
