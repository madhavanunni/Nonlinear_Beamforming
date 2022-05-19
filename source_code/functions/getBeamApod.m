%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to generate apodization weights for beam apodization              %
%                                                                        %
% This is a modified version of apodization.m function in USTB           %
% UltraSound ToolBox (https://www.ustb.no/)                              %
% Originally developed by:                                               %
%                Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)   %
%                Olivier Bernard (olivier.bernard@creatis.insa-lyon.fr)  %
%                                                                        %
% Modified for channel directive beam synthesis in nonlinear beamforming %
% by: Madhavanunni A N                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function apod = getBeamApod(tof,sysPara,window)

rx_f_number = 1.1248;
rx_aperture = tof/rx_f_number;
elem_idx = 1:sysPara.nElements;
elemIdx = single(ones(size(tof,2),1)*elem_idx).';
apod = zeros([size(tof),sysPara.nElements]);
for nrx=1:sysPara.nElements
    rx_aperture_distance = abs(sysPara.elem_pos(nrx)-sysPara.elem_pos(:,1)).*ones(size(tof));
    
    switch(window)
        case 'none'
            apod(:,:,nrx) = ones(size(rx_aperture_distance));
        case 'box'
            apod(:,:,nrx) = single(rx_aperture_distance<=rx_aperture/2);
        case 'gaussian8'
            apod(:,:,nrx) = single(rx_aperture_distance<=rx_aperture/2).*exp(-(((elemIdx-nrx).^2)/(2*8)));
        case 'hanning'
            apod(:,:,nrx) = single(rx_aperture_distance<=rx_aperture/2).*(0.5 + 0.5*cos(2*pi*rx_aperture_distance./rx_aperture));
        case 'hamming'
            apod(:,:,nrx) = single(rx_aperture_distance<=rx_aperture/2).*(0.53836 + 0.46164*cos(2*pi*rx_aperture_distance./rx_aperture));
        case 'tukey25'
            roll=0.25;
            apod(:,:,nrx) =(rx_aperture_distance<(rx_aperture/2*(1-roll))) + (rx_aperture_distance>(rx_aperture/2*(1-roll))).*(rx_aperture_distance<(rx_aperture/2)).*0.5.*(1+cos(2*pi/roll*(rx_aperture_distance./rx_aperture-roll/2-1/2)));
        case 'tukey50'
            roll=0.5;
            apod(:,:,nrx)=(rx_aperture_distance<(rx_aperture/2*(1-roll))) + (rx_aperture_distance>(rx_aperture/2*(1-roll))).*(rx_aperture_distance<(rx_aperture/2)).*0.5.*(1+cos(2*pi/roll*(rx_aperture_distance./rx_aperture-roll/2-1/2)));
        case 'tukey75'
            roll=0.75;
            apod(:,:,nrx)=(rx_aperture_distance<(rx_aperture/2*(1-roll))) + (rx_aperture_distance>(rx_aperture/2*(1-roll))).*(rx_aperture_distance<(rx_aperture/2)).*0.5.*(1+cos(2*pi/roll*(rx_aperture_distance./rx_aperture-roll/2-1/2)));
        otherwise
            error('Unknown window type. Known types are: box, gaussian8, hamming, hanning, tukey25, tukey50, tukey75.');
    end
    
    
    
end
