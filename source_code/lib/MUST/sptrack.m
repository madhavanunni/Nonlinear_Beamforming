function [di,dj,ic,jc] = sptrack(I,param)

%SPTRACK   Speckle tracking using Fourier-based cross-correlation
%   [Di,Dj] = SPTRACK(I,PARAM) returns the motion field [Di,Dj] that occurs
%   from frame#k I(:,:,k) to frame#(k+1) I(:,:,k+1).
%
%   I must be a 3-D array, with I(:,:,k) corresponding to image #k. I can
%   contain more than two images (i.e. size(I,3)>2). In such a case, an
%   ensemble correlation is used.
%
%   Try it: enter "sptrack" in the command window for an example.
%
%   Di,Dj are the displacements (unit = pix) in the IMAGE coordinate system
%   (i.e. the "matrix" axes mode). The i-axis is vertical, with values
%   increasing from top to bottom. The j-axis is horizontal with values
%   increasing from left to right. The coordinate (1,1) corresponds to the
%   center of the upper left pixel.
%   To display the displacement field, you may use: quiver(Dj,Di), axis ij
%
%   PARAM is a structure that contains the parameter values required for
%   speckle tracking (see below for details).
%
%   [Di,Dj,id,jd] = SPTRACK(...) also returns the coordinates of the points
%   where the components of the displacement field are estimated.
%
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%   1) PARAM.winsize: Size of the interrogation windows (REQUIRED)
%           PARAM.winsize must be a 2-COLUMN array. If PARAM.winsize
%           contains several rows, a multi-grid, multiple-pass interro-
%           gation process is used.
%           Examples: a) If PARAM.winsize = [64 32], then a 64-by-32
%                        (64 lines, 32 columns) interrogation window is used.
%                     b) If PARAM.winsize = [64 64;32 32;16 16], a 64-by-64
%                        interrogation window is first used. Then a 32-by-32
%                        window, and finally a 16-by-16 window are used.
%   2) PARAM.overlap: Overlap between the interrogation windows
%                     (in %, default = 50)
%   3) PARAM.iminc: Image increment (for ensemble correlation, default = 1)
%            The image #k is compared with image #(k+PARAM.iminc):
%            I(:,:,k) is compared with I(:,:,k+PARAM.iminc)
%   5) PARAM.ROI: 2-D region of interest (default = the whole image).
%            PARAM.ROI must be a logical 2-D array with a size of
%            [size(I,1),size(I,2)]. The default is all(isfinite(I),3).
%
%   NOTES:
%   -----
%   The displacement field is returned in PIXELS. Perform an appropriate
%   calibration to get physical units.
%
%   SPTRACK is based on a multi-step cross-correlation method. The SMOOTHN
%   function (see Reference below) is used at each iterative step for the
%   validation and post-processing.
%
%
%   Example:
%   -------
%   I1 = conv2(rand(500,500),ones(10,10),'same'); % create a 500x500 image
%   I2 = imrotate(I1,-3,'bicubic','crop'); % clockwise rotation
%   param.winsize = [64 64;32 32];
%   [di,dj] = sptrack(cat(3,I1,I2),param);
%   quiver(dj(1:2:end,1:2:end),di(1:2:end,1:2:end))
%   axis equal ij
%
%
%   References for speckle tracking 
%   -------------------------------
%   1) Garcia D, Lantelme P, Saloux É. Introduction to speckle tracking in
%   cardiac ultrasound imaging. Handbook of speckle filtering and tracking
%   in cardiovascular ultrasound imaging and video. Institution of
%   Engineering and Technology. 2018.
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/eti18.pdf')">PDF download</a>
%   2) Perrot V, Garcia D. Back to basics in ultrasound velocimetry:
%   tracking speckles by using a standard PIV algorithm. IEEE International
%   Ultrasonics Symposium (IUS). 2018
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/ius18.pdf')">PDF download</a>
%
%   References for smoothing 
%   -------------------------------
%   1) Garcia D, Robust smoothing of gridded data in one and higher
%   dimensions with missing values. Computational Statistics & Data
%   Analysis, 2010.
%   <a
%   href="matlab:web('http://www.biomecardio.com/pageshtm/publi/csda10.pdf')">PDF download</a>
%   2) Garcia D, A fast all-in-one method for automated post-processing of
%   PIV data. Experiments in Fluids, 2011.
%   <a
%   href="matlab:web('http://www.biomecardio.com/pageshtm/publi/expfluids10.pdf')">PDF download</a>
%   
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also SMOOTHN
%
%   -- Damien Garcia & Vincent Perrot -- 2013/02, last update: 2021/01/29
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>



%------------------------%
% CHECK THE INPUT SYNTAX %
%------------------------%

% Run the example if no input
if nargin == 0   
	RunTheExample;
    return
end

narginchk(2,2);
I = double(I);

% Image size
assert(ndims(I)==3,'I must be a 3-D array')
[M,N,P] = size(I);

% Turn off warning messages for SMOOTHN
warn01 = warning('off','MATLAB:smoothn:MaxIter'); 
warn02 = warning('off','MATLAB:smoothn:SLowerBound');
warn03 = warning('off','MATLAB:smoothn:SUpperBound');

if ~isfield(param,'winsize') % Sizes of the interrogation windows
    error(['Window size(s) (PARAM.winsize) ',...
        'must be specified in the field PARAM.'])
end
assert(ismatrix(param.winsize) && size(param.winsize,2)==2,...
    'PARAM.winsize must be a 2-column array.');
tmp = diff(param.winsize,1,1);
assert(all(tmp(:)<=0),...
    'The size of interrogation windows (PARAM.winsize) must decrease.')

if ~isfield(param,'overlap') % Overlap
    param.overlap = 50;
end
assert(isscalar(param.overlap) && param.overlap>=0 && param.overlap<100,...
    'PARAM.overlap (in %) must be a scalar in [0,100[.')
overlap = param.overlap/100;

if ~isfield(param,'ROI') % Region of interest
    ROI = all(isfinite(I),3);
else
    ROI = param.ROI;
    assert(islogical(ROI) & isequal(size(ROI),[M N]),...
        'PARAM.ROI must be a binary image the same size as I(:,:,1).')
end

I(repmat(~ROI,[1 1 P])) = NaN; % NaNing outside the ROI

if ~isfield(param,'iminc') % Step increment
    param.iminc = 1;
end
assert(param.iminc>0 && param.iminc==round(param.iminc),...
    'PARAM.iminc must be a positive integer.')
assert(param.iminc<P,...
    'PARAM.iminc must be < size(I,3).')


i = []; j = [];
m = []; n = [];
options.TolZ = .1; % will be used in SMOOTHN

for kk = 1:size(param.winsize,1)
    
    % Window centers
    ic0 = (2*i+m-1)/2;
    jc0 = (2*j+n-1)/2;
    
    % Size of the interrogation window
    m = param.winsize(kk,1);
    n = param.winsize(kk,2);
    
    % Positions (row,column) of the windows (left upper corner)
    inci = ceil(m*(1-overlap)); incj = ceil(n*(1-overlap));
    [j,i] = meshgrid(1:incj:N-n+1,1:inci:M-m+1);
    
    % Size of the displacement-field matrix
    siz = floor([(M-m)/inci (N-n)/incj])+1;
    
    % Window centers
    ic = (2*i+m-1)/2;
    jc = (2*j+n-1)/2;
    
    if kk>1
        % Interpolation onto the new grid
        di = interp2(jc0,ic0,di,jc,ic,'*cubic');
        dj = interp2(jc0,ic0,dj,jc,ic,'*cubic');
        
        % Extrapolation (remove NaNs)
        dj = rmnan(di+1i*dj,2);
        di = real(dj); dj = imag(dj);
        di = round(di); dj = round(dj);
    else
        di = zeros(siz);
        dj = di;
    end
    
    % Hanning window
    H = hanning(n)'.*hanning(m);
    
    
    C = zeros(siz); % will contain the correlation coefficients
    
    for k = 1:numel(i)
        
        %-- Split the images into small windows
        if i(k)+di(k)>0 && j(k)+dj(k)>0 &&...
                i(k)+di(k)+m-2<M && j(k)+dj(k)+n-2<N
            I1w = I(i(k):i(k)+m-1,j(k):j(k)+n-1,1:end-param.iminc,:);
            I2w = I(i(k)+di(k):i(k)+di(k)+m-1,...
                j(k)+dj(k):j(k)+dj(k)+n-1,1+param.iminc:end,:);
        else
            di(k) = NaN; dj(k) = NaN;
            continue
        end
        
        if any(isnan(I1w(:)+I2w(:)))
            di(k) = NaN; dj(k) = NaN;
            continue
        end
        
        %-- FFT-based cross-correlation
        R = fft2((I2w-mean(I2w(:))).*H).*...
            conj(fft2((I1w-mean(I1w(:))).*H));
        %- 
        R2 = R./(abs(R)+eps); % normalized x-corr
        R2 = mean(R2,3);
        R2 = ifft2(R2,'symmetric');
        C(k) = max(R2(:));
        % C = correlation coefficients: will be used as weights in SMOOTHN
        %--
        R = sum(R,3); % ensemble correlation
        R = ifft2(R,'symmetric'); % x-correlation
        
        %-- Peak detection + Displacement
        %- Circular shift:
        % R = circshift(R,floor(sizR/2)); % is too slow! Better doing this:
        m1 = floor(m/2); n1 = floor(n/2);
        R = R([m-m1+1:m 1:m-m1],[n-n1+1:n 1:n-n1]);
        % --
        % [~,idx] = max(R,[],'all','linear');
        [~,idx] = max(R(:));
        di0 = mod(idx-1,m)+1; % line number
        dj0 = (idx-di0)/m+1; % column number

        %-- Subpixel motion estimation (parabolic peak fitting)
        if di0>1 && di0<m && dj0>1 && dj0<n
            tmp = di0 + (R(di0-1,dj0)-R(di0+1,dj0))/...
                (2*R(di0-1,dj0)-4*R(di0,dj0)+2*R(di0+1,dj0));
            dj0 = dj0 + (R(di0,dj0-1)-R(di0,dj0+1))/...
                (2*R(di0,dj0-1)-4*R(di0,dj0)+2*R(di0,dj0+1));
            di0 = tmp;
        end
        
        % Correction for "circshift(R,floor(sizR/2))"
        di0 = di0-m1-1;
        dj0 = dj0-n1-1;
        
        %-- Total displacement
        di(k) = di(k)+di0;
        dj(k) = dj(k)+dj0;
        
    end
    
    %-- Weighted robust smoothing
    if kk==size(param.winsize,1), options.TolZ = 1e-3; end
    dj = smoothn({di,dj},sqrt(C),'robust',options);
    di = dj{1}; dj = dj{2};

end

if isfield(param,'ROI')
    [j,i] = meshgrid(1:N,1:M);
    ROI = interp2(j,i,ROI,jc,ic,'*nearest');
    di(~ROI) = NaN;
    dj(~ROI) = NaN;
end

% Return to previous warning states
warning(warn01.state,'MATLAB:smoothn:MaxIter'); 
warning(warn02.state,'MATLAB:smoothn:SLowerBound');
warning(warn03.state,'MATLAB:smoothn:SUpperBound');

end




function y = rmnan(x,order)

% Remove NaNs by inter/extrapolation
% see also INPAINTN
% Written by Louis Le Tarnec, RUBIC, 2012

sizx = size(x);
W = isfinite(x);

x(~W) = 0;
W = W(:); x = x(:);

missing_values = find(W==0)';

% Matrix defined by Buckley (equation 23)
% Biometrika (1994), 81, 2, pp. 247-58
d = length(sizx);
for i = 1:d
    n = sizx(i);
    e = ones(n,1);
    K = spdiags([e -2*e e],-1:1,n,n);
    K(1,1) = -1; K(n,n) = -1; %#ok
    M = 1;
    for j = 1:d
        if j==i, M = kron(K,M); end
        if j~=i
            m = sizx(j);
            I = spdiags(ones(1,m)',0,m,m);
            M = kron(I,M);
        end
    end
    if i==1, A = M; else, A = A+M; end
end
A = A^order;

% Linear system to be solved
x2 = -A*x;
x2 = x2(missing_values);
A = A(missing_values, missing_values);

% Solution
x2 = A\x2;
x(missing_values) = x2;
y = reshape(x,sizx);

end


function RunTheExample

% This example shows how to obtain the motion field of a rotating disk
% insonified with plane waves.

hp = helpdlg({['A rotating disk (diameter of 2 cm) was insonified by',...
    ' a series of 32 unsteered plane waves with a Verasonics scanner',...
    ' and a linear transducer at a PRF (pulse repetition frequency) of 10 kHz.'],...
    '',...
    'The RF signals were downsampled at 4/3 times (5 MHz) = 6.66 MHz.'},...
    'SPTRACK example');
waitfor(hp)

%%
% A rotating disk (diameter of 2 cm) was insonified by a series of 32
% unsteered plane waves with a Verasonics scanner, and a linear transducer,
% at a PRF (pulse repetition frequency) of 10 kHz. The RF signals were
% downsampled at 4/3 times (5 MHz) = 6.66 MHz. The properties of the
% linear array were:
% 
% 128 elements
% center frequency = 5 MHz
% pitch = 0.298 mm

%%
% Download the experimental RF data. The 3-D array RF contains 128 columns
% (as the transducer contained 128 elements), and its length is 32 in the
% third dimension (as 32 plane waves were transmitted).
load('PWI_disk.mat'); %#ok

%%
% Demodulate the RF signals with RF2IQ.
IQ = rf2iq(RF,param); %#ok

%%
% Create a 2.5-cm-by-2.5-cm image grid.
dx = 1e-4; % grid x-step (in m)
dz = 1e-4; % grid z-step (in m)
[x,z] = meshgrid(-1.25e-2:dx:1.25e-2,1e-2:dz:3.5e-2);

%%
% Create a Delay-And-Sum DAS matrix with DASMTX.
param.fnumber = []; % an f-number will be determined by DASMTX
M = dasmtx(1i*size(IQ),x,z,param,'nearest');

%%
% Beamform the I/Q signals.
%
% The 32 I/Q series can be beamformed simultaneously with the DAS matrix.
IQb = M*reshape(IQ,[],32);
IQb = reshape(IQb,[size(x) 32]);

%%
% Create the B-mode images with -30dB range.
I = bmode(IQb,30);

%%
% Create an ROI.
param.ROI = median(I,3)>64;

%%
% Track the speckles with SPTRACK.
param.winsize = [32 32; 24 24; 16 16]; % size of the subwindows
param.iminc = 4; % image increment
[Di,Dj,id,jd] = sptrack(I,param);

%%
% Display the motion field.
image(I(:,:,1))
colormap gray
hold on
h = quiver(jd,id,Dj,Di,3,'r');
set(h,'LineWidth',1)
hold off
title('Motion field (in pix) by speckle tracking')
axis equal off ij

end



