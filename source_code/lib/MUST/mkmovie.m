function [F,info,param] = mkmovie(varargin)

%MKMOVIE   Create movie frames and animated GIF for wave propagation
%   F = MKMOVIE(DELAYS,PARAM) simulates ultrasound RF radio-frequency
%   signals by using PFIELD and returns movie frames. The array elements
%   are excited at different time delays, given by DELAYS (in s). The
%   characteristics of the transmit and receive must be given in the
%   structure PARAM (see below for details).
%
%   Note: MKMOVIE works in a 2-D space.
%
%   >--- Try it: enter "mkmovie" in the command window for an example ---< 
%
%   By default, the ROI is of size 2L-by-2L, with L being the aperture
%   length of the array [i.e. L = pitch*(Number_of_elements-1)], and its
%   resolution is 50 pix/cm. These values can be modified through
%   PARAM.movie:
%       e.g. PARAM.movie = [ROI_width ROI_height resolution]
%       IMPORTANT: Pay attention to the units!
%                  Width and height are in cm, resolution is in pix/cm.
%
%   The output variable F is a 8-bit 3-D array that contains the frames
%   along the third dimension.
%
%   MKMOVIE uses PFIELD during transmit and receive. The parameters that
%   must be included in the structure PARAM are similar as those in PFIELD.
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%       *** MOVIE PROPERTIES ***
%   0)  PARAM.movie = [width height resolution duration fps]
%       IMPORTANT: Pay attention to the units!
%            ROI width and height are in cm, resolution is in pix/cm,
%            duration is in s, fps is frame per second.
%       The default is [2L 2L 50 15 10], with L (in cm!) being the aperture
%       length of the array [i.e. L = pitch*(Number_of_elements-1)]
%
%       *** TRANSDUCER PROPERTIES ***
%   1)  PARAM.fc: central frequency (in Hz, REQUIRED)
%   2)  PARAM.pitch: pitch of the linear array (in m, REQUIRED)
%   3)  PARAM.width: element width (in m, REQUIRED)
%       or PARAM.kerf: kerf width (in m, REQUIRED)
%       note: kerf = pitch-width 
%   4)  PARAM.radius: radius of curvature (in m)
%            The default is Inf (rectilinear array)
%   5)  PARAM.bandwidth: pulse-echo (2-way) 6dB fractional bandwidth (in %)
%            The default is 75%.
%   6)  PARAM.baffle: property of the baffle:
%            'soft' (default), 'rigid' or a scalar >= 0.
%            See "Note on BAFFLE property" in PFIELD for details
%
%       *** MEDIUM PARAMETERS ***
%   7)  PARAM.c: longitudinal velocity (in m/s, default = 1540 m/s)
%   8)  PARAM.attenuation: attenuation coefficient (dB/cm/MHz, default: 0)
%            Notes: A linear frequency-dependence is assumed.
%                   A typical value for soft tissues is ~0.5 dB/cm/MHz.
%
%       *** TRANSMIT PARAMETERS ***
%   9)  PARAM.TXapodization: transmision apodization (default: no apodization)
%   10) PARAM.TXnow: pulse length in number of wavelengths (default: 1)
%            Use PARAM.TXnow = Inf for a mono-harmonic signal.
%   11) PARAM.TXfreqsweep: frequency sweep for a linear chirp (default: [])
%            To be used to simulate a linear TX chirp.
%            See "Note on CHIRP signals" in PFIELD for details
%
%   Other syntaxes:
%   --------------
%   F = MKMOVIE(X,Z,RC,DELAYS,PARAM) also simulates backscattered echoes.
%   The scatterers are characterized by their coordinates (X,Z) and
%   reflection coefficients RC. X, Z and RC must be of same size.
%
%   [F,INFO] = MKMOVIE(...) returns image information in the structure
%   INFO. INFO.Xgrid and INFO.Zgrid are the x- and z-coordinates of the
%   image. INFO.TimeStep is the time step between two consecutive frames.
%
%   [F,INFO,PARAM] = MKMOVIE(...) updates the fields of the PARAM
%   structure.
%
%   [...] = MKMOVIE without any input argument provides an interactive
%   example designed to produce a movie using a 2.5 MHz phased-array
%   transducer.
%
%   [...] = MKMOVIE(...,OPTIONS) uses the structure OPTIONS to adjust the
%   simulations performed by PFIELD:
%   
%   OPTIONS:
%   -------
%      %-- FREQUENCY SAMPLES --%
%   1) Only frequency components of the transmitted signal in the range
%      [0,2fc] with significant amplitude are considered. The default
%      relative amplitude is -60 dB in MKMOVIE. You can change this value
%      by using the following:
%          [...] = MKMOVIE(...,OPTIONS),
%      where OPTIONS.dBThresh is the threshold in dB (default = -60).
%   ---
%      %-- FULL-FREQUENCY DIRECTIVITY --%   
%   2) By default, the directivity of the elements depends only on the
%      center frequency. This makes the calculation faster. To make the
%      directivities fully frequency-dependent, use: 
%          [...] = MKMOVIE(...,OPTIONS),
%      with OPTIONS.FullFrequencyDirectivity = true (default = false).
%   ---
%       %-- ELEMENT SPLITTING --%   
%   3)  Each transducer element of the array is split into small segments.
%       The length of these small segments must be small enough to ensure
%       that the far-field model is accurate. By default, the elements are
%       split into M segments, with M being defined by:
%           M = ceil(element_width/smallest_wavelength);
%       To modify the number M of subelements by splitting, you may adjust
%       OPTIONS.ElementSplitting. For example, OPTIONS.ElementSplitting = 1
%   ---
%       %-- WAIT BAR --%   
%   4)  If OPTIONS.WaitBar is true, a wait bar appears (only if the number
%       of frequency samples >10). Default is true.
%   ---
%   
%   CREATE an animated GIF:
%   ----------------------
%   [...] = MKMOVIE(...,FILENAME) creates a 10-fps animated GIF to the file
%   specified by FILENAME. The duration of the animated GIF is ~15 seconds.
%   You can modify the duration and fps by using PARAM.movie (see above).
%
%   Example:
%   -------
%   %-- Generate a diverging wave using a phased-array transducer
%   % Phased-array @ 2.7 MHz:
%   param = getparam('P4-2v');
%   % TX time delays for a 90-degree wide diverging wave:
%   dels = txdelay(param,0,pi/2);
%   % Scatterers' position:
%   n = 20;
%   x = rand(n,1)*8e-2-4e-2;
%   z = rand(n,1)*10e-2;
%   % Backscattering coefficient
%   RC = (rand(n,1)+1)/2;
%   % Image size (in cm)
%   param.movie = [8 10];
%   % Movie frames
%   [F,info] = mkmovie(x,z,RC,dels,param);
%   % Check the movie frames
%   figure
%   colormap([1-hot(128); hot(128)]);
%   for k = 1:size(F,3)
%       image(info.Xgrid,info.Zgrid,F(:,:,k))
%       hold on
%       scatter(x,z,5,'w','filled')
%       hold off
%       axis equal off
%       title([int2str(info.TimeStep*k*1e6) ' \mus'])
%       drawnow
%   end 
%   
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also PFIELD, SIMUS, TXDELAY.
%
%   -- Damien Garcia -- 2017/10, last update 2021/01/02
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

if nargin==0
    if nargout>0
        [F,info,param] = RunTheExample;
    else
        RunTheExample;
    end
    return
end

narginchk(2,7)
if ischar(varargin{nargin})
    gifname = varargin{nargin};
    isGIF = true;
    Nargin = nargin-1;
else
    isGIF = false;
    Nargin = nargin;
end

switch Nargin
    case 2 % mkmovie(delaysTX,param)
        delaysTX = varargin{1};
        param = varargin{2};
        x = []; z = []; RC = [];
        options = [];
    case 3 % mkmovie(delaysTX,param,options)
        delaysTX = varargin{1};
        param = varargin{2};
        options = varargin{3};
        x = []; z = []; RC = []; 
    case 5 % mkmovie(x,z,RC,delaysTX,param)
        x = varargin{1};
        z = varargin{2};
        RC = varargin{3};
        delaysTX = varargin{4};
        param = varargin{5};
        options = [];
    case 6 % mkmovie(x,z,RC,delaysTX,param,options)
        x = varargin{1};
        z = varargin{2};
        RC = varargin{3};
        delaysTX = varargin{4};
        param = varargin{5};
        options = varargin{6};
    otherwise
        error('Wrong input arguments')
end

param = IgnoreCaseInFieldNames(param);
options = IgnoreCaseInFieldNames(options);
options.CallFun = 'mkmovie';


%------------------------%
% CHECK THE INPUT SYNTAX %
%------------------------%

assert(isstruct(param),'The structure PARAM is required.')
assert(isequal(size(z),size(x),size(RC)),'X, Z and RC must be of same size.')

%-- Wait bar
if ~isfield(options,'WaitBar')
    options.WaitBar = true;
end
assert(isscalar(options.WaitBar) && islogical(options.WaitBar),...
    'OPTIONS.WaitBar must be a logical scalar (true or false).')

%-- Check if syntax errors may appear when using PFIELD
try
    opt = options;
    opt.WaitBar = false;
    [~,param] = pfield([],[],delaysTX,param,opt);
catch ME
    throw(ME)
end

%-- Movie properties
NumberOfElements = param.Nelements; % number of array elements
L = param.pitch*(NumberOfElements-1);
if ~isfield(param,'movie')
    % param.movie = [width height pix/cm duration fps]
    % NOTE: width and height are in cm
    param.movie = [200*L 200*L 50 15 10]; % default
else
    assert(isnumeric(param.movie) & numel(param.movie)>1 &...
        numel(param.movie)<6 & all(param.movie>0),...
        'PARAM.movie must contain two to five positive parameters.')
end

% default resolution = 50 pix/cm, default duration = 15 s, default fps = 10
paramMOVdefault = [NaN NaN 50 15 10];
n = numel(param.movie);
param.movie = [param.movie(:)' paramMOVdefault(n+1:5)];


%-- dB threshold (i.e. faster computation if lower)
if ~isfield(options,'dBThresh')
    options.dBThresh = -60; % default is -60dB in MKMOVIE
end
assert(isscalar(options.dBThresh) && isnumeric(options.dBThresh) &&...
    options.dBThresh<0,'OPTIONS.dBThresh must be a negative scalar.')

%-- Frequency step (scaling factor)
% The frequency step is determined automatically. It is tuned to avoid
% aliasing in the temporal domain. The frequency step can be adjusted by
% using a scaling factor. For a smoother result, you may use a scaling
% factor<1.
if ~isfield(options,'FrequencyStep')
    options.FrequencyStep = 1;
end
assert(isscalar(options.FrequencyStep) &&...
    isnumeric(options.FrequencyStep) && options.FrequencyStep>0,...
    'OPTIONS.FrequencyStep must be a positive scalar.')
if options.FrequencyStep>1
    warning('MUST:FrequencyStep',...
            'OPTIONS.FrequencyStep is >1: aliasing may be present!')
end


%-------------------------------%
% end of CHECK THE INPUT SYNTAX %
%-------------------------------%


%-- Image grid
ROIwidth = param.movie(1)*1e-2; % image width (in m)
ROIheight = param.movie(2)*1e-2; % image height (in m)
pixsize = 1/param.movie(3)*1e-2; % pixel size (in m)
xi = (pixsize/2:pixsize:ROIwidth-pixsize/2);
zi = (pixsize/2:pixsize:ROIheight-pixsize/2);
[xi,zi] = meshgrid(xi-mean(xi),zi);


%-- Frequency sampling
maxD = hypot((ROIwidth+L)/2,ROIheight); % maximum travel distance
df = 1/2/(maxD/param.c);
df = df*options.FrequencyStep;
Nf = 2*ceil(param.fc/df)+1; % number of frequency samples


%-- Run PFIELD to calculate the RF spectra
SPECT = zeros(Nf,numel(xi),'single'); % will contain the RF spectra
options.FrequencyStep = df;
options.ElementSplitting = 1;
options.RC = RC;
options.x = x;
options.z = z;
%-
% we need idx...
opt = options;
opt.x = []; opt.z = []; opt.RC = [];
opt.WaitBar = false;
[~,~,~,idx] = pfield([],[],delaysTX,param,opt);
%-
[~,~,SPECT(idx,:)] = pfield(xi,zi,delaysTX,param,options);

%-- IFFT to recover the time-resolved signals
if options.WaitBar
   wbname = 'MKMOVIE / www.biomecardio.com';
   h = waitbar(0,'MKMOVIE creates the frames...','Name',wbname);
end
%
F = SPECT; clear SPECT
F = reshape(F,Nf,size(xi,1),size(xi,2));
F = shiftdim(F,1);
%
if options.WaitBar, waitbar(1/4,h); end
%
F = cat(3,F,conj(flip(F(:,:,2:end-1),3)));
%
if options.WaitBar, waitbar(1/2,h); end
%
F = ifft(F,[],3,'symmetric');
%
if options.WaitBar, waitbar(3/4,h); end
%
F = flip(F,3);
F = F(:,:,1:round(size(F,3)/2));
F = F/max(abs(F),[],'all');
%
if options.WaitBar, waitbar(1,h); end
%
% F = abs(F).^0.6.*sign(F);
F = uint8((F+1)/2*255);
%
if options.WaitBar, close(h); end

%-- some information about the movie
if nargout>1
    info.Xgrid = xi(1,:); % in m
    info.Zgrid = zi(:,1)'; % in m
    info.TimeStep = maxD/param.c/size(F,3); % in s
end

%-- animated GIF
if isGIF
    %- the matrix f contains the scatterers
    f = zeros([size(F,1) size(F,2)]);
    if ~isempty(x)
        maxRC = max(RC(:));
        for k = 1:length(x(:))
            [~,i] = min(abs(zi(:,1)-z(k)));
            [~,j] = min(abs(xi(1,:)-x(k)));
            f(i,j) = RC(k)/maxRC;
        end
        n = 2*round(param.movie(3)/10)+1;
        f = conv2(f,blackman(n)*blackman(n)','same');
    end
    f = uint8(f*128);
    
    %- add the signature
    % Please do not remove it
    if size(f,1)>37 && size(f,2)>147
        f(end-36:end-5,6:147) = Signature/2;
    else
        f = 0;
    end
    
    map = flipud([1-hot(128); hot(128)]);
        
    %- create the GIF movie
    Tmov = param.movie(4); % movie duration in s
    fps = param.movie(5); % frame per second
    
    nk = round(size(F,3)/(Tmov*fps)); % increment
        
    for k = 1:nk:size(F,3)
        if k==1
            imwrite(F(:,:,k)+f,map,gifname,'gif','LoopCount',Inf,...
                'DelayTime',1/fps);
        else
            imwrite(F(:,:,k)+f,map,gifname,'gif','WriteMode','append',...
                'DelayTime',1/fps);
        end
    end
    
end




end


function structArray = IgnoreCaseInFieldNames(structArray)

switch inputname(1)
    case 'param'
        fieldLIST = {'attenuation','baffle','bandwidth','c','fc',...
            'fnumber','focus','fs','height','kerf','movie','Nelements',...
            'pitch','radius','RXangle','RXdelay','TXapodization',...
            'TXfreqsweep','TXnow','t0','width'};
    case 'options'
        if isstruct(structArray)
            fieldLIST = {'dBThresh','ElementSplitting',...
                'FullFrequencyDirectivity','FrequencyStep','ParPool',...
                'WaitBar'};
        else
            return
        end
end

OldFieldNames = fieldnames(structArray);
tmp = lower(OldFieldNames);
assert(length(tmp)==length(unique(tmp)),...
    ['The structure ' upper(inputname(1)),...
    ' contains duplicate field names (when ignoring case).'])

[idx,loc] = ismember(lower(fieldLIST),tmp);
idx = find(idx); loc = loc(idx);
for k = 1:length(idx)
    tmp = eval(['structArray.' OldFieldNames{loc(k)}]); %#ok
    structArray = rmfield(structArray,OldFieldNames{loc(k)});
    eval(['structArray.' fieldLIST{idx(k)} ' = tmp;'])
end

end


function [F,info,param] = RunTheExample

% 2.7 MHz phased-array
param = getparam('P4-2v');

% Location of the scatterers
figure('Name','MUST - Matlab Ultrasound Toolbox - Damien Garcia',...
    'NumberTitle','off');
plot((param.Nelements-1)*param.pitch/2*[-1 1],[0 0],'g','Linewidth',5)
hold on
axis equal ij
axis([-4e-2 4e-2 0 10e-2])
xlabel('x-position (m)')
ylabel('z-position (m)')
hp = helpdlg({'We have a 2.5 MHz phased array.',...
    'Choose five points (by clicking) in the figure.'},...
    'Scatterer selection');
waitfor(hp)
x = zeros(1,5); z = x;
for k = 5:-1:1
    if k==5
        title('Choose 5 points')
    elseif k==1
        title('A last one')
    else
        title(['Still ' int2str(k) ' points'])
    end
    [x(k),z(k)] = ginput(1);
    scatter(x(k),z(k),8,'r','filled')
end
hold off
RC = ones(1,5);

% TX time delays (diverging wave):
dels = txdelay(param,0,pi/2);

% Image size (in cm)
param.movie = [8 10];

% Movie frames
title('Calculating...')
[F,info,param] = mkmovie(x,z,RC,dels,param);

% Check the movie frames
colormap([1-hot(128); hot(128)]);
for k = 1:size(F,3)
    image(info.Xgrid,info.Zgrid,F(:,:,k))
    hold on
    plot((param.Nelements-1)*param.pitch/2*[-1 1],[0 0],'g','Linewidth',5)
    scatter(x,z,8,'r','filled')
    hold off
    axis equal ij
    axis([-4e-2 4e-2 0 10e-2])
    xlabel('x-position (m)')
    ylabel('z-position (m)')
    title([int2str(info.TimeStep*k*1e6) ' \mus'])
    drawnow
end

end


function signature = Signature

signature = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,255,255,115,0,0,0,0,0,0,50,244,255,157,0,0,126,251,77,0,0,0,0,50,244,157,0,0,0,126,255,255,255,222,17,11,233,255,255,255,255,255,255,255,247,31,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,251,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,255,255,222,17,0,0,0,0,0,166,255,255,157,0,0,126,251,77,0,0,0,0,50,244,157,0,0,126,251,77,0,0,126,157,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,251,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,251,148,251,77,0,0,0,0,11,233,155,244,157,0,0,126,251,77,0,0,0,0,50,244,157,0,11,233,195,0,0,0,0,0,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,80,255,255,255,195,0,0,126,247,81,244,255,157,0,0,80,255,255,255,195,0,0,0,0,166,255,255,255,157,0,6,200,255,255,255,255,195,0,0,80,255,255,255,195,0,0,0,0,80,255,255,251,89,233,157,0,0,0,0,0,0,126,247,42,233,255,255,157,0,0,166,247,31,0,0,0,80,251,77,0,0,0,0,80,251,89,233,195,0,0,0,0,126,247,81,244,157,0,0,126,251,77,0,0,0,0,50,244,157,0,11,233,222,17,0,0,0,0,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,80,255,115,0,0,166,115,0,126,255,251,77,0,0,0,80,255,115,0,50,244,157,0,0,166,157,0,0,80,255,115,0,0,80,251,77,0,0,0,80,255,115,0,50,244,157,0,0,80,255,115,0,50,244,255,157,0,0,0,0,0,0,126,255,251,77,0,80,255,157,0,80,255,115,0,0,0,166,222,17,0,0,0,0,80,251,77,166,247,31,0,0,6,200,195,50,244,157,0,0,126,251,77,0,0,0,0,50,244,157,0,0,166,255,195,0,0,0,0,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,6,200,195,0,0,0,0,0,0,126,255,157,0,0,0,6,200,195,0,0,0,126,247,31,0,0,0,0,0,11,233,195,0,0,80,251,77,0,0,6,200,195,0,0,0,126,247,31,6,200,195,0,0,0,50,244,157,0,0,0,0,0,0,126,255,115,0,0,0,166,222,17,6,200,195,0,0,11,233,157,0,0,0,0,0,80,251,77,50,244,115,0,0,80,251,77,50,244,157,0,0,126,251,77,0,0,0,0,50,244,157,0,0,0,166,255,255,222,17,0,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,50,244,157,0,0,0,0,0,0,126,251,77,0,0,0,50,244,115,0,0,0,80,251,77,0,0,0,0,0,6,200,195,0,0,80,251,77,0,0,50,244,115,0,0,0,80,251,77,11,233,157,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,126,251,77,0,126,251,77,0,80,251,77,0,0,0,0,0,80,251,77,6,200,222,17,0,166,222,17,50,244,157,0,0,126,251,77,0,0,0,0,50,244,157,0,0,0,0,6,200,255,255,157,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,50,244,115,0,0,0,0,0,0,126,251,77,0,0,0,50,244,255,255,255,255,255,251,77,0,6,200,255,255,255,255,195,0,0,80,251,77,0,0,50,244,255,255,255,255,255,251,77,50,244,157,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,126,251,77,0,50,244,157,0,166,222,17,0,0,0,0,0,80,251,77,0,126,251,77,11,233,157,0,50,244,157,0,0,126,251,77,0,0,0,0,50,244,157,0,0,0,0,0,0,11,233,251,77,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,50,244,157,0,0,0,0,0,0,126,251,77,0,0,0,50,244,115,0,0,0,0,0,0,6,200,222,17,0,6,200,195,0,0,80,251,77,0,0,50,244,115,0,0,0,0,0,0,11,233,157,0,0,0,11,233,157,0,0,0,0,0,0,126,247,31,0,0,0,126,247,31,0,6,200,222,29,233,157,0,0,0,0,0,0,80,251,77,0,11,233,157,126,247,31,0,50,244,157,0,0,80,251,77,0,0,0,0,50,244,115,0,0,0,0,0,0,0,80,255,115,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,11,233,195,0,0,0,0,0,0,126,251,77,0,0,0,11,233,157,0,0,0,0,0,0,50,244,157,0,0,6,200,195,0,0,80,251,77,0,0,11,233,157,0,0,0,0,0,0,6,200,195,0,0,0,50,244,157,0,0,0,0,0,0,126,251,77,0,0,0,166,222,17,0,0,126,251,148,251,77,0,0,0,0,0,0,80,251,77,0,0,166,251,223,195,0,0,50,244,157,0,0,50,244,157,0,0,0,0,126,251,77,0,0,0,0,0,0,0,126,251,77,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,126,255,115,0,0,166,157,0,126,251,77,0,0,0,0,126,255,115,0,0,80,222,17,11,233,195,0,0,166,255,195,0,0,50,244,157,0,0,0,126,255,115,0,0,80,222,17,0,126,255,115,0,50,244,255,157,0,0,0,0,0,0,126,255,251,77,0,80,255,115,0,0,0,11,233,255,222,17,0,0,0,0,0,0,80,251,77,0,0,80,255,251,77,0,0,50,244,157,0,0,0,166,255,115,0,0,50,244,157,0,0,80,247,31,0,0,50,244,195,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,126,255,255,255,195,0,0,126,251,77,0,0,0,0,0,80,255,255,255,222,17,0,0,50,244,255,255,157,166,195,0,0,0,126,255,255,195,0,0,80,255,255,255,222,17,0,0,0,126,255,255,247,37,200,157,0,0,0,0,0,0,126,222,67,244,255,255,115,0,0,0,0,0,166,255,115,0,0,0,0,0,0,0,80,251,77,0,0,6,200,222,17,0,0,50,244,157,0,0,0,0,80,255,255,255,255,115,0,0,0,0,80,255,255,255,255,157,0,0,0,0,0,80,255,115,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,166,247,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,233,195,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,126,251,77,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,200,115,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,222,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,200,115,0,0,0,0,0,0,166,222,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,222,17,0,166,222,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,200,115,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,80,222,17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
126,195,0,0,80,247,31,0,11,233,115,126,195,0,0,80,247,31,0,11,233,115,126,195,0,0,80,247,31,0,11,233,115,0,0,0,0,6,200,126,233,255,251,77,0,0,166,195,0,0,6,200,255,255,247,31,0,6,197,89,233,255,247,37,200,255,251,77,0,6,200,255,255,247,31,0,0,50,244,255,255,115,0,50,244,255,255,115,0,11,229,148,255,247,31,50,244,255,251,148,222,17,0,166,195,0,6,200,255,255,247,31,0,0,0,0,0,50,244,255,255,115,0,6,200,255,255,247,31,0,6,197,89,233,255,247,37,200,255,251,77,0
50,236,31,0,122,223,115,0,50,236,31,50,236,31,0,122,223,115,0,50,236,31,50,236,31,0,122,223,115,0,50,236,31,0,0,0,0,6,200,251,77,0,126,247,31,0,166,195,0,6,200,195,0,0,126,247,31,6,200,247,31,6,200,255,115,0,126,222,17,166,195,0,0,126,222,17,11,233,157,0,11,200,29,200,17,0,126,222,17,11,233,247,31,0,6,200,157,0,11,233,222,17,0,166,195,6,200,195,0,0,126,247,31,0,0,0,11,233,157,0,11,200,24,200,195,0,0,126,247,31,6,200,247,31,6,200,255,115,0,126,222,17
11,233,115,6,192,111,157,0,126,195,0,11,233,115,6,192,111,157,0,126,195,0,11,233,115,6,192,111,157,0,126,195,0,0,0,0,0,6,200,115,0,0,11,229,77,0,166,195,0,50,236,31,0,0,11,229,77,6,200,115,0,0,126,195,0,0,80,247,81,240,77,0,0,11,229,77,80,247,31,0,0,0,0,0,0,0,50,236,31,11,233,115,0,0,50,236,31,0,0,80,222,17,0,166,195,50,236,31,0,0,11,229,77,0,0,0,80,247,31,0,0,0,50,236,31,0,0,11,229,77,6,200,115,0,0,126,195,0,0,80,247,31
0,166,157,11,173,50,211,24,200,115,0,0,166,157,11,173,50,211,24,200,115,0,0,166,157,11,173,50,211,24,200,115,0,0,0,0,0,6,200,115,0,0,11,233,115,0,166,195,0,80,247,31,0,0,6,200,115,6,200,115,0,0,126,195,0,0,80,247,81,244,255,255,255,255,251,77,80,222,17,0,0,0,0,126,255,255,255,247,31,11,233,115,0,0,80,247,31,0,0,80,222,17,0,166,195,80,247,31,0,0,6,200,115,0,0,0,80,222,17,0,0,0,80,247,31,0,0,6,200,115,6,200,115,0,0,126,195,0,0,80,247,31
0,80,222,98,157,6,197,89,225,31,0,0,80,222,98,157,6,197,89,225,31,0,0,80,222,98,157,6,197,89,225,31,0,0,0,0,0,6,200,115,0,0,11,229,77,0,166,195,0,50,236,31,0,0,11,229,77,6,200,115,0,0,126,195,0,0,80,247,81,236,31,0,0,0,0,0,80,247,31,0,0,0,80,247,31,0,50,236,31,11,233,115,0,0,50,236,31,0,0,80,222,17,0,166,195,50,236,31,0,0,11,229,77,0,0,0,80,247,31,0,0,0,50,236,31,0,0,11,229,77,6,200,115,0,0,126,195,0,0,80,247,31
0,11,229,186,115,0,166,182,195,0,0,0,11,229,186,115,0,166,182,195,0,0,0,11,229,186,115,0,166,182,195,0,0,6,200,195,0,6,200,251,77,0,126,222,17,0,166,195,0,6,200,195,0,0,126,222,17,6,200,115,0,0,126,195,0,0,80,247,37,200,157,0,0,0,0,0,11,233,157,0,11,200,98,247,31,0,166,247,31,11,233,115,0,0,11,233,157,0,50,244,222,17,0,166,195,6,200,195,0,0,126,222,17,6,200,195,11,233,157,0,11,200,24,200,195,0,0,126,222,17,6,200,115,0,0,126,195,0,0,80,247,31
0,0,166,247,31,0,80,255,115,0,0,0,0,166,247,31,0,80,255,115,0,0,0,0,166,247,31,0,80,255,115,0,0,6,200,195,0,6,197,119,244,255,247,31,0,0,166,195,0,0,11,233,255,255,222,17,0,6,200,115,0,0,126,195,0,0,80,247,31,6,200,255,255,255,222,17,0,50,244,255,251,77,0,166,255,255,167,225,31,11,233,115,0,0,0,50,244,255,247,81,211,17,0,166,195,0,11,233,255,255,222,17,0,6,200,195,0,50,244,255,251,77,0,11,233,255,255,222,17,0,6,200,115,0,0,126,195,0,0,80,247,31
];

end


