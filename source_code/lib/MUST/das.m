function [bfSIG,param] = das(SIG,x,z,varargin)

%DAS   Delay-And-Sum beamforming of RF or I/Q signals
%   BFSIG = DAS(SIG,X,Z,DELAYS,PARAM) beamforms the RF or I/Q signals
%   stored in the array SIG, and returns the beamformed signals BFSIG. The
%   signals are beamformed at the points specified by X and Z.
%
%   Try it: enter "das" in the command window for an example.
%
%     1) SIG must be a 2-D or 3-D array. The first dimension (i.e. each
%        column) corresponds to a single RF or I/Q signal over (fast-)
%        time, with the FIRST COLUMN corresponding to the FIRST ELEMENT.
%        Several 2-D signals can be stacked along the 3rd dimension.
%     2) DELAYS are the transmit time delays (in s). One must have
%        numel(DELAYS) = size(SIG,2). If a sub-aperture was used during
%        transmission, use DELAYS(i) = NaN if element #i of the array was
%        off.
%     3) PARAM is a structure that contains the parameter values required
%        for the delay-and-sum (see below for details).
%
%   Note: SIG must be complex when DAS beamforming I/Q data
%         (i.e. SIG = complex(I,Q) = I + 1i*Q).
%
%   DAS(SIG,X,Z,PARAM) uses DELAYS = PARAM.TXdelay.
%
%   DAS(...,METHOD) specifies the interpolation method. The available
%   methods are decribed in NOTE #3 below.
%
%   [BFSIG,PARAM] = DAS(...) also returns the structure PARAM.
%
%   ---
%   NOTE #1: X- and Z-axes
%   The migrated signals are calculated at the points specified by (X,Z).
%   Conventional axes are used:
%   i) For a LINEAR array, the X-axis is PARALLEL to the transducer and
%      points from the first (leftmost) element to the last (rightmost)
%      element (X = 0 at the CENTER of the transducer). The Z-axis is
%      PERPENDICULAR to the transducer and points downward (Z = 0 at the
%      level of the transducer, Z increases as depth increases).
%   ii) For a CONVEX array, the X-axis is parallel to the chord and Z = 0
%       at the level of the chord.
%   ---
%   NOTE #2: DAS uses a standard diffraction summation (delay-and-sum). It
%   calls the function DASMTX.
%   ---
%   NOTE #3: Interpolation method
%   By default DAS uses a linear interpolation to generate the DAS matrix.
%   To specify the interpolation method, use DAS(...,METHOD), with METHOD
%   being:
%      'nearest'   - nearest neighbor interpolation
%      'linear'    - (default) linear interpolation
%      'quadratic' - quadratic interpolation
%      'lanczos3'  - 3-lobe Lanczos interpolation
%      '5points'   - 5-point least-squares parabolic interpolation
%      'lanczos5'  - 5-lobe Lanczos interpolation
%
%   The linear interpolation (it is a 2-point method) returns a matrix
%   twice denser than the nearest-neighbor interpolation. It is 3, 4, 5, 6
%   times denser for 'quadratic', 'lanczos3', '5points', 'lanczos5',
%   respectively (they are 3-to-6-point methods).
%   ---
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%   1) PARAM.fs: sample frequency (in Hz, REQUIRED)
%   2) PARAM.pitch: pitch of the transducer (in m, REQUIRED)
%   3) PARAM.fc: center frequency (in Hz, REQUIRED for I/Q signals)
%   4) PARAM.radius: radius of curvature (in m, default = Inf, linear array)
%   5) PARAM.TXdelay: transmission delays (in s, required if DELAYS is not given)
%   6) PARAM.c: longitudinal velocity (in m/s, default = 1540 m/s)
%   7) PARAM.t0: start time for reception (in s, default = 0 s)
%
%   A note on the f-number
%   ----------------------
%   The f-number is defined by the ratio (depth)/(aperture size). A null
%   f-number (PARAM.fnumber = 0) means that the full aperture is used
%   during DAS-beamforming. This might be a suboptimal strategy since the
%   array elements have some directivity.
%
%   Use PARAM.fnumber = [] to obtain an "optimal" f-number, which is
%   estimated from the element directivity (and depends on fc, bandwidth,
%   element width):
%
%   8)  PARAM.fnumber: reception f-number (default = 0, i.e. full aperture)
%   9)  PARAM.width: element width (in m, REQUIRED if PARAM.fnumber = [])
%        or PARAM.kerf: kerf width (in m, REQUIRED if PARAM.fnumber = [])
%        note: width = pitch-kerf 
%   10) PARAM.bandwidth: pulse-echo 6dB fractional bandwidth (in %)
%            The default is 75% (used only if PARAM.fnumber = []).
%
%   Advanced option for vector Doppler (Reception angle):
%   ---------------------------------------------------
%   11) PARAM.RXangle: reception angles (in rad, default = 0)   
%       This option can be used for vector Doppler. Beamforming with at
%       least two (sufficiently different) reception angles enables
%       different Doppler directions and, in turn, vector Doppler.
%       (This option is not yet available for convex arrays)
%
%   Passive imaging
%   ---------------
%   12) PARAM.passive: must be true for passive imaging (i.e. no transmit).
%       The default is false.
%
%   If you need to beamform a large series of ultrasound signals acquired
%   with a same probe and a same transmit sequence, DASMTX is recommended.
%
%
%   REFERENCES:
%   ----------
%   1) If you use DAS or DASMTX, please cite:
%      V Perrot, M Polichetti, F Varray, D Garcia. So you think you can
%      DAS? A viewpoint on delay-and-sum beamforming. Ultrasonics 111,
%      106309. <a
%      href="matlab:web('https://www.biomecardio.com/publis/ultrasonics21.pdf')">PDF here</a>
%   2) If you use PARAM.RXangle for vector Doppler, please also cite:
%      Madiena C, Faurie J, Porée J, Garcia D. Color and vector flow
%      imaging in parallel ultrasound with sub-Nyquist sampling. IEEE Trans
%      Ultrason Ferroelectr Freq Control, 2018;65:795-802. <a
%      href="matlab:web('https://www.biomecardio.com/publis/ieeeuffc18a.pdf')">PDF here</a>
%
%
%   Example:
%   -------
%   %-- Generate RF signals using a phased-array transducer
%   % Phased-array @ 2.5 MHz:
%   param = getparam('PA4-2/20');
%   % TX time delays (80-degree-wide diverging wave)
%   dels = txdelay(param,0,80/180*pi);
%   % Scatterers' position:
%   xs = [(-1:0.5:1)*4e-2 zeros(1,5)];
%   zs = [ones(1,5)*6e-2 (2:2:10)*1e-2];
%   % Backscattering coefficient
%   BSC = [ones(1,9) 0];
%   % RF signals:
%   param.fs = 4*param.fc; % sampling frequency
%   RF = simus(xs,zs,BSC,dels,param);
%   % Plot the RF signals
%   subplot(121)
%   plot((RF(:,1:7:64)/max(RF(:))+(1:10)*2)',...
%      (0:size(RF,1)-1)/param.fs*1e6,'k')
%   set(gca,'XTick',(1:10)*2,'XTickLabel',int2str((1:7:64)'))
%   title('RF signals')
%   xlabel('Element number'), ylabel('time (\mus)')
%   xlim([0 22]), axis ij
%
%   %-- Demodulation and beamforming
%   % Demodulation
%   IQ = rf2iq(RF,param);
%   % Beamforming grid
%   [th,r] = meshgrid(linspace(-40,40,128)/180*pi+pi/2,...
%      linspace(1,9,256)*1e-2);
%   [x,z] = pol2cart(th,r);
%   % Beamformed IQ
%   IQb = das(IQ,x,z,dels,param);
%   % Beamformed image
%   subplot(122)
%   pcolor(x*100,z*100,abs(IQb).^.5)
%   colormap(gray)
%   title('Gamma-compressed image')
%   xlabel('[cm]'), ylabel('[cm]')
%   shading interp, axis equal ij tight
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also DASMTX, TXDELAY, RF2IQ, GETPARAM, SIMUS.
%
%   -- Damien Garcia -- 2012/05, last update 2020/02
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>


if nargin==0
    if nargout>0
        [bfSIG,param] = RunTheExample;
    else
        RunTheExample;
    end
    return
end

assert(nargin>3,'Not enough input arguments.')
assert(nargin<7,'Too many input arguments.')

siz0 = size(x);
assert(isequal(siz0,size(z)),'X and Z must of same size.')
assert(ndims(SIG)<=3,['SIG must be a 2D or 3D array whose each column ',...
    'corresponds to an RF or IQ signal acquired by a single element']);

[nl,nc,~] = size(SIG);

%-- check if we have I/Q signals
isIQ = ~isreal(SIG);


%-- DAS matrix using DASMTX
% Matlab stores sparse matrices in compressed sparse column format. When
% called by DAS, DASMTX generates a transpose DAS matrix if it requires
% less memory.
if isstruct(varargin{end})
    varargin{end}.TransposeDASMatrix = [];
elseif nargin>4 && isstruct(varargin{end-1})
    varargin{end-1}.TransposeDASMatrix = [];
else
    error('A PARAM structure is missing!')
end
%
try
    [M,param] = dasmtx((~isIQ*1+isIQ*1i)*[nl nc],x,z,varargin{:});
catch ME
    throw(ME)
end



%-- Delay-and-Sum
SIG = reshape(SIG,nl*nc,[]);


if ~isfield(param,'isADAPTIVE') || ~param.isADAPTIVE
    
    if param.TransposeDASMatrix
        bfSIG = (SIG.'*M).';
    else
        bfSIG = M*SIG;
    end
    bfSIG = reshape(bfSIG,[siz0 size(SIG,2)]);
    
    
    
    
elseif param.isADAPTIVE % (unpublished) adaptive DAS (Beta version!)
    
    assert(isIQ && size(SIG,2)==1)
    if param.TransposeDASMatrix, M = M.'; end
    
    vec = @(x) x(:);
    SIG = sparse((1:nl*nc)',vec(repmat(1:nc,nl,1)),SIG,nl*nc,nc);
    Hyperb = full(M*SIG); % the rows contain the IQ diffraction hyperbolas
    A = angle(Hyperb);
    A(A==0) = NaN;
    A = unwrap(A,[],2)';

    N = nc; L = param.pitch*(N-1);
    idx = round(x(:)*(N-1)/L+(N+1)/2); % locations of the hyperbola vertices
    idx(idx<1) = 1; idx(idx>N) = N;
    Ai = A(idx+(0:numel(x)-1)'*nc); % phases at the hyperbola vertices
    
    Wpha = exp(1i*(A'-Ai)); % phase corrections
    Wamp = max(0,1-abs(A'-Ai)/pi); % phase-dispersion-based weights

    bfSIG = sum(Hyperb.*Wamp.*Wpha,2,'omitnan');
    bfSIG = reshape(bfSIG,siz0);
    
end


param = rmfield(param,'TransposeDASMatrix');

end


function [IQb,param] = RunTheExample

%-- Generate RF signals using a phased-array transducer

% Phased-array @ 2.7 MHz:
param = getparam('P4-2v');

% TX time delays (80-degree-wide diverging wave)
dels = txdelay(param,0,80/180*pi);

% Scatterers' position:
xs = [(-1:0.5:1)*4e-2 zeros(1,5)];
zs = [ones(1,5)*6e-2 (2:2:10)*1e-2];

% Backscattering coefficient
BSC = [ones(1,9) 0];

% RF signals:
param.fs = 4*param.fc; % sampling frequency
RF = simus(xs,zs,BSC,dels,param);

% Plot the RF signals
figure
subplot(121)
plot((RF(:,1:7:64)/max(RF(:))+(1:10)*2)',...
    (0:size(RF,1)-1)/param.fs*1e6,'k')
set(gca,'XTick',(1:10)*2,'XTickLabel',int2str((1:7:64)'))
title('RF signals')
xlabel('Element number'), ylabel('time (\mus)')
xlim([0 22]), axis ij


%-- Demodulation and beamforming

% Demodulation
IQ = rf2iq(RF,param);

% Beamforming grid
[th,r] = meshgrid(linspace(-40,40,128)/180*pi+pi/2,...
    linspace(1,9,256)*1e-2);
[x,z] = pol2cart(th,r);

% Beamformed IQ
IQb = das(IQ,x,z,dels,param);

% Beamformed image
subplot(122)
pcolor(x*100,z*100,abs(IQb).^.5)
colormap(gray)
shading interp, axis equal ij tight
title('Gamma-compressed image')
xlabel('[cm]'), ylabel('[cm]')
hold on
% position of the scatterers (in cm)
plot(xs*100,zs*100,'ro')

end


