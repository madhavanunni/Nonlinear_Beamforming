function varargout = iq2doppler(IQ,param,M,lag)

%IQ2DOPPLER   Convert I/Q data to color Doppler
%   VD = IQ2DOPPLER(IQ,PARAM) returns the Doppler velocities from the I/Q
%   time series using a slow-time autocorrelator.
%
%   PARAM is a structure that must contain the following fields:
%        a) PARAM.fc: center frequency (in Hz, REQUIRED)
%        b) PARAM.c: longitudinal velocity (in m/s, default = 1540 m/s)
%        c) PARAM.PRF (in Hz) or PARAM.PRP (in s):
%                pulse repetition frequency or period (REQUIRED)
%
%   VD = IQ2DOPPLER(IQ,PARAM,M):
%   - If M is of a two-component vector [M(1) M(2)], the output Doppler
%     velocity is estimated from the M(1)-by-M(2) neighborhood around the
%     corresponding pixel.
%   - If M is a scalar, then an M-by-M neighborhood is used.
%   - If M is empty, then M = 1.
%
%   VD = IQ2DOPPLER(IQ,PARAM,M,LAG) uses a lag of value LAG in the
%   autocorrelator. By default, LAG = 1.
%
%   [VD,VarD] = IQ2DOPPLER(...) also returns an estimated Doppler variance.
%
%   Important note:
%   --------------
%   IQ must be a 3-D complex array, where the real and imaginary parts
%   correspond to the in-phase and quadrature components, respectively. The
%   3rd dimension corresponds to the slow-time axis. IQ2DOPPLER uses a full
%   ensemble length to perform the auto-correlation, i.e. ensemble length
%   (or packet size) = size(IQ,3).
%
%
%   REFERENCE
%   ---------
%   If you find this function useful, you can cite the following paper.
%   Key references are included in the text of the function.
%
%   1) Madiena C, Faurie J, Porée J, Garcia D, Color and vector flow
%   imaging in parallel ultrasound with sub-Nyquist sampling. IEEE Trans
%   Ultrason Ferroelectr Freq Control, 2018;65:795-802.
%   <a
%   href="matlab:web('https://www.biomecardio.com/publis/ieeeuffc18a.pdf')">download PDF</a>
%
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also RF2IQ, WFILT.
%
%   -- Damien Garcia & Jonathan Porée -- 2015/01, last update: 2020/06/24
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


narginchk(2,4);
%-
if nargin==2 || isempty(M), M = [1 1]; end
if isscalar(M), M = [M M]; end
assert(all(M>0) && isequal(M,round(M)),'M must contain integers >0')
%-
if ndims(IQ)==4
    error('IQ is a 4-D array: use IQ2DOPPLER3.')
end
assert(ndims(IQ)==3,'IQ must be a 3-D array');
%-
if nargin<4
    lag = 1;
else
    assert(isscalar(lag) && lag>0 && lag==round(lag),...
        'The 4th input parameter LAG must be a positive integer')
end

%----- Input parameters in PARAM -----
%-- 1) Speed of sound
if ~isfield(param,'c')
    param.c = 1540; % longitudinal velocity in m/s
end
c = param.c;
%-- 2) Center frequency
if isfield(param,'fc')
    fc = param.fc;
else
    error(['A center frequency (fc) ',...
        'must be specified in the structure PARAM: PARAM.fc'])
end
%-- 3) Pulse repetition frequency or period (PRF or PRP)
if isfield(param,'PRF')
    PRF = param.PRF;
elseif isfield(param,'PRP')
    PRF = 1./param.PRP;
else
    error(['A pulse repetition frequency or period ',...
        'must be specified in the structure PARAM: PARAM.PRF or PARAM.PRP'])
end
if isfield(param,'PRP') && isfield(param,'PRF')
    assert(abs(param.PRF-1./param.PRP)<eps,...
        ['A conflict exists for the pulse repetition frequency & period:',13,...
        'PARAM.PRF and 1/PARAM.PRP are different!'])
end


%--- AUTO-CORRELATION METHOD ---
% Eq. 55 in Loupas et al. (IEEE UFFC 42,4;1995)
IQ1 = IQ(:,:,1:1:end-lag);
IQ2 = IQ(:,:,1+lag:1:end);

AC = sum(IQ1.*conj(IQ2),3); % ensemble auto-correlation

if ~isequal([M(1) M(2)],[1 1]) % spatial weighted average
    h = hamming(M(1))*hamming(M(2))';
    AC = imfilter(AC,h,'replicate');
end

%-- Doppler velocity
VN = c*PRF/4/fc/lag; % Nyquist velocity
varargout{1} = -VN*imag(log(AC))/pi;


%-- Doppler variance
if nargout==2
    P = sum(real(IQ).^2+imag(IQ).^2,3); % power
    if ~isequal([M(1) M(2)],[1 1]) % spatial weighted average
        P = imfilter(P,h,'replicate');
    end
    varargout{2} = 2*(c*PRF/4/fc/lag/pi)^2*(1-abs(AC)./P);
    %-- cf. Eq. 7.48 in Estimation of Blood Velocities Using Ultrasound:
    %   A Signal Processing Approach by Jørgen Arendt Jensen,
    %   Cambridge University Press, 1996
end






