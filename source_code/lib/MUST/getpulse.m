function [pulse,t] = getpulse(param,way)

%GETPULSE   Get the transmit pulse
%   PULSE = GETPULSE(PARAM,WAY) returns the one-way or two-way transmit
%   pulse with a time sampling of 1 nanosecond. Use WAY = 1 to get the
%   one-way pulse, or WAY = 2 to obtain the two-way (pulse-echo) pulse.
%
%   PULSE = GETPULSE(PARAM) uses WAY = 1.
%
%   [PULSE,t] = GETPULSE(...) also returns the time vector.
%
%   PARAM is a structure which must contain the following fields:
%   ------------------------------------------------------------
%   1) PARAM.fc: central frequency (in Hz, REQUIRED)
%   2) PARAM.bandwidth: pulse-echo 6dB fractional bandwidth (in %)
%            The default is 75%.
%   3) PARAM.TXnow: number of wavelengths of the TX pulse (default: 1)
%   4) PARAM.TXfreqsweep: frequency sweep for a linear chirp (default: [])
%                         To be used to simulate a linear TX chirp.
%
%   Example #1:
%   ----------
%   %-- Get the one-way pulse of a phased-array probe
%   % Phased-array @ 2.7 MHz:
%   param = getparam('P4-2v');
%   % One-way transmit pulse
%   [pulse,t] = getpulse(param);
%   % Plot the pulse
%   plot(t*1e6,pulse)
%   xlabel('{\mu}s')
%   axis tight
%
%   Example #2:
%   ----------
%   %-- Check the pulse with a linear chirp
%   % Linear array:
%   param = getparam('L11-5v');
%   % Modify the fractional bandwidth:
%   param.bandwidth = 120;
%   % Define the properties of the chirp
%   param.TXnow = 20;
%   param.TXfreqsweep = 10e6;
%   % One-way transmit pulse
%   [pulse,t] = getpulse(param);
%   % Plot the pulse
%   plot(t*1e6,pulse)
%   xlabel('{\mu}s')
%   axis tight
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also PFIELD, SIMUS, GETPARAM.
%
%   -- Damien Garcia -- 2020/12, last update: 2020/12/06
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

narginchk(1,2)
nargoutchk(1,2)
if nargin==1, way = 1; end % one-way is the default
assert(way==1 || way==2,'WAY must be 1 (one-way) or 2 (two-way)')

%-- Center frequency (in Hz)
assert(isfield(param,'fc'),...
    'A center frequency value (PARAM.fc) is required.')
fc = param.fc; % central frequency (Hz)

%-- Fractional bandwidth at -6dB (in %)
if ~isfield(param,'bandwidth')
    param.bandwidth = 75;
end
assert(param.bandwidth>0 && param.bandwidth<200,...
    'The fractional bandwidth at -6 dB (PARAM.bandwidth, in %) must be in ]0,200[')

%-- TX pulse: Number of wavelengths
if ~isfield(param,'TXnow')
    param.TXnow = 1;
end
NoW = param.TXnow;
assert(isscalar(NoW) && isnumeric(NoW) && NoW>0,...
    'PARAM.TXnow must be a positive scalar.')

%-- TX pulse: Frequency sweep for a linear chirp
if ~isfield(param,'TXfreqsweep') || isinf(NoW)
    param.TXfreqsweep = [];
end
FreqSweep = param.TXfreqsweep;
assert(isempty(FreqSweep) ||...
    (isscalar(FreqSweep) && isnumeric(FreqSweep) && FreqSweep>0),...
    'PARAM.TXfreqsweep must be empty (windowed sine) or a positive scalar (linear chirp).')

mysinc = @(x) sinc(x/pi); % cardinal sine
% [note: In MATLAB, sinc is sin(pi*x)/(pi*x)]

%-- FREQUENCY SPECTRUM of the transmitted pulse
if isempty(FreqSweep)
    % We want a windowed sine of width PARAM.TXnow
    T = NoW/fc; % temporal pulse width
    wc = 2*pi*fc;
    pulseSpectrum = @(w) 1i*(mysinc(T*(w-wc)/2)-mysinc(T*(w+wc)/2));
else
    % We want a linear chirp of width PARAM.TXnow
    % (https://en.wikipedia.org/wiki/Chirp_spectrum#Linear_chirp)
    T = NoW/fc; % temporal pulse width
    wc = 2*pi*fc;
    dw = 2*pi*FreqSweep;
    s2 = @(w) sqrt(pi*T/dw)*exp(-1i*(w-wc).^2*T/2/dw).*...
        (fresnelint((dw/2+w-wc)/sqrt(pi*dw/T)) +...
        fresnelint((dw/2-w+wc)/sqrt(pi*dw/T)));
    pulseSpectrum = @(w) (1i*s2(w)-1i*s2(-w))/T;
end

%-- FREQUENCY RESPONSE of the ensemble PZT + probe
% We want a generalized normal window (6dB-bandwidth = PARAM.bandwidth)
% (https://en.wikipedia.org/wiki/Window_function#Generalized_normal_window)
wB = param.bandwidth*wc/100; % angular frequency bandwidth
p = log(126)/log(2*wc/wB); % p adjusts the shape
probeSpectrum = @(w) exp(-(abs(w-wc)/(wB/2/log(2)^(1/p))).^p);
% The frequency response is a pulse-echo (transmit + receive) response. A
% square root is thus required when calculating the pressure field:
probeSpectrum = @(w) sqrt(probeSpectrum(w));
% Note: The spectrum of the pulse (pulseSpectrum) will be then multiplied
% by the frequency-domain tapering window of the transducer (probeSpectrum)

%-- frequency samples
dt = 1e-9; % time step is 1 ns
df = param.fc/param.TXnow/32;
p = nextpow2(1/dt/2/df);
Nf = 2^p;
f = linspace(0,1/dt/2,Nf);

%-- spectrum of the pulse
F = pulseSpectrum(2*pi*f).*probeSpectrum(2*pi*f).^way;

%-- pulse in the temporal domain (step = 1 ns)
tmp = [F conj(F(end-1:-1:2))];
pulse = fftshift(ifft(tmp,'symmetric'));
pulse = pulse/max(abs(pulse));

%-- keep the significant magnitudes
idx1 = find(pulse>(1/1023),1);
idx2 = find(pulse>(1/1023),1,'last');
idx = min(idx1,2*Nf-1-idx2);
pulse = pulse(end-idx+1:-1:idx);

%-- time vector
if nargout==2
    t = (0:length(pulse)-1)*dt;
end

end


function f = fresnelint(x)

% FRESNELINT Fresnel integral.
%
% J = FRESNELINT(X) returns the Fresnel integral J = C + 1i*S.
%
% We use the approximation introduced by Mielenz in
%       Klaus D. Mielenz, Computation of Fresnel Integrals. II
%       J. Res. Natl. Inst. Stand. Technol. 105, 589 (2000), pp 589-590
%

siz0 = size(x);
x = x(:);

issmall = abs(x)<=1.6;
c = zeros(size(x));
s = zeros(size(x));

% When |x| < 1.6, a Taylor series is used (see Mielenz's paper)
if any(issmall)
    n = 0:10;
    cn = [1 cumprod(-pi^2*(4*n+1)./(4*(2*n+1).*(2*n+2).*(4*n+5)))];
    sn = [1 cumprod(-pi^2*(4*n+3)./(4*(2*n+2).*(2*n+3).*(4*n+7)))]*pi/6;
    n = [n 11];
    c(issmall) = sum(cn.*x(issmall).^(4*n+1),2);
    s(issmall) = sum(sn.*x(issmall).^(4*n+3),2);    
end

% When |x| > 1.6, we use the following:
if any(~issmall)
    n = 0:11;
    fn = [0.318309844, 9.34626e-8, -0.09676631, 0.000606222, ...
        0.325539361, 0.325206461, -7.450551455, 32.20380908, ...
        -78.8035274, 118.5343352, -102.4339798, 39.06207702];
    gn = [0, 0.101321519, -4.07292e-5, -0.152068115, -0.046292605, ...
        1.622793598, -5.199186089, 7.477942354, -0.695291507, ...
        -15.10996796, 22.28401942, -10.89968491];
    fx = sum(fn.*x(~issmall).^(-2*n-1),2);
    gx = sum(gn.*x(~issmall).^(-2*n-1),2);    
    c(~issmall) = 0.5*sign(x(~issmall)) + ...
        fx.*sin(pi/2*x(~issmall).^2) - gx.*cos(pi/2*x(~issmall).^2);
    s(~issmall) = 0.5*sign(x(~issmall)) - ...
        fx.*cos(pi/2*x(~issmall).^2) - gx.*sin(pi/2*x(~issmall).^2);
end

f = reshape(c,siz0) + 1i*reshape(s,siz0);

end
