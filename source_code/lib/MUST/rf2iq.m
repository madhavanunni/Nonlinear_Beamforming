function [IQ,Fc] = rf2iq(RF,Fs,Fc,B)

%RF2IQ   I/Q demodulation of RF data
%   IQ = RF2IQ(RF,Fs,Fc) demodulates the radiofrequency (RF) bandpass
%   signals and returns the Inphase/Quadrature (I/Q) components. IQ is a
%   complex whose real (imaginary) part contains the inphase (quadrature)
%   component.
%       1) Fs is the sampling frequency of the RF signals (in Hz),
%       2) Fc represents the center frequency (in Hz).
%
%   IQ = RF2IQ(RF,Fs) or IQ = RF2IQ(RF,Fs,[],...) calculates the carrier
%   frequency.
%   IMPORTANT: Fc must be given if the RF signal is undersampled (as in
%   bandpass sampling).
%
%   [IQ,Fc] = RF2IQ(...) also returns the carrier frequency (in Hz).
%
%   RF2IQ uses a downmixing process followed by low-pass filtering. The
%   low-pass filter is determined by the normalized cut-off frequency Wn.
%   By default Wn = min(2*Fc/Fs,0.5). The cut-off frequency Wn can be
%   adjusted if the relative bandwidth (in %) is given:
%
%   IQ = RF2IQ(RF,Fs,Fc,B)
%
%   The bandwidth in % is defined by:
%     B = Bandwidth_in_% = Bandwidth_in_Hz*(100/Fc).
%   When B is an input variable, the cut-off frequency is
%     Wn = Bandwidth_in_Hz/Fs, i.e:
%     Wn = B*(Fc/100)/Fs. 
%       
%   If there is a time offset, use PARAM.t0, as explained below.
%
%   An alternative syntax for RF2IQ is the following:
%   IQ = RF2IQ(RF,PARAM), where the structure PARAM must contain the
%   required parameters:
%       1) PARAM.fs: sampling frequency (in Hz, REQUIRED)
%       2) PARAM.fc: center frequency (in Hz, OPTIONAL, required for
%            undersampled RF signals)
%       3) PARAM.bandwidth: fractional bandwidth (in %, OPTIONAL)
%       4) PARAM.t0: time offset (in s, OPTIONAL, default = 0)
%
%   Notes on Undersampling (sub-Nyquist sampling)
%   ----------------------
%   If the RF signal is undersampled, the carrier frequency Fc must be
%   specified. If a fractional bandwidth (B or PARAM.bandwidth) is given, a
%   warning message appears if harmful aliasing is suspected.
%
%   Notes:
%   -----
%   RF2IQ treats the data along the first non-singleton dimension as
%   vectors, i.e. RF2IQ demodulates along columns for 2-D and 3-D RF data.
%   Each column corresponds to a single RF signal over (fast-) time.
%   Use IQ2RF to recover the RF signals.
%
%   Method:
%   ------
%   RF2IQ multiplies RF by a phasor of frequency Fc (down-mixing) and
%   applies a fifth-order Butterworth lowpass filter using FILTFILT:
%       IQ = RF.*exp(-1i*2*pi*Fc*t);
%       [b,a] = butter(5,2*Fc/Fs);
%       IQ = filtfilt(b,a,IQ)*2;
%
%
%   Example #1: Envelope of an RF signal
%   ----------
%   % Load an RF signal sampled at 20 MHz
%   load RFsignal@20MHz.mat
%   % I/Q demodulation
%   IQ = rf2iq(RF,20e6);
%   % RF signal and its envelope
%   plot(RF), hold on
%   plot(abs(IQ),'Linewidth',1.5), hold off
%   legend({'RF signal','I/Q amplitude'})
%
%   Example #2: Demodulation of an undersampled RF signal
%   ----------
%   % Load an RF signal sampled at 20 MHz
%   % (Center frequency = 5 MHz / Bandwidth = 2 MHz)
%   load RFsignal@20MHz.mat
%   % I/Q demodulation of the original RF signal
%   Fs = 20e6;
%   IQ = rf2iq(RF,Fs);
%   % Create an undersampled RF signal (sampling at Fs/5 = 4 MHz)
%   bpsRF = RF(1:5:end);
%   subplot(211), plot(1:1000,RF,1:5:1000,bpsRF,'.-')
%   title('RF signal (5 MHz array)')
%   legend({'sampled @ 20 MHz','bandpass sampled @ 4 MHz'})
%   % I/Q demodulation of the undersampled RF signal
%   Fs = 4e6; Fc = 5e6;
%   iq = rf2iq(bpsRF,Fs,Fc);
%   % Display the IQ signals
%   subplot(212), plot(1:1000,abs(IQ),1:5:1000,abs(iq),'.-')
%   title('I/Q amplitude')
%   legend({'sampled @ 20 MHz','bandpass sampled @ 4 MHz'})
%
%
%   REFERENCE
%   ---------
%   If you find this function useful, you can cite the following paper.
%
%   1) Madiena C, Faurie J, Porée J, Garcia D, Color and vector flow
%   imaging in parallel ultrasound with sub-Nyquist sampling. IEEE Trans
%   Ultrason Ferroelectr Freq Control, 2018;65:795-802.
%   <a
%   href="matlab:web('http://www.biomecardio.com/publis/ieeeuffc18a.pdf')">download PDF</a>
%
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also IQ2DOPPLER, IQ2RF, BMODE, WFILT.
%
%   -- Damien Garcia -- 2012/01, last update: 2020/05
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


%-- Check input arguments
narginchk(2,4);
assert(isreal(RF),'RF must contain real RF signals.')
t0 = 0; % default value for time offset
if nargin==2
    if isstruct(Fs) % RF2IQ(RF,PARAM)
        param = Fs;
        param = IgnoreCaseInFieldNames(param);
        assert(isfield(param,'fs'),...
            'A sampling frequency (PARAM.fs) is required.')
        Fs = param.fs;
        if isfield(param,'bandwidth')
            B = param.bandwidth;
        else
            B = [];
        end
        if isfield(param,'fc')
            Fc = param.fc;
        else
            Fc = [];
        end
        if isfield(param,'t0') % time offset
            t0 = param.t0;
        end
    else
        Fc = []; B = [];
    end
elseif nargin==3
    B = [];
end
assert(isscalar(Fs),...
    'The sampling frequency (Fs or PARAM.fs) must be a scalar.')
assert(isempty(Fc) || isscalar(Fc),...
    'The center frequency (Fc or PARAM.fc) must be [] or a scalar.')

%-- Convert to column vector (if RF is a row vector)
wasrow = isrow(RF);
if wasrow, RF = RF(:); end

%-- Time vector
nl = size(RF,1);
t = (0:nl-1)'/Fs;
assert(isnumeric(t0) && (isscalar(t0) || isvector(t0)) &&...
    (numel(t0)==1 || numel(t0)==nl),...
    'PARAM.t0 must be a numeric scalar or vector of size = size(RF,1).')
t = t+t0(:);

%-- Seek the carrier frequency (if required)
if isempty(Fc)
    % Keep a maximum of 100 randomly selected scanlines
    Nc = size(RF,2);
    if Nc<100, idx = 1:Nc; else, idx = randperm(Nc); idx = idx(1:100); end
    % Power Spectrum
    P = sum(abs(fft(RF(:,idx))).^2,2);
    P = P(1:floor(nl/2)+1);
    % Carrier frequency
    idx = sum((0:floor(nl/2))'.*P)/sum(P);
    Fc = (idx-1)*Fs/nl;
end

%-- Normalized cut-off frequency
if isempty(B)
    Wn = min(2*Fc/Fs,0.5);
else
    assert(isscalar(B),...
        'The signal bandwidth (B or PARAM.bandwidth) must be a scalar.')
    assert(B>0 && B<200,...
        'The signal bandwidth (B or PARAM.bandwidth, in %) must be within the interval of ]0,200[.')
    B = Fc*B/100; % bandwidth in Hz
    Wn = B/Fs;
end
assert(Wn>0 && Wn<=1,'The normalized cutoff frequency is not within the interval of (0,1). Check the input parameters!')

%-- Down-mixing of the RF signals
IQ = double(RF).*exp(-1i*2*pi*Fc*t);

%-- Low-pass filter
[b,a] = butter(5,Wn);
IQ = filtfilt(b,a,IQ)*2; % factor 2: to preserve the envelope amplitude

%-- Recover the initial size (if was a vector row)
if wasrow, IQ = IQ.'; end

%-- Display a warning message if harmful aliasing is suspected
if Fs<(2*Fc+B) % the RF signal is undersampled
    fL = Fc-B/2; fH = Fc+B/2; % lower and higher frequencies of the bandpass signal
    n = floor(fH/(fH-fL));
    harmlessAliasing = any(2*fH./(1:n)<=Fs & Fs<=2*fL./(0:n-1));
    if ~harmlessAliasing
        warning('RF2IQ:harmfulAliasing',...
            'Harmful aliasing is present: the aliases are not mutually exclusive!')
    end
end

end




function structArray = IgnoreCaseInFieldNames(structArray)

fieldLIST = {'fs','fc','t0','bandwidth'};

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

