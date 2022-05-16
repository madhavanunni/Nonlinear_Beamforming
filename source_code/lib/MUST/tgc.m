function [S,C] = tgc(S)

%TGC Time-gain compensation for RF or IQ signals
%   TGC(RF) or TGC(IQ) performs a time-gain compensation of the RF or IQ
%   signals using a decreasing exponential law. Each column of the RF/IQ
%   array must correspond to a single RF/IQ signal over (fast-) time.
%
%   [~,C] = TGC(RF) or [~,C] = TGC(IQ) also returns the coefficients used
%   for time-gain compensation (i.e. new_SIGNAL = C.*old_SIGNAL)
%
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also RF2IQ, DAS.
%
%   -- Damien Garcia -- 2012/10, last update 2020/05
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


siz0 = size(S);
if isvector(S), S = S(:); end
S = reshape(S,size(S,1),[]);

if isreal(S) % we have RF signals
    C = mean(abs(hilbert(S)),2);
    % C = median(abs(hilbert(S)),2);
else % we IQ signals
    C = mean(abs(S),2);
    % C = median(abs(S),2);
end
n = length(C);
n1 = ceil(n/10);
n2 = floor(n*9/10);

% -- Robust linear fitting of log(C)
% The intensity is assumed to decrease exponentially as distance increases.
% A robust linear fitting is performed on log(C) to seek the TGC
% exponential law.
% --
% See RLINFIT for details
N = 200; % a maximum of N points is used for the fitting
p = min(N/(n2-n1)*100,100);
[slope,intercept] = rlinfit((n1:n2)',log(C(n1:n2)),p);

C = exp(intercept+slope*(1:n)');
C = C(1)./C;
S = S.*C;

S = reshape(S,siz0);
if isvector(S), C = reshape(C,siz0); end

end

function [slope,intercept] = rlinfit(x,y,p)

%RLINFIT   Robust linear regression
%   See the original RLINFIT function for details
N = numel(x);
I = randperm(N);
n = round(N*p/100);
I = I(1:n);
x = x(I); y = y(I);
C = combnk(1:n,2);
slope = median((y(C(:,2))-y(C(:,1)))./(x(C(:,2))-x(C(:,1))));
intercept = median(y-slope*x);
end

