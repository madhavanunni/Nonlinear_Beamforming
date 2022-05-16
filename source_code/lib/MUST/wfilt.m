function SIG = wfilt(SIG,method,n)

%WFILT   Wall filtering (or clutter filtering)
%   fSIG = WFILT(SIG,METHOD,N) high-pass (wall) filters the RF or I/Q
%   signals stored in the 3-D array SIG for Doppler imaging.
%    
%   The first dimension of SIG (i.e. each column) corresponds to a single
%   RF or I/Q signal over (fast-) time, with the first column corresponding
%   to the first transducer element. The third dimension corresponds to the
%   slow-time axis.
%
%   Three methods are available.
%   METHOD can be one of the following (case insensitive):
%
%   1) 'poly' - Least-squares (Nth degree) polynomial regression.
%               Orthogonal Legendre polynomials are used. The fitting
%               polynomial is removed from the original I/Q or RF data to
%               keep the high-frequency components. N (>=0) represents the
%               degree of the polynomials. The (slow-time) mean values are
%               removed if N = 0 (the polynomials are reduced to
%               constants).
%   2) 'dct'  - Truncated DCT (Discrete Cosine Transform).
%               Discrete cosine transforms (DCT) and inverse DCT are
%               performed along the slow-time dimension. The signals are
%               filtered by withdrawing the first N (>=1) components, i.e.
%               those corresponding to the N lowest frequencies (with
%               respect to slow-time).
%   3) 'svd'  - Truncated SVD (Singular Value Decomposition).
%               An SVD is carried out after a column arrangement of the
%               slow-time dimension. The signals are filtered by
%               withdrawing the top N singular vectors, i.e. those
%               corresponding to the N greatest singular values.
%   
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also IQ2DOPPLER, RF2IQ.
%
%   -- Damien Garcia -- 2014/06, last update 2020/06
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>



%-- Check the input arguments
narginchk(3,3);
assert(ndims(SIG)==3,'SIG must be a 3-D array with size(SIG,3)>=2');
assert(isscalar(n) && isnumeric(n) && isequal(abs(n),round(n)),...
    'N must be a nonnegative integer.')

siz0 = size(SIG);
N = siz0(3); % number of slow-time samples


switch lower(method)
    
    case 'poly'
        % ---------------------------------
        % POLYNOMIAL REGRESSION WALL FILTER
        % ---------------------------------
        
        % If the degree is 0, the mean is removed.
        if n==0
            SIG = SIG-mean(SIG,3);
            return
        end
        
        % We will work along columns. Each column represents a single
        % slow-time ensemble of N samples.
        SIG = reshape(shiftdim(SIG,2),N,[]);
        
        % We first create the orthogonal polynomial (Legendre) family:
        % PL contains the orthogonal polynomials i.e.
        % PL(:,1,k) = polynomial or order (k-1)
        PL = ones(N,1,n(end)+1);
        x = linspace(-1,1,N)';
        % PL(:,1,2) = x-mean(x);
        PL(:,1,2) = x; % we have mean(x) = 0
        for i = 3:n(end)+1
            tmp1 = sum(x.*PL(:,1,i-1).^2)/sum(PL(:,1,i-1).^2);
            tmp2 = sum(x.*PL(:,1,i-2).*PL(:,1,i-1))/sum(PL(:,1,i-2).^2);
            PL(:,1,i) = (x-tmp1).*PL(:,1,i-1)-tmp2*PL(:,1,i-2);
        end
        
        % ------
        % Notes:
        % 1) If we want the Legendre polynomials to range within
        %    [-1,1], use: PL = PL./PL(N,1,:);
        % 2) PL is of size:
        %    (number of slow-time samples)-by-1-by-(order+1)
        % 3) C = coefficients associated to the polynomials
        %   (least-squares sense). y = sum(P.*C,2) contains the fitting
        %   polynomials.
        % ------
        
        C = sum(PL.*SIG)./sum(PL.^2); % polynomial coefficients
        lfSIG = sum(PL.*C,3); % low-pass filtered SIG
        SIG = SIG - lfSIG; % high-pass filtered SIG
        SIG = reshape(permute(SIG,[3 2 1]),siz0);
        
        
    case 'dct'
        % -------------------------------------
        % DISCRETE COSINE TRANSFORM WALL FILTER
        % -------------------------------------
        
        assert(n>0,'N must be >0 with the ''dct'' method.')
        
        % If the degree is 0, the mean is removed.
        if n==1
            SIG = SIG-mean(SIG,3);
            return
        end
        
        % We will work along columns. Each column represents a single
        % slow-time ensemble of N samples.
        SIG = reshape(shiftdim(SIG,2),N,[]);
        
        D = dctmtx(N); % DCT matrix
        % lfSIG = D(1:n,:)'*D(1:n,:)*SIG;  % low-pass filtered SIG
        % SIG = SIG - lfSIG; % high-pass filtered SIG
        SIG = D(n+1:N,:)'*D(n+1:N,:)*SIG; % high-pass filtered SIG
        SIG = reshape(permute(SIG,[3 2 1]),siz0);
        
        
    case 'svd'
        % ----------------------------------------
        % SINGULAR VALUE DECOMPOSITION WALL FILTER
        % ----------------------------------------
        
        assert(n>0,'N must be >0 with the ''svd'' method.')
        
        % Each column represents a column-rearranged frame.
        SIG = reshape(SIG,[],N);
        
        [U,S,V] = svd(SIG,'econ'); % SVD decomposition
        % lfSIG = U(:,1:n)*S(1:n,1:n)*V(:,1:n)'; % low-pass filtered SIG
        % SIG = SIG - lfSIG; % high-pass filtered SIG
        SIG = U(:,n+1:N)*S(n+1:N,n+1:N)*V(:,n+1:N)'; % high-pass filtering
        SIG = reshape(SIG,siz0);
        
        
    otherwise
        error('METHOD must be ''poly'', ''dct'', or ''svd''.')
end


    
