function [x,z,z0] = impolgrid(siz,zmax,width,param)

%IMPOLGRID   Polar-type grid for ultrasound images
%   IMPOLGRID returns a polar-type (fan-type) grid expressed in Cartesian
%   coordinates. This is a "natural" grid (before scan-conversion) used
%   when beamforming signals obtained with a cardiac phased array or a
%   convex array.
%
%   [X,Z] = IMPOLGRID(SIZ,ZMAX,WIDTH,PARAM) returns the X,Z coordinates of
%   the fan-type grid of size SIZ and angular width WIDTH (in rad) for a
%   phased array described by PARAM. The maximal Z (maximal depth) is ZMAX.
%
%   [X,Z] = IMPOLGRID(SIZ,ZMAX,PARAM) returns the X,Z coordinates of
%   the fan-type grid of size SIZ and angular width WIDTH (in rad) for a
%   convex array described by PARAM. For a convex array, PARAM.radius is
%   not Inf. The maximal Z (maximal depth) is ZMAX.
%
%   If SIZ is a scalar M, then the size of the grid is [M,M].
%
%   [X,Z,Z0] = IMPOLGRID(...) also returns the z-coordinate of the grid
%   origin. Note that X0 = 0.
%
%   Units: X,Z,Z0 are in m. WIDTH must be in rad.
%
%   PARAM is a structure which must contain the following fields:
%   ------------------------------------------------------------
%   1) PARAM.pitch: pitch of the array (in m, REQUIRED)
%   2) PARAM.Nelements: number of elements in the transducer array (REQUIRED)
%   3) PARAM.radius: radius of curvature (in m, default = Inf, linear array)
%
%
%   Examples:
%   --------
%   %-- Generate a focused pressure field with a phased-array transducer
%   % Phased-array @ 2.7 MHz:
%   param = getparam('P4-2v');
%   % Focus position:
%   xf = 2e-2; zf = 5e-2;
%   % TX time delays:
%   dels = txdelay(xf,zf,param);
%   % 60-degrees wide grid:
%   [x,z] = impolgrid([100 50],10e-2,pi/3,param);
%   % RMS pressure field:
%   P = pfield(x,z,dels,param);
%   % Scatter plot of the pressure field:
%   figure
%   scatter(x(:)*1e2,z(:)*1e2,5,20*log10(P(:)/max(P(:))),'filled')
%   colormap jet, axis equal ij tight
%   xlabel('cm'), ylabel('cm')
%   caxis([-20 0])
%   c = colorbar;
%   c.YTickLabel{end} = '0 dB';
%   % Image of the pressure field:
%   figure
%   pcolor(x*1e2,z*1e2,20*log10(P/max(P(:))))
%   shading interp
%   colormap hot, axis equal ij tight
%   xlabel('[cm]'), ylabel('[cm]')
%   caxis([-20 0])
%   c = colorbar;
%   c.YTickLabel{end} = '0 dB';
%
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also DAS, DASMTX, PFIELD.
%
%   -- Damien Garcia -- 2020/05, last update: 2021/01/09
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


assert(nargin>2,'Not enough input arguments.')
assert(nargin<5,'Too many input arguments.')
if nargin==3
    param = width;
end

assert(numel(siz)==1 || numel(siz)==2,'SIZ must be [M,N] or M.')
if numel(siz)==1; siz = [siz siz]; end
assert(all(siz>0) && all(siz==round(siz)),...
    'SIZ components must be positive integers.')

assert(isscalar(zmax) && zmax>0,...
    'ZMAX must be a positive scalar.')

assert(isstruct(param),'PARAM must be a structure.')

%-- Pitch (in m)
if ~isfield(param,'pitch')
    error('A pitch value (PARAM.pitch) is required.')
end
p = param.pitch;

%-- Number of elements
if isfield(param,'Nelements')
    N = param.Nelements;
else
    error('The number of elements (PARAM.Nelements) is required.')
end

%-- Radius of curvature (in m)
% for a convex array
if ~isfield(param,'radius')
    param.radius = Inf; % default = linear array
end
R = param.radius;
isLINEAR = isinf(R);

if ~isLINEAR && nargin==4
        warning('MUST:impolgrid',...
            'The parameter WIDTH is ignored with a convex array.')
end

%-- Origo (x0,z0)
% x0 = 0;
if isLINEAR
    L = (N-1)*p; % array width
    z0 = -L/2*(1+cos(width))/sin(width);
else
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    % https://en.wikipedia.org/wiki/Circular_segment
    z0 = -d;
end

%-- Image polar grid
if isLINEAR
    R = hypot(L/2,z0);
    [th,r] = meshgrid(...
        linspace(width/2,-width/2,siz(2))+pi/2,...
        linspace(R+2*p,-z0+zmax,siz(1)));
    [x,z] = pol2cart(th,r);
else
    [th,r] = meshgrid(...
        linspace(atan2(L/2,d),atan2(-L/2,d),siz(2))+pi/2,...
        linspace(R+2*p,-z0+zmax,siz(1)));
    [x,z] = pol2cart(th,r);
end
z = z+z0;

