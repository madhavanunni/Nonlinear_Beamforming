function varargout = vplot(x,y,u,v,varargin)

%VPLOT Vector plot.
%   VPLOT(X,Y,U,V) plots randomly distributed velocity vectors as colored
%   wedges (default shapes). The vector directions and locations are
%   determined by interpolation from the components (U,V) and coordinates
%   (X,Y). The arrays X,Y,U,V must all be the same size. By default, the
%   shape colors are determined by the velocity amplitudes.
%
%   Use VPLOT(...,'nonrandom') to discard the random distribution, i.e to
%   plot velocity vectors with components (U,V) at the points (X,Y).
%
%   VPLOT(X,Y,U,V,C) uses C to determine the shape colors. C must be of
%   same size as U and V. If C = [], the wedge colors are determined by the
%   velocity amplitudes. If C is a single character string chosen from the
%   list 'r','g','b','c','m','y','w','k', or an RGB row vector triple,
%   [r g b], the shape is filled with the constant specified color.
%
%   By default, VPLOT automatically scales the shapes, and their areas are
%   proportional to the velocity amplitudes. The colors are based on the
%   current colormap and the color scaling is determined by the range of C,
%   or by the current setting of CAXIS.
%
%   The shapes and their size can be modified by the structure OPTIONS:
%   VPLOT(X,Y,U,V,OPTIONS) or VPLOT(X,Y,U,V,C,OPTIONS)
%   ---
%   The structure OPTIONS can contain the following fields:
%       -----------------
%       OPTIONS.shape:      'wedge' (default), 'waterdrop', 'arrowhead' or
%                           'heart'
%       OPTIONS.s:          Automatically scales the shapes and then
%                           stretches them by OPTIONS.s (default = 1). Use
%                           OPTIONS.s = 0 to plot the arrows without the
%                           automatic scaling (as in QUIVER).
%       OPTIONS.size:       Must be 'area' or 'height'. Adjusts the shape
%                           size by making shape area or shape height
%                           proportional to the velocity amplitude
%                           (default = 'area').
%       OPTIONS.separate:   Level of separation between the shapes (avoid
%                           overlapping shapes). It must be a positive row
%                           couple of the form [s1 s2]. The 1st and 2nd
%                           components are related to the long and short
%                           axes of the ellipse that surrounds the shape.
%                           Large values = large separation
%                           (default = [1 1]).
%                           This option is not used with 'nonrandom'.
%       -----------------
%
%   P = VPLOT(...) returns the patch object that contains the data for all
%   the vector markers. Use P to query and modify properties of the patch
%   object after it is created. For a list of properties and descriptions,
%   see "Patch Properties" in Matlab Documentation.
%
%   NOTE:
%   ----
%   You can also define your own shape by using OPTIONS.shape = [X Y],
%   where X and Y are two column vectors of same length. A filled polygon
%   is created using the elements of X and Y as the coordinates for each
%   vertex (as in PATCH).
%
%   EXAMPLES:
%   --------
%   %--- Example #1
%   load wind
%   figure
%   vplot(x(:,:,7),y(:,:,7),u(:,:,7),v(:,:,7))
%   colormap(flipud(hot))
%   axis equal off
%
%   %--- Example #2
%   [X,Y] = meshgrid(-2:.05:2);
%   Z = X.*exp(-X.^2-Y.^2);
%   [DX,DY] = gradient(Z,.05,.05);
%   figure, subplot(121)
%   vplot(X,Y,DX,DY)
%   colormap(1-hot)
%   hold on
%   contour(X,Y,Z,'k','Linewidth',1.5)
%   axis equal, box on
%   set(gca,'XTick',[],'YTick',[])
%   % now, on a regular grid:
%   subplot(122)
%   opt.size = 'height';
%   vplot(X(1:3:end,1:3:end),Y(1:3:end,1:3:end),...
%       DX(1:3:end,1:3:end),DY(1:3:end,1:3:end),...
%       opt,'nonrandom')
%   hold on
%   contour(X,Y,Z,'k','Linewidth',1.5)
%   axis equal, box on
%   set(gca,'XTick',[],'YTick',[])
%
%   %--- Example #3
%   [x,y] = meshgrid(linspace(-1,1,32));
%   u = cos(x*pi+pi/2).*cos(y*pi);
%   v = sin(x*pi+pi/2).*sin(y*pi);
%   Q = curl(x,y,u,v);
%   opt.shape = 'waterdrop';
%   figure
%   vplot(x,y,u,v,Q,opt)
%   colormap([hot;1-hot])
%   axis equal off
%
%   %--- Example #4
%   [x,y] = meshgrid(0:0.2:2,0:0.2:2);
%   u = cos(x).*y;
%   v = sin(-x).*y;
%   opt.shape = 'heart';
%   figure
%   vplot(x,y,u,v,opt)
%   colormap(flipud(spring))
%   axis equal off
%
%   %--- Example #5
%   [x,y] = meshgrid(linspace(-1,1,32));
%   u = cos(x*pi+pi/2).*cos(y*pi);
%   v = sin(x*pi+pi/2).*sin(y*pi);
%   %-- define your own shape
%   t = linspace(0,2*pi,100)';
%   opt.shape = [1/3*(3*cos(t)+cos(3*t)) 3*sin(t)-sin(3*t)];
%   figure
%   vplot(x,y,u,v,opt)
%   colormap winter
%   axis equal off
%
%   See also QUIVER, PATCH
%
%   -- Damien Garcia -- 2016/12, update: 2019/07
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

narginchk(4,7)
nargoutchk(0,1)

assert(isequal(size(x),size(y),size(u),size(v)),...
    'X, Y, U and V must be of same size.')

if nargin==4
    C = [];
    options = [];
    isRandomlyDistributed = true;
elseif nargin==5
    if isstruct(varargin{1}) % vplot(x,y,u,v,options)
        options = varargin{1};
        C = [];
        isRandomlyDistributed = true;
    elseif strcmpi(varargin{1},'nonrandom') % vplot(x,y,u,v,'nonrandom')
        isRandomlyDistributed = false;
        C = [];
        options = [];
    else % vplot(x,y,u,v,C)
        C = varargin{1};
        options = [];
        isRandomlyDistributed = true;
    end
elseif nargin==6
    if strcmpi(varargin{2},'nonrandom') % vplot(x,y,u,v,C/options,'nonrandom')
        if isstruct(varargin{1}) % vplot(x,y,u,v,options,'nonrandom')
            options = varargin{1};
            C = [];
        else % vplot(x,y,u,v,C,'nonrandom')
            C = varargin{1};
            options = [];
        end
        isRandomlyDistributed = false;
        
    else % vplot(x,y,u,v,C,opt)
        C = varargin{1};
        options = varargin{2};
        isRandomlyDistributed = true;
    end
else
    assert(strcmpi(varargin{3},'nonrandom'),...
        'The distribution option must be ''nonrandom'' ')
    C = varargin{1};
    options = varargin{2};
    isRandomlyDistributed = false;
end

%-- Color specifications
if isscalar(C) && any(strcmpi(C,{'r','g','b','c','m','y','w','k'}))
    % color defined by a single character string
    isConstantColor = true;
elseif isnumeric(C) && isequal(size(C),[1 3])
    % color defined by an RGB triplet
    isConstantColor = true;
elseif isempty(C)
    % color defined by the vector amplitudes
    isConstantColor = false;
else
    % color defined by the array C
    assert(isequal(size(C),size(x)),...
        ['If the color scaling is determined by the range of C, ',...
        'C, X and Y must be of same size.'])
    isConstantColor = false;
    C = C(:);
end

if issparse(u), u = full(u); end
if issparse(v), v = full(v); end
V = complex(u,v);
clear u v
% vectorize
x = x(:); y = y(:); V = V(:);

%-- Default values for OPTIONS
%- Scaling
if ~isfield(options,'s')
    options.s = 1;
end
s = options.s;
assert(isscalar(s) & isnumeric(s) & s>=0,...
    'OPTIONS.s must be nonnegative (i.e. >=0)')
%- Area or Height size?
if ~isfield(options,'size')
    options.size = 'area';
end
assert(any(strcmpi(options.size,{'area','height'})),...
    'OPTIONS.size must be ''area'' or ''height''.')
%- Separation factor
if ~isfield(options,'separate')
    options.separate = [1 1];
end
assert(isnumeric(options.separate) &...
    isequal(size(options.separate),[1 2]) &...
    all(options.separate>0),...
    'OPTIONS.separate must be a positive row couple.')
%- Shape
if ~isfield(options,'shape')
    options.shape = 'wedge';
end
if isnumeric(options.shape)
    assert(ismatrix(options.shape) & size(options.shape,2)==2,...
        ['When defining a personalized shape, OPTIONS.shape must be '...
        'a two-column matrix.'])
else
    assert(any(strcmpi(options.shape,...
        {'wedge','waterdrop','arrowhead','heart','arrow'})),...
        ['OPTIONS.shape must be ',...
        '''wedge'', ''waterdrop'', ''heart'', ''arrow'' or ''arrowhead''.'])
end

%-- Advanced options
if isfield(options,'rng')
    rng(options.rng)
end
% number of initial random points
if isfield(options,'N')
    N = options.N;
else
    N = 50000;
end

% x- and y-ranges
minx = min(x); rangex = max(x)-minx;
miny = min(y); rangey = max(y)-miny;

Vmax = max(abs(V(:)));

if isRandomlyDistributed
    
    % random distribution of centers
    xi = rand(N,1)*rangex + minx;
    yi = rand(N,1)*rangey + miny;
    F = scatteredInterpolant(x,y,V,'linear','none');
    V = F(xi,yi);
    
    % discard the smallest speeds
    % (since the corresponding shapes would be barely visible)
    idx = ~isfinite(V);
    if strcmpi(options.size,'area')
        idx = idx | (abs(V) < Vmax/100);
    else
        idx = idx | (abs(V) < Vmax/10000);
    end
    xi(idx) = []; yi(idx) = [];
    V(idx) = [];
    
    % vector direction
    TH = atan2(imag(V),real(V));
    if ~isConstantColor
        if ~isempty(C)
            F = scatteredInterpolant(x,y,C,'linear','none');
            Ci = F(xi,yi);
        end
    end
    
    x = xi; clear xi
    y = yi; clear yi
    
else
    N = numel(V);
    TH = atan2(imag(V),real(V));
end

%-- Automatic scaling then stretch by S
if s>0, s = s*max(rangex,rangey)/10; end

%-- Elementary shape
% note: a = (short axis)/(long axis) ratio
if isnumeric(options.shape)
    % Personalized shape
    x0 = options.shape(:,1);
    y0 = options.shape(:,2);
    y0 = y0-min(y0);
    x0 = x0/max(y0);
    y0 = y0/max(y0);
    y0 = y0-0.5;
    a = range(x0);
elseif strcmpi(options.shape,'wedge')
    % Wedge
    a = 1/4;
    x0 = a*[1 -1 0]'/2;
    y0 = [-1 -1 1]'/2;
elseif strcmpi(options.shape,'waterdrop')
    % Waterdrop
    a = 1/4;
    t = linspace(0,2*pi,50)';
    x0 = a*(1-sin(t)).*cos(t)/(3*sqrt(3)/2);
    y0 = 0.5*(sin(t)-1)+0.5;
elseif strcmpi(options.shape,'arrow')
    % Arrow
    a = 1/4;
    x0 = a*[0 -5 -2 -2 2 2 5]'/10;
    y0 = [2 1 1 -2 -2 1 1]'/4;
elseif strcmpi(options.shape,'arrowhead')
    % Arrowhead
    a = 1/4;
    x0 = a*[-1 0 1 0]'/2;
    y0 = [-1 -2 -1 2]'/4;
    % y0 = [-2 -1 -2 2]'/4;
elseif strcmpi(options.shape,'heart')
    % Heart
    a = 1;
    t = linspace(0,2*pi,50)';
    x0 = a*sin(t).^3/2;
    y0 = (-13*cos(t)+5*cos(2*t)+2*cos(3*t)+cos(4*t)+11.9232)/28.9226-.5;
end

xp = NaN(length(x0),N);
yp = NaN(length(x0),N);

if isConstantColor
    Cp = C;
else
    Cp = NaN(N,1);
end

f = options.separate;
idx = [];


if isRandomlyDistributed
    
    for k = 1:N
        
        TH(idx) = [];
        TH1 = TH(1);
        V(idx) = [];
        speed = abs(V(1)); % speed
        
        %-- height (h) of the shape
        if strcmpi(options.size,'area')
            if s==0
                h = speed;
            else
                h = speed/Vmax*s;
            end
        elseif strcmpi(options.size,'height')
            if s==0
                h = sqrt(speed);
            else
                h = sqrt(speed/Vmax)*s;
            end
        end
        
        %-- shape coordinates
        xk = h*(-x0*sin(TH1) + y0*cos(TH1)) + x(1);
        yk = h*(x0*cos(TH1) + y0*sin(TH1)) + y(1);
        %- coordinates of the vertex (for PATCH function)
        xp(:,k) = xk; yp(:,k) = yk;
        if isConstantColor
            % do nothing
        elseif isempty(C)
            Cp(k) = speed;
        else
            Ci(idx) = [];
            Cp(k) = Ci(1);
        end
        
        %-- ellipse surrounding the shape
        xe = -(x-x(1))*sin(TH1) + (y-y(1))*cos(TH1);
        ye = (x-x(1))*cos(TH1) + (y-y(1))*sin(TH1);
        L = h*f(1); % long axis
        l = a*h*f(2); % small axis
        % Reject the points inside this ellipse (avoid overlapping shapes)
        idx = (xe.^2/l^2 + ye.^2/L^2) < 1;
        
        x(idx) = []; y(idx) = [];
        if isempty(x), break, end
        
    end
    
else
    
    speed = abs(V);
    
    %-- height (h) of the shape
    if strcmpi(options.size,'area')
        if s==0
            h = speed;
        else
            h = speed/Vmax*s;
        end
    elseif strcmpi(options.size,'height')
        if s==0
            h = sqrt(speed);
        else
            h = sqrt(speed/Vmax)*s;
        end
    end
    
    %- coordinates of the vertex (for PATCH function)
    if verLessThan('matlab','9.2')
        xp = bsxfun(@times,h',(-x0*sin(TH') + y0*cos(TH')));
        xp = bsxfun(@plus,xp,x');
        yp = bsxfun(@times,h',(x0.*cos(TH') + y0*sin(TH')));
        yp = bsxfun(@plus,yp,y');
    else % BSXFUN is no longer required
        xp = h'.*(-x0*sin(TH') + y0*cos(TH')) + x';
        yp = h'.*(x0*cos(TH') + y0*sin(TH')) + y';
    end
    
    if isConstantColor
        % do nothing
    elseif isempty(C)
        Cp = speed;
    else
        Cp = C;
    end
    
end

if nargout==1
    varargout{1} = patch(xp,yp,Cp,'EdgeColor','none');
else
    patch(xp,yp,Cp,'EdgeColor','none')
end
