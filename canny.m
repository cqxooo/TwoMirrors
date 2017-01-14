function E = canny(I, S, U, L)

% E = canny(I, S, U, L) performs a canny edge detection on the input
% image I (sy x sx) but convolving with a gaussian filter with
% sigma S and finding local maxima of the intensity gradients.
% Hystersis is then applied to the maxima located with an upper
% threshold U and a lower threshold L.
% Default S = 1, U = 30, L = 10.  
%
% For an image of size (sy,sx):
% X-extend is from 0.5 to sx.5
% Y-extend is from 0.5 to sy.5
%          X
%     0-------->
%     |       
%     |
%   Y |
%     |
%     v 
%
% The output of the canny edge detector is a matrix E
% of the following format:
%   E(:,:,1) = edge image : 1 = weak edgel, 
%                           2 = connected weak edgel and
%                           3 = strong edgel
%   E(:,:,2) = sub-pixel x-coordinates
%   E(:,:,3) = sub-pixel y-coordinates
%   E(:,:,4) = gradients
%   E(:,:,5) = orientations
%
 
% written by Kenneth Wong.
% last modified on 12-Feb-2000.
  
switch(nargin)
 case 1
  S = 1;
  U = 30;
  L = 10;
 case 2
  U = 30;
  L = 10;
 case 3
  if S<=0
    S = 1;
  end
  if U<=0
    U = 1;
  end
  L = U/3;
end

  
[sy, sx] = size(I);
E = zeros(sy,sx,5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Smooth the image with a gaussian filter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma = S;
sig2 = 2*sigma*sigma;
GaussianDieOff = 0.0001;

% possible widths
pw = 1:50;
width = max(find(exp(-(pw.*pw)/sig2)>GaussianDieOff));
if isempty(width)
   width = 1;
end

% define the gaussian filter with suitable size
t = (-width:width);
gau = exp(-(t.*t)/sig2);
scale = sum(gau);
gau = gau./scale;

% convolve the image with the gaussian filter in x and y directions
Is = conv2(gau, gau, I, 'same');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Find the gradient of the smpoohed image %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the gradient vectors
Ix = conv2(Is, [-1 0 1] , 'same');
Iy = conv2(Is, [-1 0 1]', 'same');

% find the magitudes and orientations of the gradient vectors
mag = sqrt(Ix.^2 + Iy.^2);
direc = atan2(Iy, Ix);

E(:,:,4) = mag;
E(:,:,5) = direc;

% normalize
%magmax = max(mag(:));
%if magmax > 0
%   mag = mag / magmax;   
%end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Non-maximal suppression %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                         The X marks the pixel in question, and each
%        (2)   (3)        of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45 
%    (1)|         |(4)     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight 
%       |         |       divisions, but for the non-maximum supression  
%     4 |         | 1     we are only worried about 4 of them since we 
%       O----O----O       use symmetric points about the center pixel.
%         3     2         

idxLocalMax = [];
XMap = E(:,:,2);
YMap = E(:,:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% guadrant 1
idx = find(((direc>0 & direc<=pi/4) | ...
   (direc>-pi & direc<=-3*pi/4)) & mag>L);

if ~isempty(idx)
   v = mod(idx, sy);
   extIdx = find(v==1 | v==0 | idx<=sy | (idx>(sx-1)*sy));
   idx(extIdx) = [];
end

if ~isempty(idx)
   ixv = Ix(idx);
   iyv = Iy(idx);
   gradmag = mag(idx);
   
   idx4 = idx-sy;
   idx6 = idx+sy;

   % Do the linear interpolations
   d = iyv./ixv;
   gradmag1 = mag(idx4).*(1-d) + mag(idx4-1).*d;
   gradmag2 = mag(idx6).*(1-d) + mag(idx6+1).*d; 

   max_ = find(gradmag>gradmag1 & gradmag>=gradmag2);
  
   idx = idx(max_);
   d = d(max_);
   gradmag = gradmag(max_);
   gradmag1 = gradmag1(max_);
   gradmag2 = gradmag2(max_);

   idxLocalMax = [idxLocalMax; idx];
   sec_direc = sqrt(1+d.^2);
   delta = sec_direc.*(gradmag1-gradmag2)./(2*(gradmag1+gradmag2-2*gradmag));
   
   dx = delta./sec_direc;
   XMap(idx) = floor((idx-1)/sy)+1 + dx;
   YMap(idx) = cycle(idx,sy) + dx.*d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% guadrant 2
idx = find(((direc>pi/4 & direc<=pi/2) | ...
   (direc>-3*pi/4 & direc<=-pi/2)) & mag>L);

if ~isempty(idx)
   v = mod(idx, sy);
   extIdx = find(v==1 | v==0 | idx<=sy | (idx>(sx-1)*sy));
   idx(extIdx) = [];
end

if ~isempty(idx)
   ixv = Ix(idx);
   iyv = Iy(idx);
   gradmag = mag(idx);
   
   idx2 = idx+1;
   idx8 = idx-1;

   % Do the linear interpolations
   d = ixv./iyv;
   gradmag1 = mag(idx8).*(1-d) + mag(idx8-sy).*d;
   gradmag2 = mag(idx2).*(1-d) + mag(idx2+sy).*d; 

   max_ = find(gradmag>gradmag1 & gradmag>=gradmag2);
  
   idx = idx(max_);
   d = d(max_);
   gradmag = gradmag(max_);
   gradmag1 = gradmag1(max_);
   gradmag2 = gradmag2(max_);

   idxLocalMax = [idxLocalMax; idx];
   sec_direc = sqrt(1+d.^2);
   delta = sec_direc.*(gradmag1-gradmag2)./(2*(gradmag1+gradmag2-2*gradmag));
   
   dy = delta./sec_direc;
   XMap(idx) = floor((idx-1)/sy)+1 + dy.*d;
   YMap(idx) = cycle(idx,sy) + dy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% guadrant 3
idx = find(((direc>pi/2 & direc<=3*pi/4) | ...
   (direc>-pi/2 & direc<=-pi/4)) & mag>L);

if ~isempty(idx)
   v = mod(idx, sy);
   extIdx = find(v==1 | v==0 | idx<=sy | (idx>(sx-1)*sy));
   idx(extIdx) = [];
end

if ~isempty(idx)
   ixv = Ix(idx);
   iyv = Iy(idx);
   gradmag = mag(idx);
   
   idx2 = idx+1;
   idx8 = idx-1;

   % Do the linear interpolations
   d = -ixv./iyv;
   gradmag1 = mag(idx8).*(1-d) + mag(idx8+sy).*d; 
   gradmag2 = mag(idx2).*(1-d) + mag(idx2-sy).*d; 

   max_ = find(gradmag>gradmag1 & gradmag>=gradmag2);
  
   idx = idx(max_);
   d = d(max_);
   gradmag = gradmag(max_);
   gradmag1 = gradmag1(max_);
   gradmag2 = gradmag2(max_);

   idxLocalMax = [idxLocalMax; idx];
   sec_direc = sqrt(1+d.^2);
   delta = sec_direc.*(gradmag1-gradmag2)./(2*(gradmag1+gradmag2-2*gradmag));
   
   dy = delta./sec_direc;
   XMap(idx) = floor((idx-1)/sy)+1 - dy.*d;
   YMap(idx) = cycle(idx,sy) + dy;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% guadrant 4
idx = find(((direc>3*pi/4 & direc<=pi) | ...
   (direc>-pi/4 & direc<=0)) & mag>L);

if ~isempty(idx)
   v = mod(idx, sy);
   extIdx = find(v==1 | v==0 | idx<=sy | (idx>(sx-1)*sy));
   idx(extIdx) = [];
end

if ~isempty(idx)
   ixv = Ix(idx);
   iyv = Iy(idx);
   gradmag = mag(idx);
   
   idx4 = idx-sy;
   idx6 = idx+sy;

   % Do the linear interpolations
   d = -iyv./ixv;
   gradmag1 = mag(idx6).*(1-d) + mag(idx6-1).*d; 
   gradmag2 = mag(idx4).*(1-d) + mag(idx4+1).*d;

   max_ = find(gradmag>gradmag1 & gradmag>=gradmag2);
  
   idx = idx(max_);
   d = d(max_);
   gradmag = gradmag(max_);
   gradmag1 = gradmag1(max_);
   gradmag2 = gradmag2(max_);

   idxLocalMax = [idxLocalMax; idx];
   sec_direc = sqrt(1+d.^2);
   delta = sec_direc.*(gradmag1-gradmag2)./(2*(gradmag1+gradmag2-2*gradmag));
   
   dx = delta./sec_direc;
   XMap(idx) = floor((idx-1)/sy)+1 - dx;
   YMap(idx) = cycle(idx,sy) + dx.*d;
end

E(:,:,2) = XMap;
E(:,:,3) = YMap;

if isempty(idxLocalMax)
   warning('No edge detected!');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Threshold the edgels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E(idxLocalMax) = 1;
idxStrong = idxLocalMax(mag(idxLocalMax)>U);
E(idxStrong) = 3;


%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Hysteresis %
%%%%%%%%%%%%%%%%%%%%%%
EMap = E(:,:,1);

while (1)
   idx = find(EMap == 1);
   
   if ~isempty(idx)
      % 4-neighbours
      N4 = EMap(idx-sy)>1 | EMap(idx+sy)>1 | EMap(idx-1)>1 | EMap(idx+1)>1;
      % diagonal-neighbours
      ND = EMap(idx-sy-1)>1 | EMap(idx-sy+1)>1 | EMap(idx+sy-1)>1 | EMap(idx+sy+1)>1; 
      % look for weak edgels 8-connected to strong edgels/connected weak edgels
      idxConnected = idx(N4==1 | ND==1);
      if ~isempty(idxConnected)
         EMap(idxConnected) = 2;
      else
         break;
      end
   else
      break;
   end
end

E(:,:,1) = EMap;

return
