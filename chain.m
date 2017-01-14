function C = chain(E, L, T, P)

% C = chain(E, L, T, P) links the edgels resulted from the canny detector
% to form open and closed chains C. T is the maximium deviation of 
% the angle in degree from 90 between the gradient vector and the 
% vector formed by the 2 linked edgels and P is the maximium difference
% between the angles in degree of the gradient vectors of the 2 linked 
% edgels. L is the minimum length of a chain to be included.
% Default L = 3, T = 60 and P = 30.
%
% C is a ?x4 with the following format
% [max_idx  no_chains   no_edgels      0     ] <---header of structure
% [length  closed_flag    xmax       ymax    ] <---header of each chain
% [  x         y        magitude  orientation]
% [  x         y        magitude  orientation] <---body of each chain
% [  x         y        magitude  orientation]
%    .         .          .           .
%    .         .          .           .  
%    .         .          .           .
% [length  closed_flag    xmax       ymax    ] <---header of each chain
% [  x         y        magitude  orientation]
% [  x         y        magitude  orientation] <---body of each chain
% [  x         y        magitude  orientation]
%    .         .          .           .
%    .         .          .           .  
%    .         .          .           .
%
% max_idx is the index to row of E where the longest chain starts.
% It starts counting from the 2nd row of E and so if the first chain
% is the longest chain, max_idx would be 1 instead of 2.
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
% The matrix E is the output of canny.m and
% has the following format:
%   E(:,:,1) = edge image : 1 = weak edgel, 
%                           2 = connected weak edgel and
%                           3 = strong edgel
%   E(:,:,2) = sub-pixel x-coordinates
%   E(:,:,3) = sub-pixel y-coordinates
%   E(:,:,4) = gradients
%   E(:,:,5) = orientations
%
% see kenDist.
%

% written by Kenneth Wong.
% last modified on 12-Feb-2000.

switch(nargin)
case 1
  L = 3;
  T = 60;
  P = 30;
case 2
  T = 60;
  P = 30;
case 3
  P = 30;  
end

no_edgels = 0;
no_chains = 0;
max_idx   = 0;
max_len   = 0;
C = [];

[sy, sx, layers] = size(E);

I = E(:,:,1);
X = E(:,:,2);
Y = E(:,:,3);
M = E(:,:,4);
O = E(:,:,5);

idx = find(I>1);

NEXT = zeros(sy,sx);
PREV = zeros(sy,sx);
FLAG = zeros(sy,sx);

sin_T = sin((90-T)*pi/180);
o_P = P*pi/180;
pi2 = 2*pi;

neighbour = [idx+sy idx++sy+1 idx+1 idx-sy+1 idx-sy idx-sy-1 idx-1 idx+sy-1];

for i = 1:length(idx)
  p = idx(i);
  direc1 = O(p);
  pt1 = [X(p); Y(p)];
  
  % the 8-neighbour
  n8 = neighbour(i,:);
  n8 = n8(I(n8)>1);
  
  links = [];
  for j = 1:length(n8)
    
    n = n8(j);
    direc2 = O(n);
    pt2 = [X(n); Y(n)];
    
    % calculate the orientation of the vector from pt1 to pt2
    vec = [pt2-pt1;0]; vec = vec/norm(vec);
    grad = [cos(direc1); sin(direc1); 0];
    
    sin_t = cross(grad,vec);
    
    % calculate the difference between the orienatations
    % of the gradient vectors at pt1 and pt2
    o_diff = abs(direc1-direc2);
    if o_diff > pi
      o_diff = pi2 - o_diff;
    end
    
    % the linked edgel must be in the direction 90-T deg to 90+T deg
    % from the gradient vector clock-wisely, and the angle between
    % the gradient vectors must be greater than P degree.
    if (sin_t(3) > sin_T & ...
        o_diff < o_P & ...
      PREV(n) == 0)
      links = [links; n];
    end
  end %% for j = 1:length(n8) 
  
  % find the closest possible edgel
  if ~isempty(links)
    [d, n] = min(kenDist([X(links) Y(links)],[X(p) Y(p)]));
    n = links(n);
    PREV(n) = p;
    NEXT(p) = n;
  end
  
end %% for i = 1:length(idx)

% check for open chains
for i = 1:length(idx)
  p = idx(i);
  
  if (PREV(p)==0 & NEXT(p)>0)
    n = p;
    oc = [];
    
    while(n>0)
      oc = [oc;[X(n) Y(n) M(n) O(n)]];
      FLAG(n) = 1;
      n = NEXT(n);
    end
    
    len = size(oc,1);
    if len >= L
      C = [C; [len 0 sx sy]; oc];
      no_chains = no_chains+1;
      no_edgels = no_edgels+len;
      
      if len > max_len
        max_len = len;
        max_idx = size(C,1)-len;
      end  
    end
    
  end %% if (PREV(p)==0 & NEXT(p)>0)
end %% for i = 1:length(idx)

% check for closed chains
for i = 1:length(idx)
  p = idx(i);
  
  if (FLAG(p)==0 & PREV(p)>0 & NEXT(p)>0)
    cc = [X(p) Y(p) M(p) O(p)];
    FLAG(p) = 1;
    n = NEXT(p);
    
    while(n~=p)
      cc = [cc;[X(n) Y(n) M(n) O(n)]];
      FLAG(n) = 1;
      n = NEXT(n);
    end
    
    len = size(cc,1);
    if len >= L
      C = [C; [len 1 sx sy]; cc];
      no_chains = no_chains+1;
      no_edgels = no_edgels+len;
      
      if len > max_len
        max_len = len;
        max_idx = size(C,1)-len;
      end  
    end
    
  end %% if (FLAG(p)==0 & PREV(p)>0 & NEXT(p)>0)
end %% for i = 1:length(idx)

C = [[max_idx no_chains no_edgels 0]; C];

return
