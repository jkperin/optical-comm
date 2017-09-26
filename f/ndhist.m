function qHist = ndhist(varargin)
%function qHist = ndhist(mData, vEdge1, vEdge2, ..., vEdgen)
%Finds n-dimensional histogram data cube
%**Arguments**
%   mData   : [# Data points, # dimension] table of data points that you want to calcuate its histogram
%   vEdge*  : Edge for each dimension
%**Return value**
%   qHist   : [length(vEdge1)-1,length(vEdge2)-1,...length(vEdgen)-1]

%{
if nargin<1, error('Not enough input arguments.'); end
vSizeData = size(varargin{1});

nPoints = vSizeData(1); 
nDim    = vSizeData(2);
clear vSizeData;

if (nargin<=nDim)
    error('Not enough input arguments.');
end %if (nargin==2)

%%%%%%%%%%%%%%%%%%%%%%%%
% Edge -> dimension of the histogram & address calculation multiplier
%%%%%%%%%%%%%%%%%%%%%%%%
%Given all the edges are given

%Allocation: Length of edge for each dimension - 1 == size of histogram in each dimension
vNEdges_Minus_1 = ones(nDim,1);

%Allocation: Address calculation multiplier. Just copies.
vnMul = vNEdges_Minus_1;

%Address calculation multiplier
%Matlab matrix
% [11 12]
% [21 22];
%Matlab matrix in memory
% [11]  1   = 1 + 0 x 2
% [21]  2   = 2 + 0 x 2
% [12]  3   = 1 + 1 x 2
% [22]  4   = 2 + 1 x 2
%                     ^ multiplier for dim 2 == size of dim 1

%Matlab matrix
% [11 12]
% [21 22]
% [31 32];
%Matlab matrix in memory
% [11]  1   = 1 x 1 + 0 x 3
% [21]  2   = 2 x 1 + 0 x 3
% [31]  3   = 3 x 1 + 0 x 3
% [12]  4   = 1 x 1 + 1 x 3
% [22]  5   = 2 x 1 + 1 x 3
% [32]  6   = 3 x 1 + 1 x 3
%                 ^       ^ multiplier for dim 2 == size of dim 1
% here, vnMul = [1 3]';

%%%%%%%%%%%%%%%%%%%%%%%%
%Dimension loop - build vEdges
%%%%%%%%%%%%%%%%%%%%%%%%

iDim = 1;
vNEdges_Minus_1(iDim) = length(varargin{iDim+1})-1; %size of dim 1

for iDim = 2:nDim
    vNEdges_Minus_1(iDim) = length(varargin{iDim+1})-1; %size of iDim
    
    vnMul(iDim) = vNEdges_Minus_1(iDim-1) * vnMul(iDim-1);
    %multiplier for iDim = Pi(size(qHist,i),i=1..(iDim-1))
    
end %for iDim = 1:nDim

pack;
qHist = zeros(vNEdges_Minus_1');

%%%%%%%%%%%%%%%%%%%%%%%%
% Point Loop
%%%%%%%%%%%%%%%%%%%%%%%%
for iPoint = 1:nPoints
    %extract one point from the data
    vThis = varargin{1}(iPoint,:);
    
    %Buffer for index
    vIndexThis = zeros(1,nDim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %Dimension loop
    %%%%%%%%%%%%%%%%%%%%%%%%
    for iDim = 1:nDim
        %Index for this dimension
        
        iIndexThis = FindIndex(varargin{iDim+1}, ...
            vNEdges_Minus_1(iDim), ...
            vThis(iDim));
        
        if (iIndexThis>vNEdges_Minus_1(iDim))
            error('ndhist:index out of range');
        end
        
        %if ~exist(Index for this dimension)
        if isempty (iIndexThis)
            vIndexThis = [];
            break;
        end
        vIndexThis(iDim) = iIndexThis;
    end
    
    if ~isempty(vIndexThis)
        iHere = (vIndexThis-1) * vnMul + 1; %vnMul = zeros(nDim,1)
        qHist(iHere) = qHist(iHere) + 1;
    end
    
end %for iPoint = 1:nPoints

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Histogram matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allocation
%   vNEdges = [(nEdge-1)'s]
%   qHist = zeros(vNEdges);
%Addressing
%   iHere = i1 + (nEdge1Minus1) * i2 + (nEdge1Minus1) * (nEdge1Minus2) * i3 + ... 
%         = [-i*-]*[-nEdge*Minus1-]';
%   qHist(iHere) = qHist(iHere) + 1;
%Histogram Example for vThis
%   vIndexThis = zeros(1,nDim);
%   for iDim = 1:nDim
%       vIndexThis(iDim) = FindIndex(vEdges(vpEdgeStart(iDim), vpEdgeEnd(iDim)), ...
%                           vNEdges_Minus_1(iDim), ...
%                           vThis(iDim))
%   end
%   iHere = vIndexThis * vnMul; %vnMul = zeros(nDim,1)
%   qHist(iHere) = qHist(iHere) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Bin edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%One long vector including all edges 
%   vEdges(sum(length(edges)),1)
%Another vector indicating starting point of each edges
%   vpEdgeStart = zeros(nDim,1);
%   vpEdgeEnd   = zeros(nDim,1);
%Addressing
%   vEdges(vpEdges(iDim)+iEdge)
%Extracting one edge
%   vEdges(vpEdgeStart(iDim):vpEdgeEnd(iDim))
%Query
%   max(find(vEdges(vpEdgeStart(iDim):vpEdgeEnd(iDim)) < vThis(iDim)));


%}function qHist = ndhist(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iEdge = FindIndex(vEdge, nUB, rData)
%{
%   nUB == length(vEdge) - 1;
iEdge = min(max(find(vEdge < rData)),nUB);
%}