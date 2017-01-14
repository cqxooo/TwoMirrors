%RANSAC Implements base algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    [vMask, Model] = RANSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold )
%    -----------------------------------------------------------------------------------
%    Arguments:
%           mData - matrix of data, where each column-vector is point
%           ModelFunc - handle to Model Creating function. It must create a
%                   model from nSampLen column-vectors organized in
%                   matrix
%           nSampLen - number of point for ModelFunc
%           ResidFunc - handle to Residuum calculating function. As
%                   argument this function takes model, calculated by
%                   ModelFunc, and matrix of data (all or maybe part of it)
%           nIter - number of iterations for RANSAC algorithm
%           dThreshold - threshold for residuum
%    Return:
%           vMask - 1s set for inliers, and 0s for outliers
%           Model - approximate model for this data


function [vMask, Model] = RANSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold )

% Cheking arguments
if length(size(mData)) ~=2
    error('Data must be organized in column-vecotors massive');
end

nDataLen = size(mData, 2);

% Initialization
Model = NaN;
vMask = zeros([1 nDataLen]);
nMaxInlyerCount = -1;

% Main cycle
for i = 1:nIter
    
    % 1. Sampling
    SampleMask = zeros([1 nDataLen]);
    
    % Takes nSampleLen different points
    while sum( SampleMask ) ~= nSampLen
        ind = ceil(nDataLen .* rand(1, nSampLen - sum(SampleMask)));  
        SampleMask(ind) = 1;
    end    
    Sample = find( SampleMask );
    
    
    % 2. Creating model
    ModelSet = feval(ModelFunc, mData(:, Sample));
    
    for iModel = 1:size(ModelSet, 3)
        
        CurModel = ModelSet(:, :, iModel);
        
        % 3. Model estimation
        CurMask = ( abs( feval(ResidFunc, CurModel, mData)) < dThreshold );
        nCurInlyerCount = sum(CurMask);
    
        % 4. The best is selected
        if nCurInlyerCount > nMaxInlyerCount

            % Save some parameters
            nMaxInlyerCount = nCurInlyerCount;
            vMask = CurMask;
            Model = CurModel;
        end
    end
end

return; 
%END of RANSAC