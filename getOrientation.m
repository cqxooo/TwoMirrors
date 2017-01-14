function RR = getOrientation(polyBoundaryVecCell, K, lss, lhs)
RR={};
for imLoop=1:length(polyBoundaryVecCell)
    RR{imLoop} = rectify(polyBoundaryVecCell{imLoop}, K, lss(:,imLoop), lhs(:,imLoop));
end
