function RR = getOrientation3(polyBoundaryVecCell, K, cps)
for imLoop=1:length(polyBoundaryVecCell)    
    i = cps(:,2*imLoop-1);
    j = cps(:,2*imLoop);
    r2 = inv(K)*imag(i);
    r3 = cross(inv(K)*real(i), inv(K)*imag(i));
    r1 = inv(K)*real(i);
    RR{imLoop} = [r1/norm(r1) r2/norm(r2) r3/norm(r3)];    
end