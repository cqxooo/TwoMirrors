function RR = getOrientation2(polyBoundaryVecCell, K, lss, lhs, vxs)
x0=K(1:2,3);
for imLoop=1:length(polyBoundaryVecCell)    
    ls = lss(:,imLoop);
    lh = lhs(:,imLoop);
    vx = vxs(:,imLoop);
    lp = [lh(2); -lh(1); -lh(2)*x0(1)+lh(1)*x0(2)];
    vy = cross(lp, ls);
    vz = cross(lh, ls);
    r2 = inv(K)*vy;
    r3 = inv(K)*vz;
    r1 = inv(K)*vx;
    RR{imLoop} = [r1/norm(r1) r2/norm(r2) r3/norm(r3)];
end