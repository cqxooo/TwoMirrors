function writePrjMat(F, PrjMat)
fid = fopen(F,'w');
rows = size(PrjMat,1);
for i = 1:rows;
    fprintf(fid, '%10.14e %10.14e %10.14e %10.14e \n', PrjMat(i,:));
end
fprintf(fid, '\n');
fclose(fid);

return;