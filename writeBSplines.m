function writeBSplines(F,C,D)
% F gives the output BSpline filename
% C gives the control points of the BSpline
% D tells whether the BSpline is closed(1) or open(0)
% This function assumes that the text file only
%  contains a single snake.

fid = fopen(F,'w');
file_type = 'SnakeCtrlPtFile';
n_snakes = 1;
n_cps = size(C,1);
fprintf(fid, '%s ', file_type);
fprintf(fid, '%d\n', n_snakes);
fprintf(fid, '%d %d\n', D, n_cps);

for j = 1:n_cps;
    fprintf(fid, '%10.14f %10.14f\n', C(j, [2 1]));
end

fclose(fid);

return;