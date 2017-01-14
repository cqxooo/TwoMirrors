function sp = loadBSplines(F)
  
%  sp = loadTxtSpl(F) return a N+1x2 matrix containing
%  the control points of a B-Spline snake. F is the
%  filename of the text file containing the snake information.
%  sp(1,1) gives the number of control points and
%  sp(1,2) tells where the snake is closed(1) or open(0).
%  sp(?,1) and sp(?,2) gives the x and y coordinates of the
%  control points. This function assumes that the text file only
%  contains a single snake.

sp = [0 0];

fid = fopen(F, 'r');
if (fid==-1)
  return
end

header = fscanf(fid,'%s',1);
n_snakes = fscanf(fid,'%d\n',1);

if (n_snakes<1)
  warning('The text file has less than 1 B-Spline snake!');
  return
end

%disp(sprintf('%d B-Spline snake(s) in %s', n_snakes, F));

closed = fscanf(fid,'%d',1);
n_cpoints =  fscanf(fid,'%d\n',1);

sp_head = [n_cpoints closed];

sp = fscanf(fid,'%f',[2 n_cpoints])';
sp = sp(:,[1 2]);

sp = [sp_head; sp];

fclose(fid);

return
