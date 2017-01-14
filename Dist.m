function W = Dist( line, Data )
%DIST Calculates distances to the line

W = line * [Data; ones([1 size(Data, 2)])];