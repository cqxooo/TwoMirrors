function getOverallEdge(pbv, row, col, filename)
data = zeros(row,col);
for i =1: length(pbv)
   X = round(pbv{i})';
   tdata = zeros(row,col);
   for j = 1:size(X,1)
        tdata(X(j,2),X(j,1)) = 1;     
   end
   for j = 1:row
        Idx = find(tdata(j,:));
        if size(Idx,2)~=0
            for k = Idx(1):Idx(end)
                tdata(j,k) = 1 ;   
            end
        end
   end
   data = data+tdata;
end
imData = zeros(row,col,3);
%     imData(:,:,1)=imData(:,:,1)+100; 
%     imData(:,:,2)=imData(:,:,2)+100;    
%     imData(:,:,3)=imData(:,:,3)+255;
    for j=1:row
        for k=1:col
            if data(j,k)~=0
                imData(j,k,1)=255;
                imData(j,k,2)=255;
                imData(j,k,3)=255;                
            end
        end
    end
imData =uint8(imData);
imwrite(imData,filename);
