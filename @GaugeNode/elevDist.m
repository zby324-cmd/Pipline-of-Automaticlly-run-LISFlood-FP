function elevDist(this,allBetween,fileDEM)
%% This function computes the horizontal distance between the current node and its parent.
% The path2parent field has to be filled before calling this function
%% input parameters
% allBetweeh: true|false: output results bewteen each pixel (including the
% current node and the parent| only compute the results between the current
% node and its parent
% geoTrans: The 6-GT coefficient defining the CS
% isGeographic: if the projection is in length|lat,lon
%%
if isempty(this.parent)
    disp('skip root');
    return;
end
if allBetween
    if ~isempty(this.Path2Parent)
        rows=[this.row;this.Path2Parent(:,1);this.Parent.row];
        cols=[this.col;this.Path2Parent(:,2);this.Parent.col];
    else
        rows=[this.row;this.Parent.row];
        cols=[this.col;this.Parent.col];
    end
else
    error('not implemented')
end
maxRow=max(rows);
r0=min(rows);
maxCol=max(cols);
c0=min(cols);
rs=maxRow-r0+1;
cs=maxCol-c0+1;
DEM=ReadRaster(fileDEM,r0,c0,rs,cs);
ind=sub2ind([rs,cs],rows-r0+1,cols-c0+1);
elev=DEM(ind);
this.vDist=elev(1:end-1)-elev(2:end);
end