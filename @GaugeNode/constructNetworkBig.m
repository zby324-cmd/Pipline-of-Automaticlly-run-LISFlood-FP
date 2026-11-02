function [outlets,lst2Remove]=constructNetworkBig(lstGauges,fileFDR,bufSize,thErr)
%% This function assumes the FDR file is bigger than the memory can afford
% The algorithm is to find the parent of each station then link
n=length(lstGauges);
[~,rows,cols]=RasterInfo(fileFDR);
rowAll=[lstGauges.row];
colAll=[lstGauges.col];
outlets=[];
lst2Remove=false(n,1);
for i=1:n
    rowC=lstGauges(i).row;
    colC=lstGauges(i).col;
    foundPar=false;
    fprintf('searching parent for node %d...',lstGauges(i).STCD);
    bCannotConn=false;
    while ~foundPar && ~bCannotConn
        row0=max(1,rowC-bufSize);
        col0=max(1,colC-bufSize);
        row1=min(rows-row0+1,2*bufSize+1);
        col1=min(cols-col0+1,2*bufSize+1);
        FDR=ReadRaster(fileFDR,row0,col0,row1,col1);
        rowCL=rowC-row0+1;
        colCL=colC-col0+1;
        if rowCL<1 || rowCL>row1 || colCL<1 || colCL>col1
            disp('incomplete hydrogrpahy')
            break;
        end
        % h=figure;
        % imagesc(FDR);
        % savefig(h,'/raid-13/SFS/xinyis/retrieval/flood_depth/ohio/new_polygons3/FDR.fig')
        iMove=0;
        while ~foundPar
            [nextRow,nextCol]=GetNextCellV(FDR(rowCL,colCL),rowCL,colCL);
            iMove=iMove+1;
            if iMove>thErr
                bCannotConn=true;
                disp('erroneous hydrography')
                lst2Remove(i)=true;
                break;
            end
            if nextRow<1 ||nextRow>row1 || nextCol<1 || nextCol>col1
                rowC=row0+nextRow-1;
                colC=col0+nextCol-1;
                % fprintf('centered to a new block');
                break;
            end
            if nextRow==rowCL && nextCol==colCL% cannot find the next downstream cell. It is an outlet
                foundPar=true;
                outlets=[outlets;lstGauges(i)];
                fprintf('the node is an outlet.\n');
                break;
            end
            mask=(rowAll-row0+1)==nextRow & (colAll-col0+1)==nextCol;
            if sum(mask)==0 || sum(mask)>1
                rowCL=nextRow;
                colCL=nextCol;
            else
                nodePar=lstGauges(mask);
                lstGauges(i).Parent=nodePar;
                if isempty(nodePar.children)
                    nodePar.children=lstGauges(i);
                else
                    nodePar.children(end+1)=lstGauges(i);
                end
                fprintf('found.\n');
                foundPar=true;
                break;
            end
        end
    end
end
end
