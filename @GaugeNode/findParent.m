function findParent(this,fileFDR,blockSize)
%% Algortihm Description
% Find the parent node
%% Input
% iNode
% lstGauges
% fileFDR: the Flow Direction map, which can be big
% the: The buffer size for each blockReading
%% Output
path=[];
if isempty(this.Parent)
    disp([num2str(this.STCD) ' is an outlet']);
    return;
end
[~,rows,cols]=RasterInfo(fileFDR);
rowC=this.row;
colC=this.col;

rowPar=this.Parent.row;
colPar=this.Parent.col;
foundPar=false;
while ~foundPar
    row0=max(1,rowC-blockSize);
    col0=max(1,colC-blockSize);
    row1=min(rows-row0+1,2*blockSize+1);
    col1=min(cols-col0+1,2*blockSize+1);
    FDR=ReadRaster(fileFDR,row0,col0,row1,col1);
    rowCL=rowC-row0+1;
    colCL=colC-col0+1;
    if rowCL<1 || rowCL>row1 || colCL<1 || colCL>col1
        disp('incomplete hydrogrpahy')
        break;
    end
    while ~foundPar
        [nextRow,nextCol]=GetNextCellV(FDR(rowCL,colCL),rowCL,colCL);
        if nextRow<1 ||nextRow>row1 || nextCol<1 || nextCol>col1
            rowC=row0+nextRow-1;
            colC=col0+nextCol-1;
            break;
        end
        if rowPar-row0+1~=nextRow || colPar-col0+1~=nextCol
            rowCL=nextRow;
            colCL=nextCol;
            path=[path;nextRow+row0-1,col0+nextCol-1];
        else
            foundPar=true;
            break;
        end
    end
end
this.path2parent=path;
end