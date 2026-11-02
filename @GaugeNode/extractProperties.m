function [rows,cols,hDistTotal,vDistTotal,values]=extractProperties(node, fieldNameExist,fieldName)
%% This function recursively until all branches meet a filled node
if node.(fieldNameExist)
    values=nan(1,length(fieldName));
    for ifield=1:length(fieldName)
        values(ifield)=node.(fieldName{ifield});
    end
    rows=node.row;
    cols=node.col;
    hDistTotal=sum(node.hDist);
    vDistTotal=sum(node.vDist);
else
    rows=[];
    cols=[];
    values=[];
    hDistTotal=[];
    vDistTotal=[];
    hDist2Par=sum(node.hDist);
    vDist2Par=sum(node.vDist);
    for c=1:length(node.children)
        [rowsC,colsC,hDc,vDc,valuesC]=extractProperties(node.children(c),... 
            fieldNameExist,fieldName);
        values=[values;valuesC];
        rows=[rows;rowsC];
        cols=[cols;colsC];
        hDistTotal=[hDistTotal;hDc+hDist2Par];
        vDistTotal=[vDistTotal;vDc+vDist2Par];
    end
end
end