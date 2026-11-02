function [vos,maxID]=confluence(this,fileFDR,maxID,blockSize)
%% Algorithm Description: 
% This function adds virtual nodes (VN) at confluence location. 
% Confluence location is the pixel where multiple branches separates.
% confluence nodes are pivot for interpolation of water level.
% 1) Identify the flow paths from a parent to all children.
% 2) compute the overlapped parts between any two paths and use a table to
% store the results
% 3) For a given child, select all siblings with overlapped paths then find the common overlapping paths. 
% 4) Selected siblings are removed from the table
% 5) Create a VN at the last pixel of the overlappinp path and insert it
% between the parent and these selected children.
% 6) recursively call the function using the new VN as the parent until no VN is needed. 
vos=[];
children=this.Children;
nChildren=length(children);
commonPaths=cell(nChildren,nChildren);
paths=cell(nChildren,1);
%% step 1)
[~,rows,cols,geoTrans,proj]=RasterInfo(fileFDR);
for c=1:nChildren
    child=children(c);
    child.findParent(fileFDR,blockSize);
    if ~isempty(child.Path2Parent) % if a child is next to the parent, the path is empty.
        paths{c}=sub2ind([rows,cols],child.Path2Parent(:,1),child.Path2Parent(:,2));
    end
end
%% step 2)
for c1=1:nChildren
    path1=paths{c1};
    for c2=c1+1:nChildren
        path2=paths{c2};
        isCommon=ismember(path1,path2);
        commonPaths{c1,c2}=path1(isCommon);
        commonPaths{c2,c1}=commonPaths{c1,c2};
    end
end
%% step 3)
grouped=false(nChildren,1);
allConfNodes=[];
for c1=1:nChildren
    if grouped(c1)
        continue;
    end
    lenCPath=Inf;
    confNodes=c1;
    for c2=1:nChildren
        if c1==c2 || grouped(c2)
            continue;
        end
        %% 4)
        if ~isempty(commonPaths{c1,c2})
            if lenCPath>length(commonPaths{c1,c2})
                lenCPath=length(commonPaths{c1,c2});
                cPath=commonPaths{c1,c2};
            end
            confNodes=[confNodes;c2];
            grouped(c2)=true;
        end
    end
    if ~isinf(lenCPath)
        allConfNodes=[allConfNodes;confNodes];
        %% 5)
        [rowVN,colVN]=ind2sub([rows,cols],cPath(1));
        [Y,X]=RowCol2Proj(geoTrans,rowVN,colVN);
        vnc=GaugeNode(maxID+1,X,Y,proj);
        maxID=maxID+1;
        vnc.Parent=this;
        [row2Par,col2Par]=ind2sub([rows,cols],cPath(2:end));
        vnc.Path2Parent=[row2Par,col2Par];
        lenCPath=length(cPath);
        vnc.GetRowCol(proj,geoTrans);
        children(c1).Path2Parent=children(c1).Path2Parent(1:end-lenCPath);
        vnc.children=children(c1);
        children(c1).Parent=vnc;
        for c2=2:length(confNodes)
            children(confNodes(c2)).Path2Parent=children(confNodes(c2)).Path2Parent(1:end-lenCPath);
            vnc.children(end+1)=children(confNodes(c2));
            children(confNodes(c2)).Parent=vnc;
        end
        this.children(end+1)=vnc;% append the virtual node to the last child
        %% 6)
        [vosC,maxID]=vnc.confluence(fileFDR,maxID,blockSize);
        vos=[vos;vnc;vosC];
    end
    grouped(c1)=true;
end
if ~isempty(allConfNodes)
    this.children(allConfNodes)=[];
end
end