function cleanKeyNodes(dirOut,fileShp,fileFDR,fileDEM,fileFAC,thErr)
%% Algorithm description: 
% 0) Keep only one node in each grid by selecting the one with the largest
% area in PW.
% 1) Construct the tree structure.
% 2) for each node, search for its 5 downstream nodes (unique) and 5
% upstream nodes (always pick the brach with the largest FAC)
    % 3) Fit a equation WL=DEM+b where DEM is the elevation and WL is the water level and b is the only scalar parameter
    % 4) reject the current node if it yields a big error in the previous
    % equation
addpath('../MEX-2.3.0-matlab-2021')
addpath('../river_topology/')
addpath('../geomorphology/common');
addpath('../blockwise')
nBuf=2;
th=1;
GDALLoad();
S=shaperead(fileShp);    
proj=readProj(fileShp);
ID=[S.ID];
X=[S.X];
Y=[S.Y];
areaScene=[S.AreaScene]';
areaRatio=[S.AreaRatio]';
area=areaScene.*areaRatio;
nodes=GaugeNode(ID,X,Y,proj);
[~,rows,cols,geoTransFDR,projFDR]=RasterInfo(fileFDR);
nodes.GetRowCol(projFDR,geoTransFDR);

%% 0)
ind=sub2ind([rows,cols],[nodes.row],[nodes.col]);
ind=ind';
areaMax=accumarray(ind,area,[],@max,[],true);
[indU,~,areaMax]=find(areaMax);
ind2Remove=[];
for iind=1:length(indU)
    indu=indU(iind);
    indGroup=find(ind==indu);
    if length(indGroup)==1
        continue;
    end
    areaI=area(indGroup);
    ind2Remove=[ind2Remove;indGroup(areaI<areaMax(iind))];
end
nodes(ind2Remove)=[];
% retrieve the area of each polygon
% keep only one at each stream pixel
%% 1)
outlets=nodes.constructNetworkBig(fileFDR,500);
[groundElev,fac]=extractHydro(fileDEM,fileFAC,X,Y,5000);
stageMax=[S.StageMax];
stageMin=[S.StageMin];
stageMean=(stageMax+stageMin)/2;
nNodes=length(nodes);
mask2Modify=false(nNodes,1);
for n=1:nNodes
    error=markErr(nodes(n),nBuf,thErr,stageMax,stageMin,stageMean,groundElev,fac,ID);
    mask2Modify(n)=error;
end
STCDValid=[nodes(~mask2Modify).STCD];
maskValid=ismember(ID,STCDValid);
fileShpValid=fullfile(dirOut,'stages_valid.shp');
fileShp2Modify=fullfile(dirOut,'stages_2Modify.shp');
SValid=S(maskValid);
shapewrite(SValid,fileShpValid);
SMoify=S(~maskValid);
shapewrite(SMoify,fileShp2Modify);
prj=readProj(fileShp);
writeProj(fileShpValid,prj);
writeProj(fileShp2Modify,prj)
end
function error=markErr(node,nBuf,th,maxStageAll,minStageAll,meanStageAll,elevGroundAll,facAll,STCDAll)
error=false;
neighbors=[];
nodeCur=node;
STCD=[];
for i=1:nBuf% go downstream
    if ~isempty(nodeCur.Parent)
        parent=nodeCur.Parent;
        neighbors{end+1}=parent;
        nodeCur=parent;
        STCD=[STCD;parent.STCD];
    else
        break;
    end
end
neighbors{end+1}=node;
indCenter=length(neighbors);
STCD=[STCD;node.STCD];
for i=1:nBuf% go upstream (follows the max branch)
    if isempty(nodeCur.Children)
        break;
    end
    children=nodeCur.Children;
    STCDc=[children.STCD];
    [mask,loc]=ismember(STCDc,STCDAll);
    facC=facAll(loc(mask));
    [~,maxChild]=max(facC);
    child=children(maxChild);
    neighbors{end+1}=child;
    STCD=[STCD;child.STCD];
end
if length(neighbors)<2*nBuf-1
    return;
end
nBuf=length(STCD);
% STCD=[neighbors.STCD];
[mask,loc]=ismember(STCD,STCDAll);
meanStage=meanStageAll(loc(mask));
meanStage=meanStage(:);
maxStage=maxStageAll(loc(mask));
minStage=minStageAll(loc(mask));
elevGround=elevGroundAll(loc(mask));
elevGround=elevGround(:);
maxStageCenter=maxStage(indCenter);
minStageCenter=minStage(indCenter);
X=[ones(nBuf,1),elevGround];
b=regress(meanStage,X);
estMeanCenter=[1,elevGround(indCenter)]*b;

if minStageCenter-estMeanCenter>th || estMeanCenter-maxStageCenter> th
    error=true;
end
end
function [elev,fac]=extractHydro(fileDEM,fileFAC,X,Y,blockSize)
[nR,nC,r0,c0,nRows,nCols]=GetBlock(fileDEM,blockSize,0,0,false);
[~,~,~,geoTrans]=RasterInfo(fileDEM);
[row,col]=Proj2RowCol(geoTrans,Y,X);
elev=nan(length(X),1);
fac=nan(length(X),1);
for iR=1:nR
    for iC=1:nC
        mask=row>=r0(iR) & row<(r0(iR)+nRows(iR)) & col>=c0(iC) & col<(c0(iC)+nCols(iC));
        rowIn=row(mask);
        colIn=col(mask);
        prog=((iR-1)*nC+iC)/(nR*nC);
        fprintf(1,'%d%%...',int32(prog*100));
        if isempty(rowIn)
            continue;
        end
        DEM=ReadRaster(fileDEM,r0(iR),c0(iC),nRows(iR),nCols(iC));
        FAC=ReadRaster(fileFAC,r0(iR),c0(iC),nRows(iR),nCols(iC));
        ind=sub2ind([nRows(iR),nCols(iC)],rowIn-r0(iR),colIn-c0(iC));
        elev(mask)=DEM(ind);
        fac(mask)=FAC(ind);
        if prog<0.1
            fprintf(1,'\b\b\b\b\b');
        elseif prog<1
            fprintf(1,'\b\b\b\b\b\b');
        else
            fprintf(1,'\b\b\b\b\b\b\b');
        end
    end
end
end
