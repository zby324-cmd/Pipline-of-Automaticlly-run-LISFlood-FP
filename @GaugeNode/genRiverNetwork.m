function outlet=genRiverNetwork(shapefile,STCDOut,fdrFile,demFile,dirQ,dirMask,dirOut)
addpath('E:\Dropbox\code\MFiles\commonHydro\')
GDALLoad();
NoDataVal=-9999;
[xLoc,yLoc,proj,IDs]=readShapeLoc(shapefile,0);
%% construct gauges
nodes=GaugeNode(IDs,xLoc,yLoc,proj);
%% remove gagues without records
for i=length(nodes):-1:1
    nodes(i).ImportFlow(dirQ,NoDataVal);
    fileMask=[dirMask,nodes(i).STCD,'_mask.tif'];
    nodes(i).SetArea(fileMask);
    if ~nodes(i).HasObs
        nodes(i)=[];
    end
end
indexOut=strcmpi({nodes.STCD},STCDOut);
outlet=nodes(indexOut);
nodes(indexOut)=[];
[fdr,geoTrans,projFDR]=ReadRaster(fdrFile);
nodes.GetRowCol(projFDR,geoTrans);
outlet.GetRowCol(projFDR,geoTrans);
outlet.constructNetwork(fdr,nodes);
% global outletOI goi arrParaGauges arrSubGauges
shapeOut=[dirOut,'out.shp'];
delete(shapeOut)
CreateShape(shapeOut,'STCD',str2OTFType('OFTString'), 'errSub',str2OTFType('OFTReal'),'errSum',str2OTFType('OFTReal'));
for i=1:length(nodes)
    STCDOI=nodes(i).STCD;
%     STCDOI='01144000';
%     STCDOI='01161000';
    figureOut=[dirOut,STCDOI,'.fig'];
    disp(['processing', STCDOI ,'...']);
    [outletOI,goi,arrParaGauges,arrSubGauges]=outlet.GetAdjGauges(STCDOI);
    if isempty(outletOI) || isempty(goi)
        disp('skipped, outlet information is not found.')
        continue;
    end
    [beta,areaRatioSub,areaRatioSum]=estimate_adj_coef(outletOI,[arrParaGauges,arrSubGauges],goi,fdrFile,demFile,figureOut,...
        'timeUnit','d','coeM',88.4334,'expM',0.7139,'coeR',1.9048,'GM',1.364364,'STCDSubOI',arrSubGauges);
    AddAPoint(shapeOut,nodes(i).X, nodes(i).Y,'STCD',STCDOI,'errSub',num2str(areaRatioSub),'errSum',num2str(areaRatioSum))
end
end