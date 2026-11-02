function outlet=genRiverNetwork(shapefile,STCDOut,fdrFile,demFile,dirQ)
GDALLoad();
NoDataVal=-9999;
[xLoc,yLoc,proj,IDs]=readShapeLoc(shapefile,0);
%% construct gauges
nodes=GaugeNode(IDs,xLoc,yLoc,proj);
%% remove gagues without records
for i=length(nodes):-1:1
    nodes(i).ImportFlow(dirQ,NoDataVal);
    if ~nodes(i).HasObs
        nodes(i)=[];
    end
end
[fdr,geoTrans,projFDR]=ReadRaster(fdrFile);
nodes.GetRowCol(projFDR,geoTrans);
indexOut=strcmpi({nodes.STCD},STCDOut);
outlet=nodes(indexOut);
nodes(indexOut)=[];
outlet.constructNetwork(fdr,nodes)
STCDOI='01172010';
[outlet,goi,arrParaGauges,arrSubGauges]=outlet.GetAdjGauges(STCDOI);
beta=estimate_adj_coef(outlet,[arrParaGauges,arrSubGauges],goi,fdrFile,demFile,...
    'timeUnit','d','coeM',88.4334,'expM',0.7139,'coeR',1.9048,'GM',1.364364,'STCDSubOI',arrSubGauges);
end