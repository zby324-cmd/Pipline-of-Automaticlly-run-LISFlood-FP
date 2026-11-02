function [outlet,goi,arrParaGauges,arrSubGauges]=GetAdjGauges(this,STCDOI)
goi=this.findSub(STCDOI);
if isempty(goi)
    outlet=[];
    arrParaGauges=[];
    arrSubGauges=[];
    return;
end
outlet=goi.parent;
arrParaGauges=outlet.children;
indexSelf=strcmpi({arrParaGauges.STCD},STCDOI);
arrParaGauges(indexSelf)=[];
arrSubGauges=goi.children;
end