function constructNetwork(this,FDR,lstGauges)
%% construct river topology
%% input 
% FDR: the flow direction map
% lstGauges: last downstream gauges (used in recursive calls)
idR=[0  -1 -1 -1 0 1 1 1];
idC=[-1 -1  0  1 1 1 0 -1];
% initialize the upstream pixel by this current pixel
uR=this.row;uC=this.col; 
while ~isempty(uR)
    % pop an upstream pixel
    uRow=uR(1);
    uCol=uC(1);
    uR(1)=[];
    uC(1)=[];
    for ifdr=1:GaugeNode.nFDR
        if uRow+idR(ifdr)<1 || uRow+idR(ifdr)>size(FDR,1) || uCol+idC(ifdr)<0 ||uCol+idC(ifdr)>size(FDR,2)
            continue;
        end
        if FDR(uRow+idR(ifdr),uCol+idC(ifdr))==GaugeNode.fdrCode(ifdr)% the pixel is a new upstream pixel
            uRN=uRow+idR(ifdr);
            uCN=uCol+idC(ifdr);
            hasGauge=false;
            for iG=1:length(lstGauges)
                if lstGauges(iG).row==uRN && lstGauges(iG).col==uCN
                    if isempty(this.children)
                        this.children=lstGauges(iG);
                    else
                        this.children(end+1)=lstGauges(iG);
                    end
                    this.children(end).parent=this;
                    disp([this.children(end).STCD,' is connected']);
                    hasGauge=true;
                    break;
                end
            end
            if ~hasGauge% if no gauge is found, push this new upstream pixel to the upstream pixel list
                uR(end+1)=uRN;
                uC(end+1)=uCN;
            else % a gauge is found at this new upstream pixel, add it as a child and call make it construct its own network
                this.children(end).constructNetwork(FDR,lstGauges);
            end
        end
    end
end
end