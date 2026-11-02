classdef GaugeNode<dynamicprops
    properties(Constant)
        fdrCode=2.^(0:7);
        nFDR=8;
    end
    properties(Access=public)
        STCD;%station code
        X;Y;% (X,Y) in the given coordinate system
%         elevMap;% elevation difference to each outlet
%         distMap;% hydrological distance to each outlet
        row,col;
    end
    properties(SetAccess=private)
       proj;
       parent;
       children;% sub gauges
       path2parent;%(row,col)of the path from the child node to the parent. 
       % The (row,col) of the child and the parent are excluded
       Qsim;
       Qobs;
       dates;
       Area;
       Level;
       dt;
       hDist;
       vDist;
    end
    properties(Dependent)
        HasObs;
        nRec;
        Parent;
        Children;
        Path2Parent;
        HasChildren;
        HasParent;
        HDist;
        VDist;
    end
    methods
        function this=GaugeNode(STCD,X,Y,proj)
            %construction function to construct a gauge
            if nargin==0
                this.STCD=NaN;
                this.X=NaN;
                this.Y=NaN;
                this.proj='';
                return;
            end
            if length(X)>1
                n = length(X);
                this(n,1) = GaugeNode;
                if ~iscell(STCD)
                    STCD=num2cell(STCD);
                end
                [this.STCD]=deal(STCD{:});
                Xc=num2cell(X);
                Yc=num2cell(Y);
                [this.X]=Xc{:};
                [this.Y]=Yc{:};
                [this.proj]=deal(proj);
            elseif length(X)==1
                this.STCD=STCD;
                this.X=X;
                this.Y=Y;
                this.proj=proj;
            end
        end
        
        function GetRowCol(these,projFDR,geoTransFDR)
           [XFDR,YFDR]=ProjTransform(these(1).proj,projFDR,[these.X],[these.Y]);
           [r,c]=Proj2RowCol(geoTransFDR,YFDR,XFDR);
           rc=num2cell(r);
           cc=num2cell(c);
           [these.row]=rc{:};
           [these.col]=cc{:};
        end
        function SetArea(this,value)
            if ischar(value)
                [mask,geoTrans,projMask]=ReadRaster(value);
                bGCS=IsGeographic(projMask,geoTrans);
                this.Area=RasterArea(mask,geoTrans,bGCS);
            elseif isnumeric(value)
               areac=num2cell(value);
               [this.Area]=areac{:}; 
            end
        end
        
        function subNetwork=findSub(this,STCDOI)
            subNetwork=[];
            if strcmpi(STCDOI,this.STCD)
                subNetwork=this;
                return;
            end
            if ~isempty(this.children)
                for i=1:length(this.children)
                    isubNetwork=this.children(i).findSub(STCDOI);
                    if ~isempty(isubNetwork)
                        subNetwork=isubNetwork;
                        return;
                    end
                end
            else
                subNetwork=[];
            end
        end
        function GetLevel(this)
            if isempty(this.children)
                this.Level=1;
            else
                for i=1:length(this.children)
                    this.children(i).GetLevel();
                end
                this.Level=max([this.children.Level])+1;
            end
        end  
        
        function value=get.HasObs(this)
            if isempty(this.Qobs)
                disp('Flow data has not been imported yet. Please call "importFlow" first.');
                value=false;
                return;
            end
            ratio=sum(~isnan(this.Qobs))/((this.dates(end)-this.dates(1))/this.dt);
            if ratio<0.3
                disp([this.STCD,' has a ratio of ',num2str(ratio)])
                value=false;
            else
                value=true;
            end
        end
        function value=get.nRec(this)
            if isempty(this.Qsim)
                disp('Flow data has not been imported yet. Please call "importFlow" first.');
                value=[];
                return;
            end
            value=length(this.Qsim);
        end
        function value=get.Parent(this)
            value=this.parent;
        end
        function set.Parent(this,value)
            this.parent=value;
        end
        function value=get.Children(this)
            value=this.children;
        end
        function value=get.HasChildren(this)
            value=~isempty(this.children);
        end
        function value=get.Path2Parent(this)
            value=this.path2parent;
        end
        function value=get.HasParent(this)
            value=~isempty(this.parent);
        end
        function set.Path2Parent(this,value)
            this.path2parent=value;
        end
        function value=get.HDist(this)
            value=this.hDist;
        end
        function value=get.VDist(this)
            value=this.vDist;
        end
        
        loopFrmARoot(this,funcHandle,varargin);        
        [rows,cols,hDist,vDist,values]=extractProperties(node, fieldNameExist,fieldName);
        [outlets,lst2Remove]=constructNetworkBig(lstGauges,fileFDR,bufSize,thErr);
        [vos,maxID]=confluence(this,fileFDR,maxID,blockSize);
        outlets=findOutlets(this,FDR);
        ImportFlow(this,dirQ,fmt,dateCol,simCol,obsCol,NoDataVal);
        ImportFlow2(this,dirQ,fmt,NoDataVal,dt);
        [outlet,goi,arrParaGauges,arrSubGauges]=GetAdjGauges(this,STCDOI);
        findParent(this,fileFDR,blockSize);
        hydroDist(this,allBetween,geoTrans,isGeographic);
        elevDist(this,allBetween,fileDEM);
    end
    methods(Static)
        val=accumproperty(nodeDown,nodeUp,propertyName,funcHandle);
    end
end