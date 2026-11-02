function outlets=findOutlets(this,FDR)
%% find all outlet gauges without any parent(downstream gauges)
dR=[0  1 1 1 0 -1 -1 -1];
dC=[+1 1  0  -1 -1 -1 0 1];
if length(this)==1
    outlets=this;
    warning('Please input an array of gauges. Otherwise, there is not point of finding');
    return;
end
if length([this.row])<length(this)
    error('run GetRowCol before using this function')
end
isOutlets=true(length(this),1);
rowAll=[this.row];
colAll=[this.col];
for i=1:length(this)
    r0=this(i).row;c0=this(i).col;
    while FDR(r0,c0)>0
        dirCurrent=FDR(r0,c0);
        dirCurrent=log(dirCurrent)/log(2);
        r0=r0+dR(dirCurrent+1);
        c0=c0+dC(dirCurrent+1);
        hasAnotherGauge=any(rowAll==r0 & colAll==c0);
        if hasAnotherGauge
            isOutlets(i)=false;
            break;
        end
    end
end
outlets=this(isOutlets);
end