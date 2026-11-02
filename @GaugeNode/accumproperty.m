function val=accumproperty(nodeDown,nodeUp,propertyName,funcHandle)
bFound=false;
nodeCur=nodeUp;
val=0;
while ~bFound
    parent=nodeUp.Parent;
    if parent.STCD==nodeDown.STCD)
        bFound=true;
    end
    val=val+funcHandle(nodeCur.(propertyName));
    nodeCur=parent;
end
end