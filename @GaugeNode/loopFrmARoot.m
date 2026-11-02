function loopFrmARoot(this,funcHandle,varargin)
%% This function loops over all nodes in a tree from the root to leaf to run user-specified function.
%% usage, a user can call this function only from the root node
% By default, children result can have dependence on the parent.
% this: The handle of a node
% func: The user-specified function handle
funcHandle(this,varargin{:});
for c=1:length(this.children)
    child=this.children(c);
    child.loopFrmARoot(funcHandle,varargin{:});
end
end