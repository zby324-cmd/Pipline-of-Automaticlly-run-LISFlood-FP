function sr_e_ind=findGAREvt(sr,er,sh_e,eh_e,LSP,Q,prec)
%% find the start of group A rainfall event of the flood event e
%% input
% sr: the start of all rainfall events
% er: the end of all rainfall events
% sh_e: the start of the given flood event 
% LSP: empirical lag between the start of rain fall and flood event
% Q: flow in mm/h, which is normalized by basin area to match the unit of precipitation
%% sr_e matched index of the start of rain fall event
%% mainbody
% flow from the start of the event to the end of earliest possible rainfall start  
%% The time line of all precipitation events within the range
sr_e_cand_ind=find(sh_e-sr>=0 & sh_e-sr<min([sh_e-eh_e,3*LSP]));
sr_e_cand=sr(sr_e_cand_ind);
er_e_cand=er(sr_e_cand_ind);
if isempty(sr_e_cand) || any(isnan(Q(sr_e_cand(1):sh_e)))%this function requires the flow data to be complete
    sr_e_ind=[];
%     disp('missing precipiation or flow data');
    return;
end
np=length(sr_e_cand);
p_accum_cand=zeros(np,1);
q_accum_cand=zeros(np,1);
%% flow accumulation from the start time of every precipitation
% and precipiation event volume accumulation
for ip=1:np
    p_accum_cand(ip)=sum(prec(sr_e_cand(ip):er_e_cand(ip)));
    q_accum_cand(ip)=sum(Q(sr_e_cand(ip):sh_e));
end
r_cand=q_accum_cand./p_accum_cand;
ind1=find(r_cand<=0.8,1,'first');
if isempty(ind1)
    sr_e_ind=[];
    return;
end
indStart=find(sr==sr_e_cand(ind1));
indEnd=sr_e_cand_ind(end);
sr_e_ind=(indStart:indEnd)';
end