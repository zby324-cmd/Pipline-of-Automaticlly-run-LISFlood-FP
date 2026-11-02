function ImportFlow2(this,dirQ,fmt,NoDataVal,dt)
QfileSub=[dirQ,this.STCD,'_obs.csv'];
 this.dt=dt;
 if exist(QfileSub,'file')~=2
     disp(['observed file of site ', num2str(this.STCD),'does not exist. Site is removed']);
     return;
 end
 T=readtable(QfileSub);
 this.dates=datenum(num2str(T.Date),fmt);
 this.Qobs=T.Discharge;
 if exist('NoDataVal','var')==1
     this.Qobs(this.Qobs==NoDataVal |this.Qobs<0)=NaN;
 end
end