function ImportFlow(this,dirQ,fmt,dateCol,simCol,obsCol,NoDataVal)
 QfileSub=[dirQ,this.STCD,'_obs.csv'];
 if exist(QfileSub,'file')~=2
     disp(['observed file of site ', num2str(this.STCD),'does not exist. Site is removed']);
     return;
 end
 T=readtable(QfileSub);
 [rowsRaw,colsRaw]=size(T);
 [rowsData,colsData]=size(data);
 if colsRaw>colsData
     diff=colsRaw-colsData;
     simCol=simCol-diff;
     obsCol=obsCol-diff;
     if rowsData<rowsRaw
         rDStart=2;
     else
         rDStart=1;
     end
     this.dates=datenum(raw(rDStart:end,dateCol),fmt);
 else
     this.dates=datenum(num2str(data(:,dateCol)),fmt);
 end
 if ~isempty(simCol)
     this.Qsim=data(:,simCol);
 end
 if ~isempty(obsCol)
    this.Qobs=data(:,obsCol);
    this.Qobs(this.Qobs==NoDataVal |this.Qobs<0)=NaN;
 end
 
end