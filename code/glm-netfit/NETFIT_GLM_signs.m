%SCRIPT TO ESTIMATE W FROM SAMPLE OF SPIKE HISTORIES USING DALE PRIOR
fdale=1;                                % enable execution of GLM_dale

dalepr=zeros(N,1);                      % IDENTIFY COLUMN SIGNS
tmp=sum(Weights_glm,3);
for k=1:N
  ind=tmp(:,k)>0; 
  S2pl=sum(tmp(ind,k).^2);              % sum-square for positive weights
  ind=tmp(:,k)<0; 
  S2mi=sum(tmp(ind,k).^2);              % sum-square for negative wieghts
  
  dalepr(k)=sign(S2pl-S2mi);            % sign for column k
end

%clean the mess
clear tmp k ind S2pl S2mi