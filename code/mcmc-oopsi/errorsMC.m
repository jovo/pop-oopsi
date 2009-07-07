function rep=errorsMC(P,est,tgt,v)
%Calculates mean error in positions, and mean fractions of false positive 
%and negative spikes from collection of "targets" spike trains and their 
%corresponding "estimators".
%REP=errorsMC(P,ESTIMATORS,TARGETS[,V]) computes THE above measures by 
%aligning estimated and target spike trains using bipartite matching 
%solved with Hungarian algorithm. Both "targets" and "estimators" should 
%be cell-arrays of spike times. For example, targets may be cell-array of 
%100 arrays specifying times of spikes in each actual spike train, and 
%estimators may be cell-array of 100 arrays specifying times of spikes in 
%estimated trains. P is a structure of parameters containing following 
%fields:
% P.mismatch_penalty - specifying distance between spikes that are 
%                      considered mismatched, and 
% P.expected_slack   - specifying an upper expected bound on the fraction 
%                      of mismatched spikes. 
%v specifies verbosity. Y. Mishchenko 2009 Columbia Un

DM=P.mismatch_penalty^2;          %spikes discrepancy penalty 
slack_f=P.expected_slack;         %expected fraction of mismatches
if(nargin<4 || isempty(v)) v=1; end%verbosity

errs=zeros(3,length(est));
for i=1:length(est)
  %prepare spike trains and match costs
  X=tgt{i}(:);                    %true spike train
  Y=est{i}(:);                    %estimator spike train
  A=-2*repmat(X,1,length(Y)).*repmat(Y',length(X),1);
  A=A+repmat(X,1,length(Y)).^2;
  A=A+repmat(Y',length(X),1).^2;  %matrix of distances for matching
  
  slack_n=ceil(slack_f*length(X));%add slack=mismatches  
  A=[A,repmat(DM,length(X),slack_n)];
  slack_n=ceil(slack_f*length(Y));  
  A=[A;repmat(DM,slack_n,size(A,2))];
  
  %compute optimal assignments with Hungarion alg
  [a,b]=assignmentoptimal(A);     %Hungarion alg
  a=a(1:length(X)); a(a>length(Y))=0;
  x=(X(a>0)-Y(a(a>0))).^2;        %position errors
  dn_lost=length(X)-sum(a>0);     %mismatches
  dn_false=length(Y)-sum(a>0);
  errs(1,i)=sqrt(sum(x(x<DM))/sum(x<DM));%mean position error
  errs(2,i)=(dn_lost)/length(X);  %spikes lost
  errs(3,i)=(dn_false)/length(X); %spikes false
  if(v && mod(i,10)==0) fprintf('\b\b\b%3i',i); end
end
if(v) fprintf('\n'); end

rep=[];
rep.position_error  = mean(errs(1,:));
rep.false_negatives = mean(errs(2,:));
rep.false_positives = mean(errs(3,:));
rep.data=errs;
