%SCRIPT EVALUATES OUTER EM-LOOP TERMINATION CONDITIONS
eps=0.1;        % target percentile accuracy
cntMax=5;       % max loops counter

switch mode
  case {'base','iid'}          % if baseline or iid, only once !
    iflg=0;
  case {'gibbs','neal'}        % if neal or indep.-gibbs, stopping criterion
    dW=sum(abs(Weights_old(:)-Weights_glm(:)))/sum(abs(Weights_glm(:)));
    Weights_old=Weights_glm;
    cnt=cnt+1;  
  
    if(dW<eps) iflg=0; end
    if(cnt>cntMax) iflg=0; end
  otherwise
    error('Unknown main sampler type');
end

