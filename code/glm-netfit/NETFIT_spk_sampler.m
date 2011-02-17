%SCRIPT TO DRAW SAMPLE OF SPIKE HISTORIES FOR ESTIMATING W

switch(mode)
  case 'gibbs'            %PF-based indep. approx. within Gibbs
    NETFIT_sampler_indep
  case 'neal'             %neals corrected MCMC method within Gibbs
    NETFIT_sampler_neal
  case 'iid'              %PF-based independent approximation
    NETFIT_sampler_indep 
  case 'base'             %BASELINE
    NETFIT_sampler_base
  otherwise
    error('Unknown main sampler type');
end
