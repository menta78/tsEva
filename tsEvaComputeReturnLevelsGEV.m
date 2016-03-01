function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGEV( epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, returnPeriodsInYears, varargin )
  args.maxType = 'monthly';
  args = tsEasyParseNamedArgs(varargin, args);
  
  if strcmpi(args.maxType, 'monthly')
      returnPeriods = returnPeriodsInYears*12;
  elseif strcmpi(args.maxType, 'yearly')
      returnPeriods = returnPeriodsInYears;
  else
      error('computeReturnLevelsGEV: unsupported maxtype. The only supported values are monthly and yearly');
  end
      
  
% reference: Stuart Coles 2001, pag 49.
  yp = -log(1 - 1./returnPeriods);
  
% uniforming dimensions of yp, sigma, mu
  npars = length(sigma);
  nt = length(returnPeriodsInYears);
  yp = yp(ones(1,npars),:);
  sigma_ = sigma(:,ones(1,nt));
  sigmaStdErr_ = sigmaStdErr(:,ones(1,nt));
  mu_ = mu(:,ones(1,nt));
  muStdErr_ = muStdErr(:,ones(1,nt));
  if epsilon ~= 0
      %% estimating the return levels
      returnLevels = mu_ - sigma_./epsilon.*( 1 - yp.^(-epsilon) );
      %% estimating the error
      % estimating the differential of returnLevels to the parameters
      dxm_mu = 1;
      dxm_sigma = 1/epsilon * (1 - yp.^(-epsilon));
      dxm_epsilon = sigma_./epsilon.^2 .* ( 1 - yp.^(-epsilon) ) - sigma_./epsilon .* log(yp).*yp.^(-epsilon);
      
      returnLevelsErr = (  (dxm_mu.*muStdErr_).^2  +  (dxm_sigma.*sigmaStdErr_).^2  +  (dxm_epsilon.*epsilonStdErr).^2  ).^.5;
      %%
  else
      returnLevels = mu_ - sigma_.*log(yp);
      %% estimating the error
      % estimating the differential of returnLevels to the parameters
      dxm_u = 1;
      dxm_sigma = log(yp);
      
      returnLevelsErr = (  (dxm_u.*muStdErr_).^2  +  (dxm_sigma.*sigmaStdErr_).^2  ).^.5;
      %%
  end
end

