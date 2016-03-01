function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, percentile, epsilonStdErr, sigmaStdErr, thresholdStdErr, dtSample, returnPeriods )
% reference: Stuart Coles 2001, pag 81.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!
% dtSample IS NOT the time interval of the sampling, but the minimum time
% interval between 2 peaks !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  thresholdExceedProbability = 1 - percentile/100;
  %thresholdExceedProbability = 1;
  m = returnPeriods/dtSample;
  XX = m*thresholdExceedProbability;
  % uniforming dimensions of XX, sigma, mu
  npars = length(sigma);
  nt = length(returnPeriods);
  XX_ = XX(ones(1,npars),:);
  sigma_ = sigma(:,ones(1,nt));
  sigmaStdErr_ = sigmaStdErr(:,ones(1,nt));
  threshold_ = threshold(:,ones(1,nt));
  thresholdStdErr_ = thresholdStdErr(:,ones(1,nt));
  if epsilon ~= 0
      %% estimating the return levels
      returnLevels = threshold_ + sigma_./epsilon.*( (XX_).^epsilon - 1 );
      %% estimating the error
      % estimating the differential of returnLevels to the parameters
      % !! ASSUMING NON ZERO ERROR ON THE THRESHOLD AND 0 ERROR ON THE PERCENTILE.
      % THIS IS NOT COMPLETELY CORRECT BECAUSE THE PERCENTILE DEPENDS ON THE
      % THRESHOLD AND HAS THEREFORE AN ERROR RELATED TO THAT OF THE
      % THRESHOLD
      dxm_u = 1;
      dxm_sigma = 1/epsilon * (XX_.^epsilon - 1);
      dxm_epsilon = - sigma_/epsilon^2 .* ( (XX_).^epsilon - 1 ) + sigma_/epsilon .* log(XX_).*XX_.^epsilon;
      
      returnLevelsErr = (  (dxm_u.*thresholdStdErr_).^2  +  (dxm_sigma.*sigmaStdErr_).^2  +  (dxm_epsilon.*epsilonStdErr).^2  ).^.5;
      %%
  else
      returnLevels = threshold_ + sigma_.*log(XX_);
      %% estimating the error
      % estimating the differential of returnLevels to the parameters
      % !! ASSUMING NON ZERO ERROR ON THE THRESHOLD, 0 ERROR ON THE
      % PERCENTILE AND 0 ERROR ON EPSILON.
      % THIS IS NOT COMPLETELY CORRECT BECAUSE THE PERCENTILE DEPENDS ON THE
      % THRESHOLD AND HAS THEREFORE AN ERROR RELATED TO THAT OF THE
      % THRESHOLD
      dxm_u = 1;
      dxm_sigma = log(XX_);
      
      returnLevelsErr = (  (dxm_u*thresholdStdErr_).^2  +  (dxm_sigma*sigmaStdErr_).^2  ).^.5;
      %%
  end
end

