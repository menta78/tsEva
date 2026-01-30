function [returnLevels, returnLevelsErr] = tsEvaComputeReturnLevelsGPD( epsilon, sigma, threshold, epsilonStdErr, sigmaStdErr, thresholdStdErr, nPeaks, sampleTimeHorizon, returnPeriods )
% reference: Stuart Coles 2001, pag 81.
% sampleTimeHorizon and returnPeriods must be in the same units, e.g. years

% this is how XX is computed in Coles, which is more error-prone
% thresholdExceedProbability = 1 - percentile/100;
% m = returnPeriods/double(dtSample);
% XX = m*thresholdExceedProbability;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  X0 = nPeaks/sampleTimeHorizon;
  XX = X0*returnPeriods;

  npars = length(sigma);
  nt = length(returnPeriods);
  XX_ = XX(ones(1,npars),:);
  sigma_ = sigma(:,ones(1,nt));
  sigmaStdErr_ = sigmaStdErr(:,ones(1,nt));
  threshold_ = threshold(:,ones(1,nt));
  thresholdStdErr_ = thresholdStdErr(:,ones(1,nt));
 %if epsilon ~= 0
  if abs(epsilon) >= 1e-7
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
      
      returnLevelsErr = (  (dxm_u.*thresholdStdErr_).^2  +  (dxm_sigma.*sigmaStdErr_).^2  ).^.5;
      %%
  end
end

