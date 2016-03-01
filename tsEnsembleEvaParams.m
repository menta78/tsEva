function nonStatEvaParamsEnsemble = tsEnsembleEvaParams( nonStatEvaParamsArray )
  % From a cell array of non nonStatEvaParams computes the average
  % nonStatEvaParamsArray is a cell array of nonStatEvaParams, the object
  % type returned by tsEvaNonStationary
  % ASSUMING THAT ALL THE SERIES HAVE THE SAME LENGTH
  % AND UNIFORM TIME STAMPS
  
  n = length(nonStatEvaParamsArray);
  for ii = 1:n
    nsep = nonStatEvaParamsArray{ii};
    if ii == 1
      ensep(1).method = 'GEVstat';
      ensep(1).parameters = nsep(1).parameters;
      ensep(1).paramErr = nsep(1).paramErr;
      ensep(2).method = 'GPDstat';
      ensep(2).parameters = nsep(2).parameters;
      ensep(2).paramErr = nsep(2).paramErr;
      ensep(2).timeDelta = nsep(2).timeDelta;
      ensep(2).timeDeltaYears = nsep(2).timeDeltaYears;
    else
      ensep(1).parameters.epsilon = ensep(1).parameters.epsilon + nsep(1).parameters.epsilon;
      ensep(1).parameters.sigma = ensep(1).parameters.sigma + nsep(1).parameters.sigma;
      ensep(1).parameters.mu = ensep(1).parameters.mu + nsep(1).parameters.mu;
      ensep(1).paramErr.epsilonErr = ensep(1).paramErr.epsilonErr + nsep(1).paramErr.epsilonErr;
      ensep(1).paramErr.sigmaErr = ensep(1).paramErr.sigmaErr + nsep(1).paramErr.sigmaErr;
      ensep(1).paramErr.muErr = ensep(1).paramErr.muErr + nsep(1).paramErr.muErr;

      ensep(2).parameters.epsilon = ensep(2).parameters.epsilon + nsep(2).parameters.epsilon;
      ensep(2).parameters.sigma = ensep(2).parameters.sigma + nsep(2).parameters.sigma;
      ensep(2).parameters.threshold = ensep(2).parameters.threshold + nsep(2).parameters.threshold;
      ensep(2).parameters.percentile = ensep(2).parameters.percentile + nsep(2).parameters.percentile;
      ensep(2).paramErr.epsilonErr = ensep(2).paramErr.epsilonErr + nsep(2).paramErr.epsilonErr;
      ensep(2).paramErr.sigmaErr = ensep(2).paramErr.sigmaErr + nsep(2).paramErr.sigmaErr;
      ensep(2).paramErr.thresholdErr = ensep(2).paramErr.thresholdErr + nsep(2).paramErr.thresholdErr;
    end
  end

  ensep(1).parameters.epsilon = ensep(1).parameters.epsilon/n;
  ensep(1).parameters.sigma = ensep(1).parameters.sigma/n;
  ensep(1).parameters.mu = ensep(1).parameters.mu/n;
  ensep(1).paramErr.epsilonErr = ensep(1).paramErr.epsilonErr/n;
  ensep(1).paramErr.sigmaErr = ensep(1).paramErr.sigmaErr/n;
  ensep(1).paramErr.muErr = ensep(1).paramErr.muErr/n;
  ensep(2).parameters.epsilon = ensep(2).parameters.epsilon/n;
  ensep(2).parameters.sigma = ensep(2).parameters.sigma/n;
  ensep(2).parameters.threshold = ensep(2).parameters.threshold/n;
  ensep(2).parameters.percentile = ensep(2).parameters.percentile/n;
  ensep(2).paramErr.epsilonErr = ensep(2).paramErr.epsilonErr/n;
  ensep(2).paramErr.sigmaErr = ensep(2).paramErr.sigmaErr/n;
  ensep(2).paramErr.thresholdErr = ensep(2).paramErr.thresholdErr/n;
  
  nonStatEvaParamsEnsemble = ensep;
end

