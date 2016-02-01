function phandles = plotReturnLevelsGEVStationary( nonStationaryEvaParams, varargin )
% nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach
% timeIndex is the index at which the time varying analysis should be
% estimated.

statParams = nonStationaryEvaParams(1).stationaryParams;
epsilon = statParams.parameters(1);
sigma = statParams.parameters(2);
mu = statParams.parameters(3);
epsilonStdErr = abs(statParams.paramCIs(1, 1) - epsilon);
sigmaStdErr = abs(statParams.paramCIs(1, 2) - sigma);
muStdErr = abs(statParams.paramCIs(1, 3) - mu);

phandles = plotReturnLevelsGEV(epsilon, sigma, mu, epsilonStdErr, sigmaStdErr, muStdErr, varargin{:});
end

