function phandles = plotReturnLevelsGPDStationary( nonStationaryEvaParams, varargin )
% nonStationaryEvaParams, stationaryTransformData are as these returned by function nonStationaryEvaJRCApproach
% timeIndex is the index at which the time varying analysis should be
% estimated.

statParams = nonStationaryEvaParams(2).stationaryParams;
epsilon = statParams.parameters(2);
sigma = statParams.parameters(1);
threshold = statParams.parameters(3);
dtSampleYears = nonStationaryEvaParams(2).parameters.timeDeltaYears;
epsilonStdErr = abs(statParams.paramCIs(1, 2) - epsilon);
sigmaStdErr = abs(statParams.paramCIs(1, 1) - sigma);
thresholdStdErr = 0;

phandles = plotReturnLevelsGPD(epsilon, sigma, threshold, dtSampleYears, epsilonStdErr, sigmaStdErr, thresholdStdErr, varargin{:});
end

