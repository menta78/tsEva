function [nonStationaryEvaParamsEns, stationaryTransformDataEns]  = tsEnsemble( nonStatEvaParamsArray, stationaryTransformDataArray )
  % calls tsEnsembleEvaParams and tsEnsembleStatTransfData.
  % ASSUMING THAT ALL THE SERIES HAVE THE SAME LENGTH 
  % AND UNIFORM TIME STAMPS
  
  nonStationaryEvaParamsEns = tsEnsembleEvaParams(nonStatEvaParamsArray);
  stationaryTransformDataEns = tsEnsembleStatTransfData(stationaryTransformDataArray);

end

