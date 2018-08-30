function [retLev, jdist, cplParam] = tsCopulaYearExtrFit(retPeriod, retLev, retLevError, yMax, varargin)

  args.copulaFamily = 'gaussian'; % can be gaussian, t, gumbel, clayton, frank
  args = tsEasyParseNamedArgs(varargin, args);
  copulaFamily = args.copulaFamily;

  nsrs = size(yMax, 2);
  yRetPer = zeros(size(yMax))*nan;
  for isrs = 1:nsrs
    [yRetPer(:,isrs), yRetPerMax, yRetPerMin] = tsGetReturnPeriodOfLevel(retPeriod, retLev(:,isrs), retLevError(:,isrs), yMax(:,isrs), varargin{:});
  end
  yRetPer(yRetPer == inf) = .1;
  yRetPer(yRetPer == -inf) = .1;
  yRetPer(isnan(yRetPer)) = .1;

  yProb = 1 - 1./yRetPer;
  yProb(yProb <= 0) = min(yProb(yProb > 0));
  yProb(yProb >= 1) = .9999;
  yProb(isnan(yProb)) = .0001;

  nrper = length(retPeriod);
  retPerCArr = cell(nsrs, 1);
  nrperCArr = cell(nsrs, 1);
  for isrs = 1:nsrs
    retPerCArr{isrs} = retPeriod;
    nrperCArr{isrs} = nrper;
  end
  retPerMtxCArr = cell(nsrs, 1);
  [retPerMtxCArr{:}] = ndgrid(retPerCArr{:});
  retPerOutCArr = cellfun(@(rp) rp(:), retPerMtxCArr, 'UniformOutput', false);
  retPerOut = horzcat(retPerOutCArr{:});
  probOut = 1 - 1./retPerOut;

  cplParam.family = copulaFamily;
  if strcmpi(copulaFamily, 'gaussian')
    % normal copula
    rho = copulafit(copulaFamily, yProb);
    jdist = copulacdf(copulaFamily, probOut, rho);
    cplParam.rho = rho;
  elseif strcmpi(copulaFamily, 't')
    % t copula
    [rho, nu] = copulafit(copulaFamily, yProb);
    jdist = copulacdf(copulaFamily, probOut, rho, nu);
    cplParam.rho = rho;
    cplParam.nu = nu;
  elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    % one of the archimedean copulas
    [cprm, cci] = copulafit(copulaFamily, yProb);
    jdist = copulacdf(copulaFamily, probOut, cprm);
    cplParam.theta = cprm;
    cplParam.thetaCI = cci;
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end
  jdist = reshape(jdist, nrperCArr{:});
  


