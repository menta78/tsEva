function [retLev, jpdf, cplParam, jcdf] = tsCopulaYearExtrFit(retPeriod, retLev, yMax, varargin)

  args.copulaFamily = 'gaussian'; % can be gaussian, t, gumbel, clayton, frank
  args.computeCdf = false;
  args = tsEasyParseNamedArgs(varargin, args);
  copulaFamily = args.copulaFamily;
  computeCdf = args.computeCdf;

  nsrs = size(yMax, 2);
  nyr = size(yMax, 1);
  nretPer = length(retPeriod);
  szsRetLev = size(retLev);
  yRetPer = zeros(size(yMax))*nan;
  fakeErr = zeros(size(retLev));
  
  ndimRetLev = length(szsRetLev);
  
  if ndimRetLev == 2
    % this is a stationary set of return levels
    for isrs = 1:nsrs
      [yRetPer(:,isrs), ~, ~] = tsGetReturnPeriodOfLevel(retPeriod, retLev(:,isrs), fakeErr(:,isrs), yMax(:,isrs), varargin{:});
    end
    yRetPer(yRetPer == inf) = .1;
    yRetPer(yRetPer == -inf) = .1;
    yRetPer(isnan(yRetPer)) = .1;
  else
    % this is a non-stationary set of return levels
    if (size(szsRetLev, 1) ~= nyr) || (size(szsRetLev, 2) ~= nretPer) || (size(szsRetLev, 3) ~= nsrs)
      error('tsCopulaYearExtrFit: for non-stationary, retLev should be dimensioned as (nTime x nReturnPeriod x nSeries)');
    end
  end

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
  cplParam.familyId = tsCopulaGetFamilyId(copulaFamily);
  if strcmpi(copulaFamily, 'gaussian')
    % normal copula
    rho = copulafit(copulaFamily, yProb);
    jpdf = copulapdf(copulaFamily, probOut, rho);
    if computeCdf
      jcdf = copulacdf(copulaFamily, probOut, rho);
    end
    cplParam.rho = rho;
  elseif strcmpi(copulaFamily, 't')
    % t copula
    [rho, nu] = copulafit(copulaFamily, yProb);
    jpdf = copulapdf(copulaFamily, probOut, rho, nu);
    if computeCdf
      jcdf = copulacdf(copulaFamily, probOut, rho, nu);
    end
    cplParam.rho = rho;
    cplParam.nu = nu;
  elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    % one of the archimedean copulas
    [cprm, cci] = copulafit(copulaFamily, yProb);
    jpdf = copulapdf(copulaFamily, probOut, cprm);
    if computeCdf
      jcdf = copulacdf(copulaFamily, probOut, cprm);
    end
    cplParam.theta = cprm;
    cplParam.cci = cci;
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end
  jpdf = reshape(jpdf, nrperCArr{:});
  if computeCdf
    jcdf = reshape(jcdf, nrperCArr{:});
  else
    jcdf = [];
  end
  


