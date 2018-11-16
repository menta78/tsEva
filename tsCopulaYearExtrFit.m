function [retLev, copulaParam, yRetPer, yProb] = tsCopulaYearExtrFit(retPeriod, retLev, yMax, varargin)

  args.copulaFamily = 'gaussian'; % can be gaussian, t, gumbel, clayton, frank
  args = tsEasyParseNamedArgs(varargin, args);
  copulaFamily = args.copulaFamily;

  nsrs = size(yMax, 2);
  nyr = size(yMax, 1);
  nretPer = length(retPeriod);
  szsRetLev = size(retLev);
  yRetPer = zeros(size(yMax))*nan;
  fakeErr = zeros(size(retLev));
  
  ndimRetLev = length(szsRetLev);
  
  if ndimRetLev == 2
    % this is a stationary set of return levels
    if (szsRetLev(1) ~= nretPer) || (szsRetLev(2) ~= nsrs)
      error('tsCopulaYearExtrFit: for stationary, retLev should be dimensioned as (nReturnPeriod x nSeries)');
    end
    for isrs = 1:nsrs
      %[yRetPer(:,isrs), ~, ~] = tsGetReturnPeriodOfLevel(retPeriod, retLev(:,isrs), fakeErr(:,isrs), yMax(:,isrs), varargin{:});
      yRetPer(:,isrs) = interp1(retLev(:,isrs), retPeriod, yMax(:,isrs));
    end
  elseif ndimRetLev == 3
    % this is a non-stationary set of return levels
    if (szsRetLev(1) ~= nyr) || (szsRetLev(2) ~= nretPer) || (szsRetLev(3) ~= nsrs)
      error('tsCopulaYearExtrFit: for non-stationary, retLev should be dimensioned as (nTime x nReturnPeriod x nSeries)');
    end
    for isrs = 1:nsrs
      for iyr = 1:nyr
        %[yRetPer(iyr,isrs), ~, ~] = tsGetReturnPeriodOfLevel(retPeriod, retLev(iyr,:,isrs), fakeErr(iyr,:,isrs), yMax(iyr,isrs), varargin{:});
        yRetPer(iyr,isrs) = interp1(retLev(iyr,:,isrs), retPeriod, yMax(iyr,isrs));
      end
    end
  else
    error('tsCopulaYearExtrFit: retLev should be 2-dim for stationary, 3-dim for non-stationary');
  end
  
  yRetPer(yRetPer == inf) = .1;
  yRetPer(yRetPer == -inf) = .1;
  yRetPer(isnan(yRetPer)) = .1;

  yProb = 1 - 1./yRetPer;
  yProb(yProb <= 0) = min(yProb(yProb > 0));
  yProb(yProb >= 1) = .9999;
  yProb(isnan(yProb)) = .0001;

  copulaParam.family = copulaFamily;
  copulaParam.familyId = tsCopulaGetFamilyId(copulaFamily);
  copulaParam.nSeries = nsrs;
  if strcmpi(copulaFamily, 'gaussian')
    % normal copula
    rho = copulafit(copulaFamily, yProb);
    copulaParam.rho = rho;
  elseif strcmpi(copulaFamily, 't')
    % t copula
    [rho, nu] = copulafit(copulaFamily, yProb);
    copulaParam.rho = rho;
    copulaParam.nu = nu;
  elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    % one of the archimedean copulas
    [cprm, cci] = copulafit(copulaFamily, yProb);
    copulaParam.theta = cprm;
    copulaParam.cci = cci;
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end
  


