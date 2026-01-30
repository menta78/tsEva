function [monteCarloRsmpl, resampleProb, resampleRetPer] = tsCopulaYearExtrRnd(retPeriod, retLev, copulaParam, nResample, varargin)

  args.logExtrapRetlev = true;
  args = tsEasyParseNamedArgs(varargin, args);
  logExtrapRetlev = args.logExtrapRetlev;

  copulaFamily = copulaParam.family;
  if strcmpi(copulaFamily, 'gaussian')
    resampleProb = copularnd(copulaFamily, copulaParam.rho, nResample);
  elseif strcmpi(copulaFamily, 't')
    resampleProb = copularnd(copulaFamily, copulaParam.rho, copulaParam.nu, nResample);
  elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    resampleProb = copularnd(copulaFamily, copulaParam.theta, nResample);
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end


  resampleRetPer = 1./(1-resampleProb);
  
  retLev((retLev == inf) | (retLev == -inf)) = nan;
  
  monteCarloRsmpl = ones(size(resampleRetPer))*nan;
  ncpl = size(resampleRetPer, 2);
  for ivar = 1:ncpl
    monteCarloRsmpl(:,ivar) = tsInterp1Extrap(retPeriod, retLev(:,ivar), resampleRetPer(:,ivar), logExtrapRetlev);
    monteCarloRsmpl(isnan(monteCarloRsmpl(:,ivar)),ivar) = nanmean(monteCarloRsmpl(:,ivar));
  end


