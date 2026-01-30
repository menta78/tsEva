function [jpdf, jcdf] = tsCopulaYearExtrDistribution(retPeriod, copulaParam, varargin)

  args.computeCdf = false;
  args = tsEasyParseNamedArgs(varargin, args);
  computeCdf = args.computeCdf;

  nsrs = copulaParam.nSeries;
  
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

  
  copulaFamily = copulaParam.family;
  
  if strcmpi(copulaFamily, 'gaussian')
    % normal copula
    rho = copulaParam.rho;
    jpdf = copulapdf(copulaFamily, probOut, rho);
    if computeCdf
      jcdf = copulacdf(copulaFamily, probOut, rho);
    end
  elseif strcmpi(copulaFamily, 't')
    % t copula
    rho = copulaParam.rho;
    nu = copulaParam.nu;
    jpdf = copulapdf(copulaFamily, probOut, rho, nu);
    if computeCdf
      jcdf = copulacdf(copulaFamily, probOut, rho, nu);
    end
  elseif strcmpi(copulaFamily, 'gumbel') || strcmpi(copulaFamily, 'clayton') || strcmpi(copulaFamily, 'frank')
    % one of the archimedean copulas
    theta = copulaParam.theta;
    jpdf = copulapdf(copulaFamily, probOut, theta);
    if computeCdf
      jcdf = copulacdf(copulaFamily, probOut, theta);
    end
  else
    error(['copulaFamily not supported: ' copulaFamily]);
  end
  jpdf = reshape(jpdf, nrperCArr{:});
  if computeCdf
    jcdf = reshape(jcdf, nrperCArr{:});
  else
    jcdf = [];
  end
