function [EVmeta,EVdata,isValid] = tsEVstatistics(pointData, varargin)
% Evangelos Voukouvalas, Michalis Vousdoukas 2015
% gevMaxima can be annual or monthly. annual by default

% in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
% interval

  defvals={[],[],.68};
  minvars=1;
  EVdata = [];
  minGEVSample = 10;

  args.alphaCI = .95;
  args.gevMaxima = 'annual';
  args.gevType = 'GEV' % can be 'GEV' or 'Gumbel'
  args.evdType = {'GEV', 'GPD'};
  args = tsEasyParseNamedArgs(varargin, args);
  gevMaxima = args.gevMaxima;
  gevType = args.gevType;
  alphaCI = args.alphaCI;
  evdType = args.evdType;

  isValid = true;

  %% Basic data

  Tr=[5,10,20,50,100,200,500,1000];

  EVmeta.Tr=Tr;

  [nyears]=size(pointData.annualMax);

  %% stationary GEV

  imethod=1;
  methodname='GEVstat';

  paramEsts=nan(1,3);
  paramCIs=nan*ones(2,3);
  rlvls=nan(1,length(Tr));

  if ~isempty(find(strcmpi(evdType, 'GEV'), 1)) && ~isempty(pointData.annualMax)
    if strcmpi(gevMaxima, 'annual')
      tmpmat=pointData.annualMax(:);
    elseif strcmpi(gevMaxima, 'monthly')
      tmpmat=pointData.monthlyMax(:);
    else
      error(['tsEVstatistics: invalid gevMaxima type: ' gevMaxima]);
    end

    iIN=~isnan(tmpmat);

    if sum(iIN)>=minGEVSample

      tmp=tmpmat(iIN);

      if strcmpi(gevType, 'GEV')
        [paramEsts,paramCIs]=gevfit(tmp, alphaCI);
      elseif strcmpi(gevType, 'Gumbel')
        paramEsts(1) = 0;
        paramCIs(:,1) = 0;
        [paramEsts(:,[3,2]),paramCIs(:,[3,2])]=evfit(-tmp, alphaCI);
        paramEsts(3) = -paramEsts(3);
        paramCIs([1,2], 3) = -paramCIs([2,1], 3);
      else
        error(['tsEVstatistics: invalid gevType: ' gevType '. Can be only GEV or Gumbel']);
      end
      % paramEsts(jj,1): shape param
      % paramEsts(jj,2): scale param
      % paramEsts(jj,3): location param

      % the second parameter returned by gevfit is the 95% confidence
      % interval

      rlvls(1,:) = gevinv(1-1./Tr,paramEsts(1,1),paramEsts(1,2),paramEsts(1,3));

    else 

      disp('Skipping...')
      isValid = false;

    end

    EVdata(imethod).method=methodname;
    EVdata(imethod).values=rlvls;
    EVdata(imethod).parameters=paramEsts;
    EVdata(imethod).paramCIs = paramCIs;

  else
    EVdata(imethod).method = methodname;
    EVdata(imethod).values = [];
    EVdata(imethod).parameters = [];
    EVdata(imethod).paramCIs = [];
  end

  %% stationary GPD

  imethod=2;
  methodname='GPDstat';

  paramEstsall=nan(1,6);
  paramCIs=nan*ones(2,2);
  rlvls=nan(1,length(Tr));

  try
    if ~isempty(find(strcmpi(evdType, 'GPD'), 1))
      ik = 1;
      d1=pointData.POT(ik).peaks-pointData.POT(ik).threshold;

      [paramEsts,paramCIs]=gpfit(d1, alphaCI);
      % shape parameter
      ksi=paramEsts(1);
      % scale parameter
      sgm=paramEsts(2);

      if ksi < -.5
        % computing anyway the confidence interval (in a rough way)
        probs = [alphaCI/2; 1-alphaCI/2];
        [~, acov] = gplike([ksi sgm], d1);
        se = sqrt(diag(acov))';

        % Compute the CI for k using a normal distribution for khat.
        kci = norminv(probs, ksi, se(1));
        % VERY ROUGHT: minimizing the lower boundary of kci to -1
        kci(kci < -1) = -1;

        % Compute the CI for sigma using a normal approximation for
        % log(sigmahat), and transform back to the original scale.
        % se(log(sigmahat)) is se(sigmahat) / sigmahat.
        lnsigci = norminv(probs, log(sgm), se(2)./sgm);

        paramCIs = [kci exp(lnsigci)];

      end

      % paramCIs: 95% confidence interval


      paramEstsall(ik,:)=[sgm ksi pointData.POT(ik).threshold length(d1) length(pointData.POT(ik).peaks) pointData.POT(ik).percentile];

      rlvls(ik,:) = pointData.POT(ik).threshold+(sgm/ksi).*((((length(d1)/length(pointData.POT(ik).peaks))*(1./Tr)).^(-ksi))-1);

    else 

      EVdata(imethod).method = methodname;
      EVdata(imethod).values = [];
      EVdata(imethod).parameters = [];
      EVdata(imethod).paramCIs = [];

    end
  catch err
      disp(getReport(err));
      paramEstsall(ik,:)=[0 0 0 0 0 0];
      rlvls(ik,:) = zeros(1, length(Tr));
      isValid = false;
  end

  EVdata(imethod).method=methodname;
  EVdata(imethod).values=rlvls;
  EVdata(imethod).parameters=paramEstsall;
  EVdata(imethod).paramCIs = flipdim(paramCIs, 2);

end

