function [EVmeta,EVdata,isValid]=tsEVstatistics(pointData, varargin)
% Evangelos Voukouvalas, Michalis Vousdoukas 2015
% gevMaxima can be annual or monthly. annual by default

% in a gaussian approximation alphaCI~68% corresponds to 1 sigma confidence
% interval
defvals={[],[],.68};
minvars=1;
EVdata = [];

args.alphaCI = .95;
args.gevMaxima = 'annual';
args = tsEasyParseNamedArgs(varargin, args);
gevMaxima = args.gevMaxima;
alphaCI = args.alphaCI;

isValid = true;

%% Basic data

Tr=[5,10,20,50,100,200,500,1000];

EVmeta.Tr=Tr;
try
  EVmeta.lon=pointData.lon;
  EVmeta.lat=pointData.lat;
catch
  disp('lon lat not found. Cannot set them to metadata.');
end

[npoints,nyears]=size(pointData.annualMax);

%% stationary GEV

imethod=1;
methodname='GEVstat';

paramEsts=nan(npoints,3);
paramCIs=nan*ones(2,3);
rlvls=nan(npoints,length(Tr));

criterio=zeros(npoints,1);

if length(EVdata)<imethod
    
    for jj=1:npoints;

        disp(['GEV - ' num2str(jj/npoints*100) '%'])

        if strcmpi(gevMaxima, 'annual')
            tmpmat=pointData.annualMax(jj,:);
        elseif strcmpi(gevMaxima, 'monthly')
            tmpmat=pointData.monthlyMax(jj,:);
        else
            error(['tsEVstatistics: invalid gevMaxima type: ' gevMaxima]);
        end

        iIN=~isnan(tmpmat);

        if sum(iIN)>=20

            criterio(jj)=1;

            tmp=tmpmat(iIN);

            [paramEsts(jj,1:3),paramCIs]=gevfit(tmp, alphaCI);
            % paramEsts(jj,1): shape param
            % paramEsts(jj,2): scale param
            % paramEsts(jj,3): location param
            
            % the second parameter returned by gevfit is the 95% confidence
            % interval

            rlvls(jj,:) = gevinv(1-1./Tr,paramEsts(jj,1),paramEsts(jj,2),paramEsts(jj,3));

        else 

            disp('Skipping...')
            isValid = false;

        end


    end

    EVdata(imethod).method=methodname;
    EVdata(imethod).values=rlvls;
    EVdata(imethod).parameters=paramEsts;
    EVdata(imethod).paramCIs = paramCIs;
    
end

%% stationary GPD

imethod=2;
methodname='GPDstat';

paramEstsall=nan(npoints,6);
paramCIs=nan*ones(2,2);
rlvls=nan(npoints,length(Tr));

if length(EVdata)<imethod

    for ik=1:npoints

        disp(['GPD - ' num2str(ik/npoints*100) '%'])

        try
            if criterio(ik)==1

                d1=pointData.POT(ik).peaks-pointData.POT(ik).threshold;

                [paramEsts,paramCIs]=gpfit(d1, alphaCI);
                % shape parameter
                ksi=paramEsts(1);
                % scale parameter
                sgm=paramEsts(2);
                % paramCIs: 95% confidence interval


                paramEstsall(ik,:)=[sgm ksi pointData.POT(ik).threshold length(d1) length(pointData.POT(ik).peaks) pointData.POT(ik).percentile];

                rlvls(ik,:) = pointData.POT(ik).threshold+(sgm/ksi).*((((length(d1)/length(pointData.POT(ik).peaks))*(1./Tr)).^(-ksi))-1);

                else 

                disp('Skipping...')

            end
        catch err
            disp(getReport(err));
            paramEstsall(ik,:)=[0 0 0 0 0 0];
            rlvls(ik,:) = zeros(1, length(Tr));
            isValid = false;
        end
    end

    EVdata(imethod).method=methodname;
    EVdata(imethod).values=rlvls;
    EVdata(imethod).parameters=paramEstsall;
    EVdata(imethod).paramCIs = flipdim(paramCIs, 2);
end

% 
% %% stationary GEV (Wafo)
% 
% imethod=3;
% methodname='GEVstat_wafo';
% 
% paramEsts=nan(npoints,3);
% rlvls=nan(npoints,length(Tr));
% 
% criterio=zeros(npoints,1);
%     
% if length(EVdata)<imethod
% 
%     for iw=1:npoints;
% 
%         disp(['GEV_wafo - ' num2str(iw/npoints*100) '%'])
% 
%         tmpmat=pointData.annualMax(iw,:);
% 
%         iIN=~isnan(tmpmat);
% 
%         if sum(iIN)>=20
% 
%             criterio(iw)=1;
% 
%             tmp=tmpmat(iIN);
% 
% 
%             tmp_wafo=fitgev(tmp,'method','ml','plotflag',0);
% 
%             paramEsts(iw,1:3)=tmp_wafo.params;
%             rlvls(iw,:) = invgev(1./Tr,tmp_wafo,'lowertail',false,'proflog',true);
%         else 
% 
%             disp('Skipping...')
% 
%         end
% 
% 
%     end
% 
%     EVdata(imethod).method=methodname;
%     EVdata(imethod).values=rlvls;
%     EVdata(imethod).parameters=paramEsts;
% 
% end
% 
% %%  stationary GPD (Wafo)
% 
% imethod=4;
% methodname='GPDstat_wafo';
% 
% paramEstsall=nan(npoints,3);
% rlvls=nan(npoints,length(Tr));
% 
% if length(EVdata)<imethod
% 
%     for jw=1:npoints
% 
%         disp(['GPD_wafo - ' num2str(jw/npoints*100) '%'])
% 
%         if criterio(jw)==1
% 
%     %         d1=pointData.POT(jw).peaks-pointData.POT(jw).threshold;
%             lambda = numel(pointData.POT(jw))/length(pointData.years);
%             phat = fitgenpar(pointData.POT(jw).peaks, 'fixpar',[nan,nan,pointData.POT(jw).threshold], 'method','ml');
%             [xr,~,~] = invgenpar(1./(lambda*Tr),phat,'lowertail',false,'alpha', 0.05);    
%     %         
%     %         [paramEsts,paramCIs]=gpfit(d1);
%     %         ksi=paramEsts(1);
%     %         sgm=paramEsts(2);
% 
%             paramEstsall(jw,:)=[phat.params numel(pointData.POT(jw)) length(pointData.years)]
%             
%             rlvls(jw,:) = xr;
% 
%             else 
% 
%             disp('Skipping...')
% 
%             aphat(jw)=phat;
%             
%         end
%     end
% 
%     EVdata(imethod).method=methodname;
%     EVdata(imethod).values=rlvls;
%     EVdata(imethod).parameters=paramEstsall;
%     
% end
% 
% EVmeta.GPDwafo_phat=aphat;