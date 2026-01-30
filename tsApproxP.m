function [Pval]=tsApproxP(N,copulaFamily,rho,nu,snSample,s2Sample)

f = waitbar(0,'1','Name','Approximating p value...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);


 N=1000;
    Snk=[];
    for k=1:N
       
        if strcmpi(copulaFamily,'gaussian')
            psur=copularnd(copulaFamily,rho,N);
            Ux = tsPseudoObservations(psur);
            rhox = copulafit(copulaFamily, Ux);
            Yx = copulacdf(copulaFamily, Ux, rhox);
        elseif strcmpi(copulaFamily,'t')
            psur=copularnd(copulaFamily,rho,nu,N);
            Ux = tsPseudoObservations(psur);
            [rhox,nux] = copulafit(copulaFamily, Ux);
            Yx = copulacdf(copulaFamily, Ux, rhox, nux);
        elseif strcmpi(copulaFamily,'clayton')|| strcmpi(copulaFamily,'frank') || strcmpi(copulaFamily,'gumbel')
            if s2Sample<3
                psur=copularnd(copulaFamily,rho,N);
            
            else
                error('copula not detected')

            end
            Ux = tsPseudoObservations(psur);
            alphax = copulafit(copulaFamily, Ux);
            Yx = copulacdf(copulaFamily, Ux, alphax);

        end

        snk= sum((tsEmpirical(Ux) - Yx) .^ 2);
        Snk=[Snk;snk];

        waitbar(k/N,f,sprintf('%2.3f',snSample))
    end

    b=length(find(Snk>=snSample));
    Pval=(0.5+b)/(N+1);

delete(f)
end