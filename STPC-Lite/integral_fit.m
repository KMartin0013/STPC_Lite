function [ESTtotal_CC,ESTalphavarall_CC,ESTtotalparams_CC,...
        ESTtotalparamerrors_CC,ESTtotalfit_CC,ESTalphavar_CC]=integral_fit(N,ESTthedates,functionintegrals,ESTsignaldelta,Cab)

if length(N)==1
    loop_range=1:N; % The total number
else
    loop_range=N; % selective number
end

for i=loop_range

    % Since Int should have units of (fn * m^2), need to go from fractional
    % sphere area to real area.  If the fn is surface density, this output is
    % in kilograms.  Then change the units from kg to Gt in METRIC tons

    %         [eigfunINT_CC] = integratebasis(CC{i},TH,N);

    %     eigfunINT = eigfunINT*4*pi*6370000^2/10^3/10^9;
    %     functionintegrals = eigfunINT;

    % Now multiply by the appropriate slepcoffs to get the months
    % This becomes alpha by months
    %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
    %functimeseries = sleptdelta(:,1:N)';

    % Here do the total sum of the data
    %m  Original code for integrating of all slepian functions
    % >total=eigfunINT*sleptdelta(:,1:N)'; (Line 526 in in slept2resid)
    %m  Revised code for each slepian frunction
    %m because >functionintegrals=eigfunINT; (Line 518 in in slept2resid)

    ESTtotal_CC(i,:)=functionintegrals(i)*ESTsignaldelta(:,i)';

    % Get the error
    %m  change diag(Cab(1:N,1:N))' to diag(Cab(i,i))'
    thevars = diag(Cab(i,i))';
    %m  change functionintegrals to functionintegrals(i)
    ESTalphavar_CC(i) = functionintegrals(i).^2.*thevars;
    % Now the combined error with covariance
    %m  change functionintegrals to functionintegrals(i)
    ESTalphavarall_CC(i) = functionintegrals(i)*Cab(i,i)*functionintegrals(i)';

    % FITTING

    % We have uniform estimated error, which will be different than the polyfit
    % estimated residuals because ours account for sinusoidal signals.  So
    % pass the new error to our function for replacement, so
    % that the fitting confidence intervals reflect that

    % Original code for integrating of all slepian functions
    %    >[fit,delta,totalparams,paramerrors] = timeseriesfit([thedates' total'],alphavarall,1,1);
    % Revised code for each slepian frunction
    [ESTfit,ESTdelta,ESTtotalparams_CC,ESTparamerrors_CC,ESTftests_CC] = timeseriesfit([ESTthedates' ESTtotal_CC(i,:)'],ESTalphavarall_CC(i),1,1);

    % find the appropriate polyfit number
    if any(max(find(ESTftests_CC==1)))
        ESTpoly_CC(i)=max(find(ESTftests_CC==1));
    else
        ESTpoly_CC(i)=1; %default is linear regression
    end

    ESTtotalparams_CC(:,i)=ESTtotalparams_CC(:,ESTpoly_CC(i));
    % Make a matrix for the line, and 95% confidence in the fit
    ESTtotalfit_CC(:,(3*i-2):3*i) = [ESTthedates' ESTfit(:,ESTpoly_CC(i)) ESTdelta(:,ESTpoly_CC(i))];

    % Make the error valid for a year
    ESTtotalparamerrors_CC(:,i) = ESTparamerrors_CC(:,ESTpoly_CC(i))*365.25;

    %     total_CC(i,:)=total;
    %     alphavarall_CC(i,:)=alphavarall;
    %     totalparams_CC(i,:)=totalparams';
    %     totalparamerrors_CC(i,:)=totalparamerrors';
    %     totalfit_CC(:,(3*i-2):3*i)=totalfit;

end



end
