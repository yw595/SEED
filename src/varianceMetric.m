function metric = varianceMetric(normalfluxes1,obesefluxes1)

%numer = (sum(normalfluxes1.^2) - sum(normalfluxes1)^2/length(normalfluxes1) + sum(obesefluxes1.^2) - sum(obesefluxes1)^2/length(obesefluxes1));
%denom = (sum(normalfluxes1.^2) + sum(obesefluxes1.^2) - (sum(normalfluxes1)+sum(obesefluxes1))^2/(length(normalfluxes1)+length(obesefluxes1)));
numer = var(normalfluxes1)+var(obesefluxes1);
denom = var([normalfluxes1; obesefluxes1]);
if denom==0
    metric = NaN;
else
    metric = numer/denom;
end
