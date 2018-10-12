function biomass = picrustDetermineBiomass(model,fluxdist,flexibleBiomass)

if ~exist('flexibleBiomass','var')
    flexibleBiomass = 0;
end
if flexibleBiomass
  biomass = sum(abs(fluxdist( strcmp(model.subSystems,'FLEX_BIOM') )));
else
    biomass = sum(abs(fluxdist( strcmp(model.subSystems,'Exchange') | strcmp(model.subSystems,'Transport') | strcmp(model.subSystems,'') )));
end
