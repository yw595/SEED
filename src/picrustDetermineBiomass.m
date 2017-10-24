function biomass = picrustDetermineBiomass(model,fluxdist)

  biomass = sum(abs(fluxdist( strcmp(model.subSystems,'Exchange') | strcmp(model.subSystems,'Transport') | strcmp(model.subSystems,'') )));
