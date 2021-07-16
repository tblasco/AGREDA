function model = ConcatenateModels(model_agora, model_seed, so_model_name)

% ConcatenateModels.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

clc
disp('Concatenating AGORA and SEED Models');
model.modelName = so_model_name;

%% Merge all the fields of AGORA and SEED models
model.rxns = [model_agora.rxns;model_seed.rxnAbbreviation];
model.S = [model_agora.S , zeros(size(model_agora.S,1),size(model_seed.S,2))
    zeros(size(model_seed.S,1),size(model_agora.S,2)), model_seed.S];
model.lb = [model_agora.lb;model_seed.lb];
model.ub = [model_agora.ub;model_seed.ub];
model.c  = [model_agora.c;zeros(length(model_seed.rxnNames),1)];
model.mets  = [model_agora.mets;model_seed.metAbbreviation]; 
model.b  = [model_agora.b;zeros(length(model_seed.metNames),1)];
model.rules = [model_agora.rules;cell(length(model_seed.rxnNames),1)];
model.genes = model_agora.genes;
model.osense = -1;
model.csense  = [model_agora.csense; repmat('E', length(model_seed.metID),1)];
model.rxnGeneMat = [model_agora.rxnGeneMat ; sparse(length(model_seed.rxnNames), length(model_agora.genes))];
model.rxnNames = [model_agora.rxnNames;model_seed.rxnNames];
model.subSystems = [model_agora.subSystems;cell(length(model_seed.rxnNames),1)];
model.metNames = [model_agora.metNames;model_seed.metNames];
model.grRules = [model_agora.grRules;cell(length(model_seed.rxnNames),1)];
model.comments = [model_agora.comments;cell(length(model_seed.rxnNames),1)];
model.citations = [model_agora.citations;cell(length(model_seed.rxnNames),1)];
model.rxnConfidenceScores = [model_agora.rxnConfidenceScores;zeros(length(model_seed.rxnNames),1)];
model.rxnECNumbers = [model_agora.rxnECNumbers;model_seed.rxnECnumbers];
model.rxnKEGGID = [model_agora.rxnKEGGID;cell(length(model_seed.rxnNames),1)];
model.metCharges = [num2cell(model_agora.metCharges);model_seed.metcharge];
model.metFormulas = [model_agora.metFormulas;model_seed.metFormula];
model.metChEBIID = [model_agora.metChEBIID;cell(length(model_seed.metNames),1)];
model.metKEGGID = [model_agora.metKEGGID;cell(length(model_seed.metNames),1)];
model.metPubChemID = [model_agora.metPubChemID;cell(length(model_seed.metNames),1)];
model.metInChIString = [model_agora.metInChIString;cell(length(model_seed.metNames),1)];
model.metHMDBID = [model_agora.metHMDBID;cell(length(model_seed.metNames),1)];
model.metSmiles = [model_agora.metSmiles;model_seed.metSmiles];
model.description = 0;
model.disabled = cell(0);
model.modelID = model_agora.modelID;
model.taxonomy = model_agora.taxonomy;
model.metsTax = [model_agora.metsTax;cell(length(model_seed.metID),1)];
model.trRules = [model_agora.trRules; model_seed.trRules];
model.trules = [model_agora.trules; model_seed.trules];
model.rxnTaxMat = [model_agora.rxnTaxMat; model_seed.rxnTaxMat ];
model.correctedRxns = [zeros(length(model_agora.rxns),1);model_seed.correctedRxns];
model.metID = [model_agora.metID;model_seed.metID];
model.rxnID = [model_agora.rxnID;model_seed.rxnID];
model.rxnCode = [cell(length(model_agora.rxns),1);model_seed.rxnCode];
model.rxnStoichiometry = [cell(length(model_agora.rxns),1);model_seed.rxnStoichiometry];
model.rxnIsTransport = [cell(length(model_agora.rxns),1);model_seed.rxnIsTransport];
model.rxnEquation = [cell(length(model_agora.rxns),1);model_seed.rxnEquation];
model.rxnDefinition = [cell(length(model_agora.rxns),1);model_seed.rxnDefinition];
model.rxnReversibility = [cell(length(model_agora.rxns),1);model_seed.rxnReversibility];
model.rxnDirection = [cell(length(model_agora.rxns),1);model_seed.rxnDirection];
model.rxnAbstractReaction = [cell(length(model_agora.rxns),1);model_seed.rxnAbstractReaction];
model.rxnPathways = [cell(length(model_agora.rxns),1);model_seed.rxnPathways];
model.rxnAliases = [cell(length(model_agora.rxns),1);model_seed.rxnAliases];
model.rxnDeltaG = [cell(length(model_agora.rxns),1);model_seed.rxnDeltaG];
model.rxnDeltaGErr = [cell(length(model_agora.rxns),1);model_seed.rxnDeltaGErr];
model.rxnCompoundIDs = [cell(length(model_agora.rxns),1);model_seed.rxnCompoundIDs];
model.rxnStatus = [cell(length(model_agora.rxns),1);model_seed.rxnStatus];
model.rxnIsObsolete = [cell(length(model_agora.rxns),1);model_seed.rxnIsObsolete];
model.rxnLinkedReaction = [cell(length(model_agora.rxns),1);model_seed.rxnLinkedReaction];
model.rxnNotes = [cell(length(model_agora.rxns),1);model_seed.rxnNotes];
model.rxnSource = [cell(length(model_agora.rxns),1);model_seed.rxnSource];
model.rxnOntology = [cell(length(model_agora.rxns),1);model_seed.rxnOntology];
model.metMass = [cell(length(model_agora.mets),1);model_seed.metMass];
model.metSource = [cell(length(model_agora.mets),1);model_seed.metSource];
model.metInchiKey = [cell(length(model_agora.mets),1);model_seed.metInchiKey];
model.metIsCore = [cell(length(model_agora.mets),1);model_seed.metIsCore];
model.metIsObsolete = [cell(length(model_agora.mets),1);model_seed.metIsObsolete];
model.metLinkedCompound = [cell(length(model_agora.mets),1);model_seed.metLinkedCompound];
model.metIsCofactor = [cell(length(model_agora.mets),1);model_seed.metIsCofactor];
model.metDeltaG = [cell(length(model_agora.mets),1);model_seed.metDeltaG];
model.metDeltaGErr = [cell(length(model_agora.mets),1);model_seed.metDeltaGErr];
model.metPKA = [cell(length(model_agora.mets),1);model_seed.metMass];
model.metPKB = [cell(length(model_agora.mets),1);model_seed.metPKB];
model.metAbstractCompound = [cell(length(model_agora.mets),1);model_seed.metAbstractCompound];
model.metComprisedOf = [cell(length(model_agora.mets),1);model_seed.metComprisedOf];
model.metAliases = [cell(length(model_agora.mets),1);model_seed.metAliases];
model.metOntology = [cell(length(model_agora.mets),1);model_seed.metOntology];
end