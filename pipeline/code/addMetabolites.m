function model = addMetabolites(model,metabolites)

% addMetabolites.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

n_mets = size(metabolites,1);

model.metID = [model.metID; table2cell(metabolites(:,'metID'))];
model.mets  = [model.mets; table2cell(metabolites(:,'mets'))];
model.metNames = [model.metNames; table2cell(metabolites(:,'metNames'))];
model.metFormulas = [model.metFormulas; table2cell(metabolites(:,'metFormulas'))];
model.metSmiles = [model.metSmiles; table2cell(metabolites(:,'metSmiles'))];

model.b  = [model.b;zeros(n_mets,1)];
model.csense  = [model.csense; repmat('E', n_mets,1)];
model.metCharges = [model.metCharges;cell(n_mets,1)];
model.metChEBIID = [model.metChEBIID;cell(n_mets,1)];
model.metKEGGID = [model.metKEGGID;cell(n_mets,1)];
model.metPubChemID = [model.metPubChemID;cell(n_mets,1)];
model.metInChIString = [model.metInChIString;cell(n_mets,1)];
model.metHMDBID = [model.metHMDBID;cell(n_mets,1)];
model.metsTax = [model.metsTax;cell(n_mets,1)];
model.metMass = [model.metMass;cell(n_mets,1)];
model.metSource = [model.metSource;cell(n_mets,1)];
model.metInchiKey = [model.metInchiKey;cell(n_mets,1)];
model.metIsCore = [model.metIsCore;cell(n_mets,1)];
model.metIsObsolete = [model.metIsObsolete;cell(n_mets,1)];
model.metLinkedCompound = [model.metLinkedCompound;cell(n_mets,1)];
model.metIsCofactor = [model.metIsCofactor;cell(n_mets,1)];
model.metDeltaG = [model.metDeltaG;cell(n_mets,1)];
model.metDeltaGErr = [model.metDeltaGErr;cell(n_mets,1)];
model.metPKA = [model.metPKA;cell(n_mets,1)];
model.metPKB = [model.metPKB;cell(n_mets,1)];
model.metAbstractCompound = [model.metAbstractCompound;cell(n_mets,1)];
model.metComprisedOf = [model.metComprisedOf;cell(n_mets,1)];
model.metAliases = [model.metAliases;cell(n_mets,1)];
model.metOntology = [model.metOntology;cell(n_mets,1)];

% Create the exchange reactions.
id = cellfun(@strcat, repmat({'rxnAddEX'},n_mets,1), cellfun(@(x) num2str(x), ...
    num2cell((1:n_mets)'),'UniformOutput',false),'UniformOutput',false);

model.rxnID = [model.rxnID;id];
model.rxns = [model.rxns;strcat('EX_',table2cell(metabolites(:,'metNames')))];
model.rxnNames = [model.rxnNames;strcat('EX_',table2cell(metabolites(:,'metNames')))];
model.S = [[model.S;zeros(n_mets,size(model.S,2))], [zeros(size(model.S,1),n_mets);-1*eye(n_mets)]];
model.lb = [model.lb;repmat(-1000,n_mets,1)];
model.ub = [model.ub;repmat(1000,n_mets,1)];

model.c  = [model.c;zeros(n_mets,1)];
model.rxnGeneMat = [model.rxnGeneMat ; sparse(n_mets, length(model.genes))];
model.subSystems = [model.subSystems;cell(n_mets,1)];
model.rules = [model.rules;cell(n_mets,1)];
model.grRules = [model.grRules;cell(n_mets,1)];
model.rxnConfidenceScores = [model.rxnConfidenceScores;zeros(n_mets,1)];
model.rxnECNumbers = [model.rxnECNumbers;cell(n_mets,1)];
model.rxnKEGGID = [model.rxnKEGGID;cell(n_mets,1)];
model.correctedRxns = [model.correctedRxns;zeros(n_mets,1)];
model.rxnCode = [model.rxnCode;cell(n_mets,1)];
model.rxnStoichiometry = [model.rxnStoichiometry;cell(n_mets,1)];
model.rxnIsTransport = [model.rxnIsTransport;cell(n_mets,1)];
model.rxnEquation = [model.rxnEquation;cell(n_mets,1)];
model.rxnDefinition = [model.rxnDefinition;cell(n_mets,1)];
model.rxnReversibility = [model.rxnReversibility;cell(n_mets,1)];
model.rxnDirection = [model.rxnDirection;cell(n_mets,1)];
model.rxnAbstractReaction = [model.rxnAbstractReaction;cell(n_mets,1)];
model.rxnPathways = [model.rxnPathways;cell(n_mets,1)];
model.rxnAliases = [model.rxnAliases;cell(n_mets,1)];
model.rxnDeltaG = [model.rxnDeltaG;cell(n_mets,1)];
model.rxnDeltaGErr = [model.rxnDeltaGErr;cell(n_mets,1)];
model.rxnCompoundIDs = [model.rxnCompoundIDs;cell(n_mets,1)];
model.rxnStatus = [model.rxnStatus;cell(n_mets,1)];
model.rxnIsObsolete = [model.rxnIsObsolete;cell(n_mets,1)];
model.rxnLinkedReaction = [model.rxnLinkedReaction;cell(n_mets,1)];
model.rxnNotes = [model.rxnNotes;cell(n_mets,1)];
model.rxnSource = [model.rxnSource;cell(n_mets,1)];
model.rxnOntology = [model.rxnOntology;cell(n_mets,1)];
model.comments = [model.comments;cell(n_mets,1)];
model.citations = [model.citations;cell(n_mets,1)];
model.trRules = [model.trRules; cell(n_mets,1)];
model.trules = [model.trules; cell(n_mets,1)];
model.rxnTaxMat = [model.rxnTaxMat; sparse(n_mets,length(model.taxonomy))];

end