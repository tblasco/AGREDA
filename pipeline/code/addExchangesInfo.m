function model = addExchangesInfo(model,id)

% addExchangesInfo.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

% Create overall reaction information
tmp = cellfun(@(x) find(ismember(model.metID,x)), id, 'UniformOutput', false);
tmp(ismember(id,'cpd01202')) = {2391};
final_idx = cell2mat(tmp);

rxnID = cellfun(@strcat, repmat({'rxnTransport'},length(final_idx),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(final_idx))'),'UniformOutput',false),'UniformOutput',false);
rxns = rxnID;
rxnNames = rxnID;
trRules = repmat(join(model.taxonomy, ' or '), length(final_idx), 1);
tmp = join(join([repmat({'t('}, 818, 1), cellfun(@num2str, num2cell((1:818)'), 'UniformOutput', false),repmat({')'}, 818, 1)],''), ' | ');
trules = repmat(tmp, length(final_idx),1);

info = [rxnID rxns rxnNames trules trRules, model.mets(final_idx)];


n = size(info,1);
origRxns = length(model.rxns);

% Add relevant field info
model.S = [model.S, zeros(length(model.mets), n)];
for i = 1:n
    model.S(ismember(model.mets,info(i,6)),origRxns+i) = 1;
end
model.rxnID = [model.rxnID; info(:,1)];
model.rxns = [model.rxns; info(:,2)];
model.rxnNames = [model.rxnNames; info(:,3)];
model.ub = [model.ub; repmat(1000, n, 1)];
model.lb = [model.lb; repmat(-1000, n, 1)];
model.trules = [model.trules; info(:,4)];
model.trRules = [model.trRules; info(:,5)];
model.rxnTaxMat = [model.rxnTaxMat; zeros(n, size(model.rxnTaxMat,2))];
for i = 1:n
    idx = find(ismember(model.taxonomy, strsplit(info{i,5},' or ')));
    model.rxnTaxMat(origRxns+i,idx) = 1;
end

% Add rest field info
model.rules = [model.rules;cell(n,1)];
model.grRules = [model.grRules;cell(n,1)];
model.subSystems = [model.subSystems;cell(n,1)];
model.comments = [model.comments;cell(n,1)];
model.citations = [model.citations;cell(n,1)];
model.rxnECNumbers = [model.rxnECNumbers;cell(n,1)];
model.rxnKEGGID = [model.rxnKEGGID;cell(n,1)];
model.rxnOntology = [model.rxnOntology;cell(n,1)];
model.rxnSource = [model.rxnSource;cell(n,1)];
model.rxnNotes = [model.rxnNotes;cell(n,1)];
model.rxnStatus = [model.rxnStatus;cell(n,1)];
model.rxnLinkedReaction = [model.rxnLinkedReaction;cell(n,1)];
model.rxnIsObsolete = [model.rxnIsObsolete;cell(n,1)];
model.rxnCompoundIDs = [model.rxnCompoundIDs;cell(n,1)];
model.rxnDeltaG = [model.rxnDeltaG;cell(n,1)];
model.rxnDeltaGErr = [model.rxnDeltaGErr;cell(n,1)];
model.rxnAliases = [model.rxnAliases;cell(n,1)];
model.rxnPathways = [model.rxnPathways;cell(n,1)];
model.rxnAbstractReaction = [model.rxnAbstractReaction;cell(n,1)];
model.rxnDirection = [model.rxnDirection;cell(n,1)];
model.rxnReversibility = [model.rxnReversibility;cell(n,1)];
model.rxnEquation = [model.rxnEquation;cell(n,1)];
model.rxnDefinition = [model.rxnDefinition;cell(n,1)];
model.rxnIsTransport = [model.rxnIsTransport;cell(n,1)];
model.rxnStoichiometry = [model.rxnStoichiometry;cell(n,1)];
model.rxnCode = [model.rxnCode;cell(n,1)];
model.c = [model.c; zeros(n,1)];
model.rxnConfidenceScores = [model.rxnConfidenceScores; zeros(n,1)];
model.correctedRxns = [model.correctedRxns; zeros(n,1)];
model.rxnGeneMat = [model.rxnGeneMat; zeros(n, size(model.rxnGeneMat,2))];

end