function model = mergeDupRxns(model)

[~,ia] = unique(model.S','rows','stable');
D = setdiff(1:size(model.S',1),ia);
n = length(D);
for k = 1:n
    clc
    fprintf('Merging duplicated reactions. Reaction %d/%d\n',k,n)
    for i = 1:(D(k)-1)
        if isequal(model.S(:,i),model.S(:,D(k)))
            if ~isequal(model.ub(i),model.ub(D(k))) && model.ub(i)==0
                model.ub(i) = model.ub(D(k));
            end
            if isempty(model.trules{i})
                model.trules(i) = model.trules(D(k));
                tmp = split(model.trRules(D(k)),' or ');
                model.trRules(i) = model.trRules(D(k));
                model.rxnTaxMat(i,ismember(model.taxonomy, tmp)) = 1;
            else
                model.trules(i) = join(unique(split(join([model.trules(i), model.trules(D(k))],' | '),' | ')),' | ');
                tmp = unique(split(join([model.trRules(i), model.trRules(D(k))],' or '),' or '));
                model.trRules(i) = join(tmp,' or ');
                model.rxnTaxMat(i,:) = zeros(1, length(model.taxonomy));
                model.rxnTaxMat(i,ismember(model.taxonomy, tmp)) = 1;
            end
            break;
        end
    end
end

model.S(:,D) = [];
model.rxns(D) = [];
model.lb(D) = [];
model.ub(D) = [];
model.c(D) = [];
model.rules(D) = [];
model.rxnNames(D) = [];
model.subSystems(D) = [];
model.grRules(D) = [];
model.comments(D) = [];
model.citations(D) = [];
model.rxnConfidenceScores(D) = [];
model.correctedRxns(D) = [];
model.rxnECNumbers(D) = [];
model.rxnKEGGID(D) = [];
model.rxnID(D) = [];
model.rxnCode(D) = [];
model.rxnStoichiometry(D) = [];
model.rxnIsTransport(D) = [];
model.rxnEquation(D) = [];
model.rxnDefinition(D) = [];
model.rxnReversibility(D) = [];
model.rxnDirection(D) = [];
model.rxnAbstractReaction(D) = [];
model.rxnPathways(D) = [];
model.rxnAliases(D) = [];
model.rxnDeltaG(D) = [];
model.rxnDeltaGErr(D) = [];
model.rxnCompoundIDs(D) = [];
model.rxnStatus(D) = [];
model.rxnIsObsolete(D) = [];
model.rxnLinkedReaction(D) = [];
model.rxnNotes(D) = [];
model.rxnSource(D) = [];
model.rxnOntology(D) = [];
model.trules(D) = [];
model.trRules(D) = [];
model.rxnGeneMat(D,:) = [];
model.rxnTaxMat(D,:) = [];

end