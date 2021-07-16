function [model, deletedMets, deletedRxns] = Merge_mets_by_structure(model,table_similarity)

% Merge_mets_by_structure.m
% 
% Author: Francesco Balzerani, Telmo Blasco, Iñigo Apaolaza, Francisco J.
% Planes
% Email: fbalzerani@tecnun.es, tblasco@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

% Load the cpd IDs
ids = table_similarity(:,2);
n_table = size(ids,1);
to_del = cell(n_table,1);

for i = 1 : n_table
    cpd = regexp(ids{i},'cpd\d+','match');

    pos = find(ismember(model.metID,cpd));
    to_del{i,1} = pos(2:end);
    model.S(pos(1),:) = sum(model.S(pos,:));
end

% Concatenate all the position to remove the information
to_del = cat(1,to_del{:});

% Extract the information of the deleted metabolites
deletedMets = model.metID(to_del);

model.S(to_del,:) = [];
model.metAbbreviation(to_del) = [];
model.metNames(to_del) = [];
model.metcharge(to_del) = [];
model.metFormula(to_del) = [];
model.metSmiles(to_del) = [];
model.metID(to_del) = [];
model.metMass(to_del) = [];
model.metSource(to_del) = [];
model.metInchiKey(to_del) = [];
model.metIsCore(to_del) = [];
model.metIsObsolete(to_del) = [];
model.metLinkedCompound(to_del) = [];
model.metIsCofactor(to_del) = [];
model.metDeltaG(to_del) = [];
model.metDeltaGErr(to_del) = [];
model.metPKA(to_del) = [];
model.metPKB(to_del) = [];
model.metAbstractCompound(to_del) = [];
model.metComprisedOf(to_del) = [];
model.metAliases(to_del) = [];
model.metOntology(to_del) = [];

[~,~,ic] = unique(model.S','rows','stable');

[index_rxns,~] = hist(ic,unique(ic));
pos_rep = find(index_rxns>1);
n_rep = length(pos_rep);
rxns_to_remove = [];

mantainRxnID = cell(n_rep,1);
newRxnID = cell(n_rep,1);
for i = 1 : n_rep
    tmp_pos = find(ismember(ic,pos_rep(i)));
    pos_corr_rxns = tmp_pos(model.correctedRxns(tmp_pos)==0);
    if length(pos_corr_rxns) > 1
        notEmptyValues = ~cellfun(@isempty,model.rxnECnumbers(pos_corr_rxns));
        tmp = join(model.rxnECnumbers(pos_corr_rxns(notEmptyValues)),'|');
        model.rxnECnumbers(pos_corr_rxns(1)) = join(unique(split(tmp,'|')),'|');

        rxns_to_remove = [rxns_to_remove; pos_corr_rxns(2:end)];

        mantainRxnID(i) = model.rxnID(pos_corr_rxns(1));
        newRxnID{i} = model.rxnID(pos_corr_rxns(2:end));
    end
end

% Extract the information of the deleted reactions
deletedRxns = [mantainRxnID, newRxnID];

% Delete the information for the repeated reactions
model.S(:,rxns_to_remove) = [];
model.rxnAbbreviation(rxns_to_remove) = [];
model.lb(rxns_to_remove) = [];
model.ub(rxns_to_remove) = [];
model.rxnNames(rxns_to_remove) = [];
model.correctedRxns(rxns_to_remove) = [];
model.rxnECnumbers(rxns_to_remove) = [];
model.rxnID(rxns_to_remove) = [];
model.rxnCode(rxns_to_remove) = [];
model.rxnStoichiometry(rxns_to_remove) = [];
model.rxnIsTransport(rxns_to_remove) = [];
model.rxnEquation(rxns_to_remove) = [];
model.rxnDefinition(rxns_to_remove) = [];
model.rxnReversibility(rxns_to_remove) = [];
model.rxnDirection(rxns_to_remove) = [];
model.rxnAbstractReaction(rxns_to_remove) = [];
model.rxnPathways(rxns_to_remove) = [];
model.rxnAliases(rxns_to_remove) = [];
model.rxnDeltaG(rxns_to_remove) = [];
model.rxnDeltaGErr(rxns_to_remove) = [];
model.rxnCompoundIDs(rxns_to_remove) = [];
model.rxnStatus(rxns_to_remove) = [];
model.rxnIsObsolete(rxns_to_remove) = [];
model.rxnLinkedReaction(rxns_to_remove) = [];
model.rxnNotes(rxns_to_remove) = [];
model.rxnSource(rxns_to_remove) = [];
model.rxnOntology(rxns_to_remove) = [];
end