function model = mergeSupraOrganism(model, so_model_name)

% mergeSupraOrganism.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

clc
disp('Calculating the Supra Organism Model');
model.modelName = so_model_name;

%% Merge the metabolites
un_mets_rep = unique(model.mets);
n_un_mets = length(un_mets_rep);
for i = 1:n_un_mets
    clc
    disp('Calculating the Supra Organism Model');
    disp(['Metabolites: ' num2str(i) '/' num2str(n_un_mets)]);
    mask_common_mets = strcmp(un_mets_rep(i), model.mets);
    if sum(mask_common_mets) > 1
        index_common_mets = find(mask_common_mets);
        model.S(index_common_mets(1), :) = sum(model.S(index_common_mets, :));
        
        model.metsTax{index_common_mets(1)} = strjoin(model.metsTax(index_common_mets), ' or ');
        
        model.S(index_common_mets(2:end), :) = [];
        model.mets(index_common_mets(2:end), :) = [];
        model.b(index_common_mets(2:end), :) = [];
        model.csense(index_common_mets(2:end), :) = [];
        model.metNames(index_common_mets(2:end), :) = [];
        model.metCharges(index_common_mets(2:end), :) = [];
        model.metFormulas(index_common_mets(2:end), :) = [];
        model.metChEBIID(index_common_mets(2:end), :) = [];
        model.metKEGGID(index_common_mets(2:end), :) = [];
        model.metPubChemID(index_common_mets(2:end), :) = [];
        model.metInChIString(index_common_mets(2:end), :) = [];
        model.metHMDBID(index_common_mets(2:end), :) = [];
        model.metSmiles(index_common_mets(2:end), :) = [];
        model.metsTax(index_common_mets(2:end), :) = [];
    end
end

%% Merge the reactions
i = 0;
stop = 0;
rxns_to_analyze = 1:length(model.rxns);
rxns_to_delete = [];
while stop == 0
    i = i+1;
    clc
    disp('Calculating the Supra Organism Model');
    disp('Metabolites: DONE!');
    disp(['Reactions: ' num2str(i) '/' num2str(length(rxns_to_analyze))]);
    act_pos = rxns_to_analyze(i);
    act_rxn = model.S(:, act_pos);
    same_rxn = [];
    expanded_act_rxn = repmat(act_rxn, 1, size(model.S, 2));
    same_rxn = find(sum((model.S-expanded_act_rxn)~=0, 1) == 0);
    same_rxn = same_rxn(same_rxn ~= act_pos);
    
    if ~isempty(same_rxn)
        rxns_to_analyze(ismember(rxns_to_analyze, same_rxn)) = [];
        rxns_to_delete = [rxns_to_delete, same_rxn];

        model.ub(act_pos) = model.ub(act_pos)+sum(model.ub(same_rxn));
        model.lb(act_pos) = model.lb(act_pos)+sum(model.lb(same_rxn));

        model.trRules{act_pos} = strjoin(unique(model.trRules([act_pos, same_rxn])), ' or ');
        model.trules{act_pos} = strjoin(unique(model.trules([act_pos, same_rxn])), ' | ');    
    end
    if i == length(rxns_to_analyze)
        stop = 1;
    end    
end
model.rxns(rxns_to_delete) = [];
model.S(:, rxns_to_delete) = [];
model.lb(rxns_to_delete) = [];
model.ub(rxns_to_delete) = [];
model.c(rxns_to_delete) = [];
model.rules(rxns_to_delete) = [];
model.rxnGeneMat(rxns_to_delete, :) = [];
model.rxnNames(rxns_to_delete) = [];
model.subSystems(rxns_to_delete) = [];
model.grRules(rxns_to_delete) = [];
model.comments(rxns_to_delete) = [];
model.citations(rxns_to_delete) = [];
model.rxnConfidenceScores(rxns_to_delete) = [];
model.rxnECNumbers(rxns_to_delete) = [];
model.rxnKEGGID(rxns_to_delete) = []; 
model.trRules(rxns_to_delete) = [];
model.trules(rxns_to_delete) = [];

model.rxnTaxMat = sparse(length(model.rxns), length(model.taxonomy));
pos = cellfun(@regexp, model.trules, repmat({'\d*'}, length(model.trules), 1), repmat({'Match'}, length(model.trules), 1), 'UniformOutput', false);
pos = cellfun(@str2double, pos, 'UniformOutput', false);
n_rxns = length(model.rxns);
for i = 1:n_rxns
    model.rxnTaxMat(i, pos{i}) = 1;
end
end