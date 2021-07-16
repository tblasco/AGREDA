function [model,cpd_repeated,not_cyt_mets] = mergeSEEDmetabolites(model, comparison_table, so_model_name)

% mergeSEEDmetabolites.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J. Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es, fplanes@tecnun.es
% Date: 07/06/2019

clc
disp('Merging AGORA and SEED Models');
model.modelName = so_model_name;

%% Merge the metabolites

% Load AGORA and modelSEED metabolites information
cpd_id = table2cell(comparison_table(:,'SEED_ID'));
mets_agora = table2cell(comparison_table(:,'AGORA_METABOLITES'));

% Add the modelSEED's cpd id field to the AGORA's metabolites (for those that have)
n_agora = 1:length(mets_agora);
model.metID(n_agora(~cellfun(@isempty,cpd_id))) = cpd_id(~cellfun(@isempty,cpd_id));

% Identify the unique cpd id
cpd_no_rep = unique(cpd_id(~cellfun(@isempty,cpd_id)));
n_un_cpd = length(cpd_no_rep);

count = 1;
not_cyt_pos=1;
rxns_to_delete = [];
mets_to_delete = [];
for i = 1:n_un_cpd
    clc
    disp('Merging AGORA and SEED Models');
    disp('Merging the metabolites in common');
    disp(['SEED ID: ' num2str(i) '/' num2str(n_un_cpd)]);
    
% Identify the metabolites relating to the i-th cpd id
    mask_common_id = strcmp(cpd_no_rep(i), model.metID);
    index_common_mets = find(mask_common_id);

% Link AGORA and SEED metabolites
    if index_common_mets(end) > length(mets_agora)

        cyt_pos = find(~cellfun(@isempty, regexp(model.mets(index_common_mets),'\[c\]'))); % check position of cytosolic met

        if length(cyt_pos) > 1
            cpd_repeated{count,1} = model.mets(index_common_mets);
            count = count + 1;
            continue

        elseif isempty(cyt_pos)
            not_cyt_mets(not_cyt_pos,1) = model.metNames(index_common_mets(1));
            not_cyt_pos = not_cyt_pos+1;
            cyt_pos = find(~cellfun(@isempty, regexp(model.mets(index_common_mets),'\[e\]')));
        end           

        rxns_pos = find(model.S(index_common_mets(end),:)~=0);
        all_rxns_met = model.S(:,rxns_pos);
        [~,rxns_filled] = find(all_rxns_met);
        [met_counts_x_rxn,~] = hist(rxns_filled,unique(rxns_filled));
        rxns_to_delete = [rxns_to_delete, rxns_pos(find(met_counts_x_rxn==1))];

        mets_to_delete = [mets_to_delete, index_common_mets(end)];
     
% Add SEED's reaction information to AGORA metabolites
        model.S(index_common_mets(cyt_pos), :) = model.S(index_common_mets(cyt_pos),:) + model.S(index_common_mets(end), :);
        model.metMass(index_common_mets(cyt_pos), :) = model.metMass(index_common_mets(end), :);
        model.metSource(index_common_mets(cyt_pos), :) = model.metSource(index_common_mets(end), :);
        model.metInchiKey(index_common_mets(cyt_pos), :) = model.metInchiKey(index_common_mets(end), :);
        model.metIsCore(index_common_mets(cyt_pos), :) = model.metIsCore(index_common_mets(end), :);
        model.metIsObsolete(index_common_mets(cyt_pos), :) = model.metIsObsolete(index_common_mets(end), :);
        model.metLinkedCompound(index_common_mets(cyt_pos), :) = model.metLinkedCompound(index_common_mets(end), :);
        model.metIsCofactor(index_common_mets(cyt_pos), :) = model.metIsCofactor(index_common_mets(end), :);
        model.metDeltaG(index_common_mets(cyt_pos), :) = model.metDeltaG(index_common_mets(end), :);
        model.metDeltaGErr(index_common_mets(cyt_pos), :) = model.metDeltaGErr(index_common_mets(end), :);
        model.metPKA(index_common_mets(cyt_pos), :) = model.metPKA(index_common_mets(end), :);
        model.metPKB(index_common_mets(cyt_pos), :) = model.metPKB(index_common_mets(end), :);
        model.metAbstractCompound(index_common_mets(cyt_pos), :) = model.metAbstractCompound(index_common_mets(end), :);
        model.metComprisedOf(index_common_mets(cyt_pos), :) = model.metComprisedOf(index_common_mets(end), :);
        model.metAliases(index_common_mets(cyt_pos), :) = model.metAliases(index_common_mets(end), :);
        model.metOntology(index_common_mets(cyt_pos), :) = model.metOntology(index_common_mets(end), :);
          
    end
end

%Delete metabolite information
model.S(mets_to_delete,:) = [];
model.mets(mets_to_delete) = [];
model.b(mets_to_delete) = [];
model.csense(mets_to_delete) = [];
model.metNames(mets_to_delete) = [];
model.metCharges(mets_to_delete) = [];
model.metFormulas(mets_to_delete) = [];
model.metChEBIID(mets_to_delete) = [];
model.metKEGGID(mets_to_delete) = [];
model.metPubChemID(mets_to_delete) = [];
model.metInChIString(mets_to_delete) = [];
model.metHMDBID(mets_to_delete) = [];
model.metSmiles(mets_to_delete) = [];
model.metID(mets_to_delete) = [];
model.metsTax(mets_to_delete) = [];
model.metMass(mets_to_delete) = [];
model.metSource(mets_to_delete) = [];
model.metInchiKey(mets_to_delete) = [];
model.metIsCore(mets_to_delete) = [];
model.metIsObsolete(mets_to_delete) = [];
model.metLinkedCompound(mets_to_delete) = [];
model.metIsCofactor(mets_to_delete) = [];
model.metDeltaG(mets_to_delete) = [];
model.metDeltaGErr(mets_to_delete) = [];
model.metPKA(mets_to_delete) = [];
model.metPKB(mets_to_delete) = [];
model.metAbstractCompound(mets_to_delete) = [];
model.metComprisedOf(mets_to_delete) = [];
model.metAliases(mets_to_delete) = [];
model.metOntology(mets_to_delete) = [];

%Delete exchange reaction information
model.S(:,rxns_to_delete) = [];
model.rxns(rxns_to_delete) = [];
model.lb(rxns_to_delete) = [];
model.ub(rxns_to_delete) = [];
model.c(rxns_to_delete) = [];
model.rules(rxns_to_delete) = [];
model.rxnNames(rxns_to_delete) = [];
model.subSystems(rxns_to_delete) = [];
model.grRules(rxns_to_delete) = [];
model.comments(rxns_to_delete) = [];
model.citations(rxns_to_delete) = [];
model.rxnConfidenceScores(rxns_to_delete) = [];
model.correctedRxns(rxns_to_delete) = [];
model.rxnECNumbers(rxns_to_delete) = [];
model.rxnKEGGID(rxns_to_delete) = [];
model.rxnID(rxns_to_delete) = [];
model.rxnCode(rxns_to_delete) = [];
model.rxnStoichiometry(rxns_to_delete) = [];
model.rxnIsTransport(rxns_to_delete) = [];
model.rxnEquation(rxns_to_delete) = [];
model.rxnDefinition(rxns_to_delete) = [];
model.rxnReversibility(rxns_to_delete) = [];
model.rxnDirection(rxns_to_delete) = [];
model.rxnAbstractReaction(rxns_to_delete) = [];
model.rxnPathways(rxns_to_delete) = [];
model.rxnAliases(rxns_to_delete) = [];
model.rxnDeltaG(rxns_to_delete) = [];
model.rxnDeltaGErr(rxns_to_delete) = [];
model.rxnCompoundIDs(rxns_to_delete) = [];
model.rxnStatus(rxns_to_delete) = [];
model.rxnIsObsolete(rxns_to_delete) = [];
model.rxnLinkedReaction(rxns_to_delete) = [];
model.rxnNotes(rxns_to_delete) = [];
model.rxnSource(rxns_to_delete) = [];
model.rxnOntology(rxns_to_delete) = [];
model.trules(rxns_to_delete) = [];
model.trRules(rxns_to_delete) = [];
model.rxnGeneMat(rxns_to_delete,:) = [];
model.rxnTaxMat(rxns_to_delete,:) = [];

% Delete repeated reactions generated when merging metabolites
clc
disp('Merging AGORA and SEED Models');
disp('Removing the repeated reactions');

[~,~,ic] = unique(model.S','rows','stable');

[index_rxns,~] = hist(ic,unique(ic));
pos_rep = find(index_rxns>1);
n_rep = length(pos_rep);
rxns_to_remove = [];

for i = 1 : n_rep
    clc
    disp('Merging AGORA and SEED Models');
    disp('Removing the repeated reactions');
    disp(['SEED ID: ' num2str(i) '/' num2str(n_rep)]);
    
    tmp_pos = find(ismember(ic,pos_rep(i)));
    pos_corr_rxns = tmp_pos(model.correctedRxns(tmp_pos)==0);
    if length(pos_corr_rxns) > 1
        model.rxnCode(pos_corr_rxns(1)) = model.rxnCode(pos_corr_rxns(2));
        model.rxnStoichiometry(pos_corr_rxns(1)) = model.rxnStoichiometry(pos_corr_rxns(2));
        model.rxnIsTransport(pos_corr_rxns(1)) = model.rxnIsTransport(pos_corr_rxns(2));
        model.rxnEquation(pos_corr_rxns(1)) = model.rxnEquation(pos_corr_rxns(2));
        model.rxnDefinition(pos_corr_rxns(1)) = model.rxnDefinition(pos_corr_rxns(2));
        model.rxnReversibility(pos_corr_rxns(1)) = model.rxnReversibility(pos_corr_rxns(2));
        model.rxnDirection(pos_corr_rxns(1)) = model.rxnDirection(pos_corr_rxns(2));
        model.rxnAbstractReaction(pos_corr_rxns(1)) = model.rxnAbstractReaction(pos_corr_rxns(2));
        model.rxnPathways(pos_corr_rxns(1)) = model.rxnPathways(pos_corr_rxns(2));
        model.rxnAliases(pos_corr_rxns(1)) = model.rxnAliases(pos_corr_rxns(2));
        model.rxnDeltaG(pos_corr_rxns(1)) = model.rxnDeltaG(pos_corr_rxns(2));
        model.rxnDeltaGErr(pos_corr_rxns(1)) = model.rxnDeltaGErr(pos_corr_rxns(2));
        model.rxnCompoundIDs(pos_corr_rxns(1)) = model.rxnCompoundIDs(pos_corr_rxns(2));
        model.rxnStatus(pos_corr_rxns(1)) = model.rxnStatus(pos_corr_rxns(2));
        model.rxnIsObsolete(pos_corr_rxns(1)) = model.rxnIsObsolete(pos_corr_rxns(2));
        model.rxnLinkedReaction(pos_corr_rxns(1)) = model.rxnLinkedReaction(pos_corr_rxns(2));
        model.rxnNotes(pos_corr_rxns(1)) = model.rxnNotes(pos_corr_rxns(2));
        model.rxnSource(pos_corr_rxns(1)) = model.rxnSource(pos_corr_rxns(2));
        model.rxnOntology(pos_corr_rxns(1)) = model.rxnOntology(pos_corr_rxns(2));
        
        rxns_to_remove = [rxns_to_remove; pos_corr_rxns(2:end)];
        
%       Join taxonomic information of reactions
        
        pos_trules_not_empty = find(~cellfun(@isempty,model.trules(pos_corr_rxns)));

        tmp_join_trules = strjoin(model.trules(pos_corr_rxns(pos_trules_not_empty)),' | ');
        tmp_split_trules = strsplit(tmp_join_trules,' | ');
        model.trules{pos_corr_rxns(1)} = strjoin(unique(tmp_split_trules),' | ');
        tmp_join_trRules = strjoin(model.trRules(pos_corr_rxns(pos_trules_not_empty)),' or ');
        tmp_split_trRules = strsplit(tmp_join_trRules,' or ');
        model.trRules{pos_corr_rxns(1)} = strjoin(unique(tmp_split_trRules),' or ');
        model.rxnTaxMat(pos_corr_rxns(1),:) = sum(model.rxnTaxMat(pos_corr_rxns,:));
    end
end


%Delete exchange reaction information
model.S(:,rxns_to_remove) = [];
model.rxns(rxns_to_remove) = [];
model.lb(rxns_to_remove) = [];
model.ub(rxns_to_remove) = [];
model.c(rxns_to_remove) = [];
model.rules(rxns_to_remove) = [];
model.rxnNames(rxns_to_remove) = [];
model.subSystems(rxns_to_remove) = [];
model.grRules(rxns_to_remove) = [];
model.comments(rxns_to_remove) = [];
model.citations(rxns_to_remove) = [];
model.rxnConfidenceScores(rxns_to_remove) = [];
model.correctedRxns(rxns_to_remove) = [];
model.rxnECNumbers(rxns_to_remove) = [];
model.rxnKEGGID(rxns_to_remove) = [];
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
model.trules(rxns_to_remove) = [];
model.trRules(rxns_to_remove) = [];
model.rxnGeneMat(rxns_to_remove,:) = [];
model.rxnTaxMat(rxns_to_remove,:) = [];
model.rxnTaxMat(model.rxnTaxMat>1) = 1;


end