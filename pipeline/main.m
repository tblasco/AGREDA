% main.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

close all
clear all
clc

addpath(fullfile('.', 'code'));
initCobraToolbox(false)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Supra-organism model building instructions
species_name_file = 'all_models';
species_path = fullfile('.', 'input', 'AGORA_models');

%AGREDA model building instructions
ANNOTATE = true;
agora_name_file = 'supraOrganism-all_models';
seed_name_file = 'modelSEED';
peg_path = fullfile('.', 'input', 'EC_Numbers', 'myRAST');
kegg_path = fullfile('.', 'input', 'EC_Numbers', 'KEGG');
fastFVA_path = fullfile('.', 'input', 'OptionsFastcore', 'removeNotAnnotatedRxnsFastFVA');

input_metabolites = readtable(fullfile('.', 'input', 'Annotation', 'Input_metabolites.xlsx'));
comparison_table = readtable(fullfile('.', 'input', 'Annotation', 'AGORA_SEED_with_exc.xlsx'));
manual_added_ECnumbers_to_rxns = readtable(fullfile('.', 'input', 'Annotation', 'manual_added_ECnumbers_to_rxns.xlsx'));
manual_added_trRules_to_ECnumbers = readtable(fullfile('.', 'input', 'Annotation', 'manual_added_trRules_to_ECnumbers.xlsx'));
manual_added_trRules_to_rxns = readtable(fullfile('.', 'input', 'Annotation', 'manual_added_trRules_to_rxns.xlsx'));
new_metabolites = readtable(fullfile('.', 'input', 'Annotation', 'new_metabolites.xlsx'));
new_reactions = readtable(fullfile('.', 'input', 'Annotation', 'new_reactions.xlsx'));
plants = readtable(fullfile('.', 'input', 'Annotation', 'no_bacteria.xlsx'));
similarity_table = table2cell(readtable(fullfile('.', 'input', 'Annotation', 'table_similarity.xlsx')));
id_species = readtable(fullfile('.', 'input', 'Annotation', 'speciesInfo.xlsx'));

removeMets = readtable(fullfile('.','input','Annotation','removeMetabolitesSEED.xlsx'));
mergeMetsSEED = readtable(fullfile('.','input','Annotation','mergeMetabolitesSEED.xlsx'));
exMetsWithEvidence = readtable(fullfile('.','input','Annotation','targetMetabolites.xlsx'));
FVA_balancing_path = fullfile('.','input','FVA_Species');

agora_exc_model_name = [agora_name_file '-with_exc'];
seed_tax_model_name = [seed_name_file '-with_tax'];
seed_tax_exc_model_name = [seed_name_file '-with_tax_exc'];
conc_model_name = 'agora_seed_concatenated';
merg_model_name = 'agora_seed_merged_expert';
minim_model_name = 'AGREDA_UNBALANCED';
final_model_name = 'AGREDA';

options.modelName = minim_model_name;
options.fastFVAName = fullfile('.', 'input', 'OptionsFastcore', ['fastFVAAgoraSeedExpertBounded_', minim_model_name]);
options.saving = 1;
if exist(fullfile('.','input', 'OptionsFastcore', ['blockedDietRxnsUp_', minim_model_name '.mat']),'file')
    load(fullfile('.','input', 'OptionsFastcore', ['blockedDietRxnsUp_', minim_model_name '.mat']))
    options.blockedDietRxnsUp = blockedDietRxnsUp;
end
if exist(fullfile('.','input', 'OptionsFastcore', ['blockedDietRxnsOut_', minim_model_name '.mat']),'file')
    load(fullfile('.','input', 'OptionsFastcore', ['blockedDietRxnsOut_', minim_model_name '.mat']))
    options.blockedDietRxnsOut = blockedDietRxnsOut;
end
options.plantRxns = table2cell(plants(:,1));

load(fullfile('.', 'input', 'speciesToMerge', [species_name_file '.mat']));

sl_model_name = ['speciesLevel-' species_name_file];
so_model_name = ['supraOrganism-' species_name_file];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Section 1: AGORA mixed-bag model
tic
if ~exist(fullfile('.', 'output', 'Models', [sl_model_name '.mat']))
    %Create the species-level model
    [model, not_included_species] = mergeSpeciesLevel(species_path, species, sl_model_name);
    if ~isempty(not_included_species)
        disp('Some of the species in ''species'' have not been included in the model. Check ''not_included_species''');
        keyboard;
    end
    save(fullfile('.', 'output', 'Models', [sl_model_name '.mat']), 'model', '-v7.3');
    
    %Create the supra-organism model
    model = mergeSupraOrganism(model, so_model_name);
    save(fullfile('.', 'output', 'Models', [so_model_name '.mat']), 'model');

else
    if ~exist(fullfile('.', 'output', 'Models', [so_model_name '.mat']))
        clc
        disp('Loading the Species Level Model');
        load(fullfile('.', 'output', 'Models', [sl_model_name '.mat']));
        
        %Create the supra-organism model
        model = mergeSupraOrganism(model, so_model_name);
        save(fullfile('.', 'output', 'Models', [so_model_name '.mat']), 'model');                
    end
end
toc
%% Section 2: AGREDA building
tic
load(fullfile('.', 'output', 'Models', [so_model_name '.mat']))

% Create i-Diet exchange reactions in the AGORA model
model_agora = Create_exchange_reaction_agora(model, input_metabolites, agora_exc_model_name);

if ~exist(fullfile('.', 'output', 'Models', [merg_model_name '.mat']))
    if ~exist(fullfile('.', 'output', 'Models', [seed_tax_model_name '.mat']))
        clc
        disp('Loading modelSEED')
        
        %Load modelSEED
        load(fullfile('.', 'output', 'Models', [seed_name_file '.mat']));

        %Merge modelSEED common metabolites
        model_seed = Merge_mets_by_structure(model,similarity_table);
        
        %Add taxonomy information
        [model_seed, species_not_included] = manageEC_pipeline(model_seed, species, peg_path, kegg_path, manual_added_ECnumbers_to_rxns, manual_added_trRules_to_ECnumbers, manual_added_trRules_to_rxns, seed_tax_model_name);
        if ~isempty(species_not_included)
            disp('For some of the species the peg file have not been downloaded. Check ''species_not_included''');
            keyboard;
        end
        save(fullfile('.', 'output', 'Models', [seed_tax_model_name '.mat']), 'model_seed');

        % Create exchange reactions in the SEED model for those metabolites
        % of i-diet that don't have one
        model_seed = Create_exchange_reaction_seed(model_seed, input_metabolites, seed_tax_exc_model_name);
        
        % Concatenate modelAGORA and modelSEED
        model_concatenated = ConcatenateModels(model_agora, model_seed, conc_model_name);
        
        % Merge AGORA and SEED models
        model_merged = mergeSEEDmetabolites(model_concatenated, comparison_table, merg_model_name);
        
        % Add information about some metabolites of interest to the network
        model_merged_with_mets = addMetabolites(model_merged, new_metabolites);
        
        % Add the manual curated pathways for i-diet metabolites
        model_merged_expert = addReactions(model_merged_with_mets, new_reactions, id_species, new_metabolites);
        
        % Make uniform the bounds to -1000 and +1000 for those reactions
        % that have lb<0 or ub>0 respectively
        model_merged_expert.lb(model_merged_expert.lb < 0) = -1000;
        model_merged_expert.ub(model_merged_expert.ub > 0) = +1000;
        
        save(fullfile('.', 'output', 'Models', [merg_model_name '.mat']), 'model_merged_expert');
               
        % Apply the FastcoreWeighted algorithm to extract a minimum functional
        % model given a selected core
        core = [find(~cellfun(@isempty,regexp(model_merged_expert.rxnID,'rxnAGORA'))); find(~cellfun(@isempty,regexp(model_merged_expert.rxnID,'rxnAdd')))];
        sense = 'both';

        [finalModel, coreAddedRxns, iDietAddedRxns] = applyFastcoreDiet(model_merged_expert, core, sense, options);
        model = removeNotAnnotatedRxns(finalModel, fastFVA_path);
        model = mergeDupRxns(model);
        
        %Annotate blank exchange reactions if required
        if ANNOTATE
            rxnsWithoutTax = find(cellfun(@isempty, model.trules));
            for i = 1:length(rxnsWithoutTax)
                met = find(model.S(:,rxnsWithoutTax(i)));
                rxns = find(model.S(met,:));
                rxns = rxns(~ismember(rxns,rxnsWithoutTax(i)));
                if length(rxns) > 1
                    tmp = cellfun(@strsplit, model.trules(rxns), repmat({' | '}, length(rxns), 1), ...
                        'UniformOutput', false);
                    tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
                    model.trules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' | ');
                    tmp = cellfun(@strsplit, model.trRules(rxns), repmat({' or '}, length(rxns), 1), ...
                        'UniformOutput', false);
                    tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
                    model.trRules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' or ');
                else
                    model.trules(rxnsWithoutTax(i)) = model.trules(rxns);
                    model.trRules(rxnsWithoutTax(i)) = model.trRules(rxns);
                end
            end
        end
        
        save(fullfile('.', 'output', 'Models', [minim_model_name '.mat']), 'model');
    else
        clc
        disp('Loading modelSEED with taxonomy information')
        load(fullfile('.', 'output', 'Models', [seed_tax_model_name '.mat']));
        
        % Create exchange reactions in the SEED model for those metabolites
        % of idiet that don't one
        model_seed = Create_exchange_reaction_seed(model_seed, input_metabolites, seed_tax_exc_model_name);
        
        % Concatenate modelAGORA and modelSEED
        model_concatenated = ConcatenateModels(model_agora, model_seed, conc_model_name);
        
        % Merge AGORA and SEED models
        model_merged = mergeSEEDmetabolites(model_concatenated, comparison_table, merg_model_name);
        
        % Add information about some metabolites of interest to the network
        model_merged_with_mets = addMetabolites(model_merged, new_metabolites);
        
        % Add the manual curated pathways for i-diet metabolites
        model_merged_expert = addReactions(model_merged_with_mets, new_reactions, id_species, new_metabolites);
        
        % Make uniform the bounds to -1000 and +1000 for those reactions
        % that have lb<0 or ub>0 respectively
        model_merged_expert.lb(model_merged_expert.lb < 0) = -1000;
        model_merged_expert.ub(model_merged_expert.ub > 0) = +1000;
        
        save(fullfile('.', 'output', 'Models', [merg_model_name '.mat']), 'model_merged_expert');
       
        % Apply the FastcoreWeighted algorithm to extract a minimum functional
        % model given a selected core
        core = [find(~cellfun(@isempty,regexp(model_merged_expert.rxnID,'rxnAGORA'))); find(~cellfun(@isempty,regexp(model_merged_expert.rxnID,'rxnAdd')))];
        sense = 'both';

        [finalModel, coreAddedRxns, iDietAddedRxns] = applyFastcoreDiet(model_merged_expert, core, sense, options);
        model = removeNotAnnotatedRxns(finalModel, fastFVA_path);
        model = mergeDupRxns(model);
        
        %Annotate blank exchange reactions if required
        if ANNOTATE
            rxnsWithoutTax = find(cellfun(@isempty, model.trules));
            for i = 1:length(rxnsWithoutTax)
                met = find(model.S(:,rxnsWithoutTax(i)));
                rxns = find(model.S(met,:));
                rxns = rxns(~ismember(rxns,rxnsWithoutTax(i)));
                if length(rxns) > 1
                    tmp = cellfun(@strsplit, model.trules(rxns), repmat({' | '}, length(rxns), 1), ...
                        'UniformOutput', false);
                    tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
                    model.trules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' | ');
                    tmp = cellfun(@strsplit, model.trRules(rxns), repmat({' or '}, length(rxns), 1), ...
                        'UniformOutput', false);
                    tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
                    model.trRules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' or ');
                else
                    model.trules(rxnsWithoutTax(i)) = model.trules(rxns);
                    model.trRules(rxnsWithoutTax(i)) = model.trRules(rxns);
                end
            end
        end
        
        save(fullfile('.', 'output', 'Models', [minim_model_name '.mat']), 'model');
   end
else
    clc
    disp('Loading the AGORA-SEED model merged with sergio reactions')
    load(fullfile('.', 'output', 'Models', [merg_model_name '.mat']));
    
    % Apply the FastcoreWeighted algorithm to extract a minimum functional
    % model given a selected core
    core = [find(~cellfun(@isempty,regexp(model_merged_expert.rxnID,'rxnAGORA'))); find(~cellfun(@isempty,regexp(model_merged_expert.rxnID,'rxnAdd')))];
    sense = 'both';
    
    [finalModel, coreAddedRxns, iDietAddedRxns] = applyFastcoreDiet(model_merged_expert, core, sense, options);
    model = removeNotAnnotatedRxns(finalModel, fastFVA_path);
    model = mergeDupRxns(model);
    
    %Annotate blank exchange reactions if required
    if ANNOTATE
        rxnsWithoutTax = find(cellfun(@isempty, model.trules));
        for i = 1:length(rxnsWithoutTax)
            met = find(model.S(:,rxnsWithoutTax(i)));
            rxns = find(model.S(met,:));
            rxns = rxns(~ismember(rxns,rxnsWithoutTax(i)));
            if length(rxns) > 1
                tmp = cellfun(@strsplit, model.trules(rxns), repmat({' | '}, length(rxns), 1), ...
                    'UniformOutput', false);
                tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
                model.trules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' | ');
                tmp = cellfun(@strsplit, model.trRules(rxns), repmat({' or '}, length(rxns), 1), ...
                    'UniformOutput', false);
                tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
                model.trRules{rxnsWithoutTax(i)} = strjoin(unique(cat(1,tmp{:})),' or ');
            else
                model.trules(rxnsWithoutTax(i)) = model.trules(rxns);
                model.trRules(rxnsWithoutTax(i)) = model.trRules(rxns);
            end
        end
    end
    
    save(fullfile('.', 'output', 'Models', [minim_model_name '.mat']), 'model');
end
toc

%% Section 3: AGREDA balancing
tic
load(fullfile('.', 'output', 'Models', [minim_model_name '.mat']))

%Correct SEED unspurious information
model = removeProblematicMets(model,table2cell(removeMets(:,'ID')));
model = mergeMets(model,table2cell(mergeMetsSEED(:,'manID')),table2cell(mergeMetsSEED(:,'remID')));

%Add exchange reactions to metabolites with enough evidence (from HMDB)
model = addExchangesInfo(model,table2cell(exMetsWithEvidence(:,'ID')));

%Perform the species balancing
model = balanceSpecies(model,FVA_balancing_path);
save(fullfile('.', 'output', 'Models', [final_model_name '.mat']),'model')
toc