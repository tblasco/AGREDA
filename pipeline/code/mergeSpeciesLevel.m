function [model, not_included_species] = mergeSpeciesLevel(species_path, species, sl_model_name)

% mergeSpeciesLevel.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

act_dir = pwd;
cd(species_path)
files = dir('*.mat');
cd(act_dir)
files = {files(:).name}';
not_included_species = species(~ismember(species, files));
included_species = species(ismember(species, files));

%% Load the models
n_included_species = length(included_species);
for i = 1:n_included_species
    clc
    disp(['Loading Models ' num2str(i) '/' num2str(n_included_species)]);
    tmp = load(fullfile(species_path, species{i}));
    all_models{i, 1} = tmp.model;
end
clear tmp

%% Merge the models
included_species_name = cellfun(@(x) x(1:end-4), included_species, 'UniformOutput', false);
for i = 1:n_included_species
    clc
    disp(['Merging Species ' num2str(i) '/' num2str(n_included_species)]);
    % Prepare the grRules field
    n_act_rxns = length(all_models{i}.rxns);
    for j = 1:n_act_rxns
        words = strsplit(all_models{i}.grRules{j},' ');
        n_act_words = length(words);
        for k = 1:n_act_words
            if length(words{k}) > 3
                words{1, k} = strcat(included_species_name{i}, '_', words{k});
            end
        end
        all_models{i}.grRules{j} = strjoin(words);
    end
    % Prepare the genes field
    all_models{i}.genes = strcat(included_species_name{i}, '_', all_models{i}.genes);
    
    % Merge the models
    if i == 1
        model = all_models{i};
        model.modelID = {model.modelID};
        model.modelName = sl_model_name;
        model.taxonomy = included_species_name(i);
        model.metsTax = repmat(included_species_name(i), length(model.mets), 1);
        model.trRules = repmat(included_species_name(i), length(model.rxns), 1);
        model.trules = repmat({['t(', num2str(i) ')']}, length(model.rxns), 1);
    else
        % Prepare the rules field
        act_model = all_models{i};
        n_genes_model = length(model.genes);
        n_genes_act_model = length(act_model.genes);
        
        gene_id_relation(:, 1) = cellstr(string(1:n_genes_act_model));
        gene_id_relation(:, 2) = cellstr(string((n_genes_model+1):(n_genes_model+n_genes_act_model)));     
        gene_id_relation = cellfun(@strcat, repmat({'x('}, n_genes_act_model, 2), gene_id_relation, 'UniformOutput', false);
        gene_id_relation = cellfun(@strcat, gene_id_relation, repmat({')'}, n_genes_act_model, 2), 'UniformOutput', false);
        n_act_genes = size(gene_id_relation, 1);
        
        for j = 1:n_act_genes
            act_model.rules = strrep(act_model.rules, gene_id_relation(j, 1), gene_id_relation(j, 2));
        end
        
        % Merge the models
        model.rxns = [model.rxns; act_model.rxns];
        model.S = [model.S, sparse(size(model.S, 1), size(act_model.S, 2)); sparse(size(act_model.S, 1), size(model.S, 2)), act_model.S];
        model.lb = [model.lb; act_model.lb];
        model.ub = [model.ub; act_model.ub];
        model.c = [model.c; act_model.c];
        model.mets = [model.mets; act_model.mets];
        model.b = [model.b; act_model.b];
        model.rules = [model.rules; act_model.rules];
        model.genes = [model.genes; act_model.genes];
        model.csense = [model.csense; act_model.csense];
        model.rxnGeneMat = [model.rxnGeneMat, sparse(size(model.rxnGeneMat, 1), size(act_model.rxnGeneMat, 2)); sparse(size(act_model.rxnGeneMat, 1), size(model.rxnGeneMat, 2)), act_model.rxnGeneMat];
        model.metNames = [model.metNames; act_model.metNames];
        model.subSystems = [model.subSystems; act_model.subSystems];
        model.grRules = [model.grRules; act_model.grRules];
        model.rxnNames = [model.rxnNames; act_model.rxnNames];
        model.comments = [model.comments; act_model.comments];
        model.citations = [model.citations; act_model.citations];
        model.rxnConfidenceScores = [model.rxnConfidenceScores; act_model.rxnConfidenceScores];
        model.rxnECNumbers = [model.rxnECNumbers; act_model.rxnECNumbers];
        model.rxnKEGGID = [model.rxnKEGGID; act_model.rxnKEGGID];
        model.metCharges = [model.metCharges; act_model.metCharges];
        model.metFormulas = [model.metFormulas; act_model.metFormulas];
        model.metChEBIID = [model.metChEBIID; act_model.metChEBIID];
        model.metKEGGID = [model.metKEGGID; act_model.metKEGGID];
        model.metPubChemID = [model.metPubChemID; act_model.metPubChemID];
        model.metInChIString = [model.metInChIString; act_model.metInChIString];
        model.metHMDBID = [model.metHMDBID; act_model.metHMDBID];
        model.metSmiles = [model.metSmiles; act_model.metSmiles];
        model.modelID = [model.modelID; {act_model.modelID}];
        model.taxonomy = [model.taxonomy; included_species_name(i)];
        model.metsTax = [model.metsTax; repmat(included_species_name(i), length(act_model.mets), 1)];
        model.trRules = [model.trRules; repmat(included_species_name(i), length(act_model.rxns), 1)];
        model.trules = [model.trules; repmat({['t(', num2str(i) ')']}, length(act_model.rxns), 1)];
        clear gene_id_relation
    end
end
end