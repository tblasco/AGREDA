function [model_seed, species_not_included] = manageEC_pipeline(model_seed, species, peg_path, kegg_path, manual_ECnumbers, manual_trRules, manual_rxntrRules, seed_tax_model_name)

% manageEC_pipeline.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

clc
disp('Adding the taxonomic information in modelSEED');
model_seed.modelName = seed_tax_model_name;

act_dir = pwd;
cd(peg_path)
files = dir('*.txt');
cd(act_dir)
files = {files(:).name}';

species = cellfun(@strrep, species, repmat({'.mat'}, length(species), 1), repmat({'.txt'}, length(species), 1), 'UniformOutput', false);
included_species = species(ismember(species, files));
n_included_species = length(included_species);
species_not_included = species(~ismember(species,files));
species_not_included = cellfun(@strrep, species_not_included, repmat({'.txt'}, length(species_not_included), 1), repmat({''}, length(species_not_included), 1), 'UniformOutput', false);

EC_n = cell(n_included_species,1);
tax = cell(n_included_species,1);
tax_name = cell(n_included_species,1);
for i = 1:n_included_species
    clc
    disp('Adding the taxonomy information in modelSEED');
    disp(['Loading organism EC numbers ' num2str(i) '/' num2str(length(included_species))]);

    tmp = readtable(fullfile(peg_path, included_species{i}));
    tmp = table2cell(tmp(:,end));

    ec_numbers = cellfun(@regexp,tmp,repmat({'([ET]C.?)([\w\-]?\.?)+'},length(tmp),1),repmat({'match'},length(tmp),1),'UniformOutput',false);
    ec_numbers = cellfun(@(x) x', ec_numbers(~cellfun(@isempty,ec_numbers)), 'UniformOutput',false);
    ec_numbers = cat(1,ec_numbers{:});
    ec_numbers = cellfun(@regexprep,ec_numbers,repmat({'EC |EC:|EC|_'},length(ec_numbers),1),repmat({''},length(ec_numbers),1), ...
        'UniformOutput',false);
    ec_numbers = ec_numbers(~cellfun(@isempty,cellfun(@strfind,ec_numbers,repmat({'.'},length(ec_numbers),1),'UniformOutput',false)));
    EC_n{i,1} = unique(ec_numbers);
    tax{i,1} = repmat({strcat('t(',num2str(i),')')},length(EC_n{i,1}),1);
    tax_name{i,1} = repmat(included_species(i),length(EC_n{i,1}),1);
end

whole_SEED_EC = cat(1,EC_n{:});
whole_SEED_tax = cat(1,tax{:});
whole_SEED_tax_name = cat(1,tax_name{:});

cd(kegg_path)
kegg_files = dir('*.txt');
cd(act_dir)
kegg_files = {kegg_files(:).name}';

including_KEGG_species = files(ismember(files,kegg_files));
corresponding_tax = find(ismember(files,kegg_files));
n_KEGG_tax = length(corresponding_tax);
KEGG_EC_n = cell(n_KEGG_tax,1);
KEGG_tax = cell(n_KEGG_tax,1);
KEGG_tax_name = cell(n_KEGG_tax,1);
for i = 1:n_KEGG_tax
    clc
    disp('Adding the taxonomy information in modelSEED');
    disp(['Loading organism EC numbers from KEGG ' num2str(i) '/' num2str(n_KEGG_tax)]);

    tmp = readtable(fullfile(kegg_path, including_KEGG_species{i}));
    tmp = table2cell(tmp(:,end));

    KEGG_EC_n{i,1} = unique(tmp);
    KEGG_tax{i,1} = repmat({strcat('t(',num2str(corresponding_tax(i)),')')},length(KEGG_EC_n{i,1}),1);
    KEGG_tax_name{i,1} = repmat(including_KEGG_species(i),length(KEGG_EC_n{i,1}),1);
end

whole_KEGG_EC = cat(1,KEGG_EC_n{:});
whole_KEGG_tax = cat(1,KEGG_tax{:});
whole_KEGG_tax_name = cat(1,KEGG_tax_name{:});

whole_EC = [whole_SEED_EC;whole_KEGG_EC];
whole_tax = [whole_SEED_tax;whole_KEGG_tax];
whole_tax_name = [whole_SEED_tax_name;whole_KEGG_tax_name];

[unique_whole_EC,~,relations] = unique(whole_EC);
n_unique_EC = length(unique_whole_EC);
EC_x_tax = cell(n_unique_EC,3);
EC_x_tax(:,1) = unique_whole_EC;
for i = 1:n_unique_EC
    clc
    disp('Adding the taxonomy information in modelSEED');
    disp(['Calculating the taxonomy of each EC number ' num2str(i) '/' num2str(n_unique_EC)]);
    index = ismember(relations,i);
    
    EC_x_tax(i,2) = {strjoin(unique(whole_tax(index)),' | ')};
    EC_x_tax(i,3) = {strjoin(unique(whole_tax_name(index)),' or ')};
end
EC_x_tax(:,3) = cellfun(@strrep, EC_x_tax(:,3), repmat({'.txt'}, n_unique_EC, 1), repmat({''}, n_unique_EC, 1), 'UniformOutput', false);

EC_found = table2cell(manual_trRules(:,'ECNumbers'));
trules_found = cell(length(EC_found),1);
trRules_found = table2cell(manual_trRules(:,'trRules'));
name_species = cellfun(@strrep, included_species, repmat({'.txt'}, length(included_species), 1), repmat({''}, length(included_species), 1), ...
    'UniformOutput', false);
n_trRules_found = length(trRules_found);
for i = 1 : n_trRules_found
    tmp = split(trRules_found(i),' or ');
    pos_species = num2cell(find(ismember(name_species,tmp)));
    if ~isempty(pos_species)
        trules_found{i,1} = strjoin(cellfun(@strcat,repmat({'t('},length(pos_species),1),cellfun(@num2str,pos_species,'UniformOutput',false),repmat({')'},length(pos_species),1),'UniformOutput',false), ' | ');
    end
end
EC_x_tax = [EC_x_tax; [EC_found trules_found trRules_found]];

idx = find(cellfun(@isempty,regexp(EC_x_tax(:,1), 'TC')));
base = repmat({'-'},length(EC_x_tax(idx,1)),4);
rxn_EC_splitted = cellfun(@strsplit, EC_x_tax(idx,1), repmat({'.'}, length(EC_x_tax(idx,1)) ,1), 'UniformOutput', false);
for i = 1:length(EC_x_tax(idx,1))
    clc
    disp('Adding the taxonomy information in modelSEED');
    disp(['Creating the reaction-taxonomy relations 1: ' num2str(i) '/' num2str(length(EC_x_tax(idx,1)))]);
    for j = 1:length(rxn_EC_splitted{i})
        base{i,j} = rxn_EC_splitted{i,1}{1,j};
    end
end

for i = 1:size(base,2)
    errors = find(~cellfun(@isempty,cellfun(@regexp,base(:,i),repmat({'\d-'}, size(base,1), 1),'UniformOutput',false)));
    if errors
        base(errors,i) = cellfun(@strrep,base(errors,i),repmat({'-'},length(errors),1),repmat({''},length(errors),1),'UniformOutput',false);
    end
end
cellRows = mat2cell(base,ones(size(base,1),1),size(base,2));
EC_x_tax(idx,1) = cellfun(@strjoin, cellRows, repmat({'.'}, length(cellRows),1), 'UniformOutput', false);

[unique_whole_EC,~,relations] = unique(EC_x_tax(:,1));
n_unique_EC = length(unique_whole_EC);
repeated_EC = unique_whole_EC(hist(relations,unique(relations))>1);
for i = 1:length(repeated_EC)
    idx = find(ismember(EC_x_tax(:,1),repeated_EC(i)));
    tmp = strjoin(EC_x_tax(idx,2), {' | '});
    EC_x_tax(idx(1),2) = {strjoin(unique(strsplit(tmp, ' | ')),' | ')};
    tmp = strjoin(EC_x_tax(idx,3), {' or '});
    EC_x_tax(idx(1),3) = {strjoin(unique(strsplit(tmp, ' or ')),' or ')};
    EC_x_tax(idx(2:end),:) = [];
end

%% Add EC numbers from the bibliography
manual_ECnumbers = [table2cell(manual_ECnumbers(:,'rxnID')) table2cell(manual_ECnumbers(:,'rxnECnumber'))];
manual_ECnumbers = manual_ECnumbers(~cellfun(@isempty, manual_ECnumbers(:,2)),:);
manual_ECnumbers(:,1) = cellfun(@strrep, manual_ECnumbers(:,1), repmat({"'"}, size(manual_ECnumbers,1), 1), repmat({''}, ...
    size(manual_ECnumbers,1), 1), 'UniformOutput', false);
manual_ECnumbers(:,2) = cellfun(@strrep, manual_ECnumbers(:,2), repmat({'EC '}, size(manual_ECnumbers,1), 1), repmat({''}, ...
    size(manual_ECnumbers,1), 1), 'UniformOutput', false);
for i = 1:size(manual_ECnumbers,1)
    if cellfun(@isempty,model_seed.rxnECnumbers(ismember(model_seed.rxnID,manual_ECnumbers{i,1})))
        model_seed.rxnECnumbers(ismember(model_seed.rxnID,manual_ECnumbers{i,1})) = manual_ECnumbers(i,2);
    else
        tmp = join([model_seed.rxnECnumbers(ismember(model_seed.rxnID,manual_ECnumbers{i,1})),manual_ECnumbers(i,2)],'|');
        model_seed.rxnECnumbers(ismember(model_seed.rxnID,manual_ECnumbers{i,1})) = join(unique(split(tmp,'|')),'|');
    end
end

%% Create trules and trRules fields
seed_rxn_ec = find(~cellfun(@isempty,model_seed.rxnECnumbers));
tmp = cellfun(@strsplit, model_seed.rxnECnumbers(seed_rxn_ec), repmat({'|'},length(model_seed.rxnECnumbers(seed_rxn_ec)),1), ...
    'UniformOutput',false);
tmp = cellfun(@(x) x', tmp, 'UniformOutput', false);
rxn_number = [];
rxn_EC = [];
n_rxn = length(seed_rxn_ec);
for i = 1:n_rxn
    rxn_EC = [rxn_EC; unique(cat(1,tmp{i}))];
    rxn_number = [rxn_number; repmat(seed_rxn_ec(i), length(unique(cat(1,tmp{i}))), 1)];
end

base = repmat({'-'},length(rxn_EC),4);
rxn_EC_splitted = cellfun(@strsplit, rxn_EC, repmat({'.'}, length(rxn_EC) ,1), 'UniformOutput', false);

for i = 1:length(rxn_EC)
    clc
    disp('Adding the taxonomy information in modelSEED');
    disp(['Creating the reaction-taxonomy relations 2: ' num2str(i) '/' num2str(length(rxn_EC))]);
    for j = 1:length(rxn_EC_splitted{i})
        base{i,j} = rxn_EC_splitted{i,1}{1,j};
    end
end

for i = 1:size(base,2)
    errors = find(~cellfun(@isempty,cellfun(@regexp,base(:,i),repmat({'\d-'}, size(base,1), 1),'UniformOutput',false)));
    if errors
        base(errors,i) = cellfun(@strrep,base(errors,i),repmat({'-'},length(errors),1),repmat({''},length(errors),1),'UniformOutput',false);
    end
end
cellRows = mat2cell(base,ones(size(base,1),1),size(base,2));
rxn_EC = cellfun(@strjoin, cellRows, repmat({'.'}, length(cellRows),1), 'UniformOutput', false);

model_seed.trules = cell(length(model_seed.rxnAbbreviation),1);
model_seed.trRules = cell(length(model_seed.rxnAbbreviation),1);
for i = 1:n_unique_EC
    pos = rxn_number(ismember(rxn_EC, EC_x_tax{i,1}));
    pos_empty = pos(cellfun(@isempty,model_seed.trules(pos)));
    pos_not_empty = pos(~cellfun(@isempty,model_seed.trules(pos)));
    
    if ~isempty(pos_empty)
        model_seed.trules(pos_empty,1) = EC_x_tax(i,2);
        model_seed.trRules(pos_empty,1) = EC_x_tax(i,3);
    end
    
    if ~isempty(pos_not_empty)
        tmp = strcat(model_seed.trules(pos_not_empty), {' | '}, EC_x_tax(i,2));
        tmp = cellfun(@unique, cellfun(@strsplit, tmp, repmat({' | '}, length(tmp), 1), 'UniformOutput', false), 'UniformOutput', false);
        model_seed.trules(pos_not_empty) = cellfun(@strjoin, tmp, repmat({' | '}, length(tmp) , 1), 'UniformOutput', false);
        tmp = strcat(model_seed.trRules(pos_not_empty), {' or '}, EC_x_tax(i,3));
        tmp = cellfun(@unique, cellfun(@strsplit, tmp, repmat({' or '}, length(tmp), 1), 'UniformOutput', false), 'UniformOutput', false);
        model_seed.trRules(pos_not_empty) = cellfun(@strjoin, tmp, repmat({' or '}, length(tmp) , 1), 'UniformOutput', false);
    end
    
end

%% Add taxonomic information from the bibliography
manual_rxntrRules = [table2cell(manual_rxntrRules(:,'rxnID')) table2cell(manual_rxntrRules(:,'trRules'))];
for i = 1:size(manual_rxntrRules,1)
    tmp = split(manual_rxntrRules{i,2},' or ');
    pos_species = num2cell(find(ismember(name_species,tmp)));
    if ~isempty(pos_species)
        if cellfun(@isempty,model_seed.trules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))))
            model_seed.trules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))) = {strjoin(cellfun(@strcat,repmat({'t('},length(pos_species),1),cellfun(@num2str,pos_species,'UniformOutput',false),repmat({')'},length(pos_species),1),'UniformOutput',false),' | ')};
            model_seed.trRules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))) = {manual_rxntrRules{i,2}};
        else
            tmp2 = join([model_seed.trules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))),{strjoin(cellfun(@strcat,repmat({'t('},length(pos_species),1),cellfun(@num2str,pos_species,'UniformOutput',false),repmat({')'},length(pos_species),1),'UniformOutput',false),' | ')}],' | ');
            model_seed.trules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))) = join(unique(split(tmp2,' | ')),' | ');
            tmp2 = join([model_seed.trRules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))),{manual_rxntrRules{i,2}}],' | ');
            model_seed.trRules(~cellfun(@isempty,regexp(model_seed.rxnID,manual_rxntrRules{i,1}))) = join(unique(split(tmp2,' | ')),' | ');
        end
    end
end

%% Create the rxnTaxMat field
model_seed.rxnTaxMat = sparse(length(model_seed.rxnID), n_included_species);
index_rxns_tax = find(~cellfun(@isempty,model_seed.trules));
pos = cellfun(@regexp, model_seed.trules(index_rxns_tax), repmat({'\d*'}, length(model_seed.trules(index_rxns_tax)), 1), repmat({'Match'}, length(model_seed.trules(index_rxns_tax)), 1), 'UniformOutput', false);
pos = cellfun(@str2double, pos, 'UniformOutput', false);
n_rxns = length(index_rxns_tax);
for i = 1:n_rxns
    clc
    disp('Adding the taxonomy information in modelSEED');
    disp(['Creating the reaction-taxonomy relations 3: ' num2str(i) '/' num2str(n_rxns)]);
    model_seed.rxnTaxMat(index_rxns_tax(i), pos{i}) = 1;
end
end