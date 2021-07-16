function model = addReactions(model,new_reactions_table,id_species,new_metabolites)

% addReactions.m
% 
% Author: Francesco Balzerani, Telmo Blasco, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

rxnEquation_file = table2cell(new_reactions_table(:,1));
Stoichiometry_file = table2cell(new_reactions_table(:,2));
lb_file = table2cell(new_reactions_table(:,3));
ub_file = table2cell(new_reactions_table(:,4));
trRules_file = table2cell(new_reactions_table(:,5));

n_rows = length(rxnEquation_file);

% Process species regarding information

ncbi_names = table2cell(id_species(:,2));
kegg_names = table2cell(id_species(:,4));
agora_names_strains = table2cell(id_species(:,5));
agora_names = table2cell(id_species(:,6));
tax_names = model.taxonomy;

% Process new reaction´s taxonomic information
trules_file = cell(n_rows,1);

for i = 1 : n_rows
    tmp = split(trRules_file(i),', ');
    for j = 1 : length(tmp)
        tmp_1 = split(tmp{j},' ');
        if length(tmp_1) == 1
            ncbi_mask = ~cellfun(@isempty,regexp(ncbi_names,['^',tmp_1{1}]));
            kegg_mask = ~cellfun(@isempty,regexp(kegg_names,['^',tmp_1{1}]));
            agora_strain_mask = ~cellfun(@isempty,regexp(agora_names_strains,['^',tmp_1{1}]));
            agora_names_mask = ~cellfun(@isempty,regexp(agora_names,['^',tmp_1{1}]));
            tax_names_mask = ~cellfun(@isempty,regexp(tax_names,['^',tmp_1{1}]));
            final_mask = ncbi_mask | kegg_mask | agora_strain_mask | agora_names_mask | tax_names_mask;
            pos_tax{j,1} = find(final_mask);
        else
            if isequal(tmp_1{2},'spp') || isequal(tmp_1{2},'strains')
                ncbi_mask = ~cellfun(@isempty,regexp(ncbi_names,['^',tmp_1{1}]));
                kegg_mask = ~cellfun(@isempty,regexp(kegg_names,['^',tmp_1{1}]));
                agora_strain_mask = ~cellfun(@isempty,regexp(agora_names_strains,['^',tmp_1{1}]));
                agora_names_mask = ~cellfun(@isempty,regexp(agora_names,['^',tmp_1{1}]));
                tax_names_mask = ~cellfun(@isempty,regexp(tax_names,['^',tmp_1{1}]));
                final_mask = ncbi_mask | kegg_mask | agora_strain_mask | agora_names_mask | tax_names_mask;
                pos_tax{j,1} = find(final_mask);
            else
                compare_name = join(tmp_1(1:2));
                if strcmp(compare_name,'Bifidobacterium lactis')
                    compare_name = 'Bifidobacterium animalis subsp. lactis';
                end
                if strcmp(compare_name,'Peptostreptococcus productus')
                    compare_name = 'Blautia producta';
                end
                ncbi_mask = ~cellfun(@isempty,regexp(ncbi_names,compare_name));
                kegg_mask = ~cellfun(@isempty,regexp(kegg_names,compare_name));
                agora_strain_mask = ~cellfun(@isempty,regexp(agora_names_strains,compare_name));
                agora_names_mask = ~cellfun(@isempty,regexp(agora_names,compare_name));
                tax_names_mask = ~cellfun(@isempty,regexp(tax_names,['^',tmp_1{1}]));
                final_mask = ncbi_mask | kegg_mask | agora_strain_mask | agora_names_mask | tax_names_mask;
                pos_tax{j,1} = find(final_mask);
            end
        end
    end
    if ~isempty(pos_tax)
        whole_tax = num2cell(cat(1,pos_tax{:}));
        trules_file{i,1} = strjoin(cellfun(@strcat,repmat({'t('},length(whole_tax),1),cellfun(@num2str,whole_tax,'UniformOutput',false),repmat({')'},length(whole_tax),1),'UniformOutput',false), ' | ');
        trRules_file{i,1} = strjoin(model.taxonomy(cell2mat(whole_tax)),' or ');
    end
    clear pos_tax
end
            
% Calculate rxnTaxMat field

rxnTaxMat_file = sparse(length(trRules_file), length(model.taxonomy));
pos = cellfun(@regexp, trules_file, repmat({'\d*'}, n_rows, 1), repmat({'Match'}, n_rows, 1), 'UniformOutput', false);
pos = cellfun(@str2double, pos, 'UniformOutput', false);

for i = 1:n_rows
    rxnTaxMat_file(i, pos{i}) = 1;
end

% Add expert knowledge reactions to the model
tmp_split = cellfun(@strsplit, rxnEquation_file, repmat({'=>'}, n_rows, 1), 'UniformOutput', false);
cpd_id_mets = cellfun(@regexp, tmp_split, repmat({'cpd(\d)+'}, n_rows, 1), repmat({'Match'}, n_rows, 1), 'UniformOutput', false);
stoichiometry_values = cellfun(@strsplit, Stoichiometry_file, repmat({';'}, n_rows,1), 'UniformOutput', false);
S_file = zeros(length(model.mets),n_rows);
for i = 1 : n_rows
    clc
    disp('Creating the new reactions');
    disp(['Reaction: ' num2str(i) '/' num2str(n_rows)]);
    
    for j = 1 : length(cpd_id_mets{i}{1})
        pos_mets = find(~cellfun(@isempty,regexp(model.metID,cpd_id_mets{i}{1}{j})));
        if length(pos_mets) > 1
            pos_cyt = ~cellfun(@isempty,regexp(model.mets(pos_mets),'\[c\]'));
            pos_mets = pos_mets(pos_cyt);
        end
        S_file(pos_mets,i) = -str2double(stoichiometry_values{i}{j});
    end

    count_prod = length(cpd_id_mets{i}{1});
    for k = 1 : length(cpd_id_mets{i}{2})
        pos_mets = find(~cellfun(@isempty,regexp(model.metID,cpd_id_mets{i}{2}{k})));
        if length(pos_mets) > 1
            pos_cyt = ~cellfun(@isempty,regexp(model.mets(pos_mets),'\[c\]'));
            pos_mets = pos_mets(pos_cyt);
        end
        count_prod = count_prod + 1;
        S_file(pos_mets,i) = str2double(stoichiometry_values{i}{count_prod});
    end

    cpd_id_mets{i}(1,1) = cellfun(@(x) x', cpd_id_mets{i}(1,1), 'UniformOutput',false);
    cpd_id_mets{i}(1,2) = cellfun(@(x) x', cpd_id_mets{i}(1,2), 'UniformOutput',false);
    cpd_id_mets{i} = cellfun(@(x) x', cpd_id_mets(i), 'UniformOutput',false);
    cpd_id_mets{i} = cat(1,cpd_id_mets{i}{1});
    cpd_id_mets{i} = cat(1,cpd_id_mets{i}{:});
end

% Complete rxnNames and rxnID fields
rxnNames_file = cell(length(trules_file),1);
rxnID_file = cell(length(trules_file),1);
for i = 1 : n_rows
    rxnNames_file{i} = ['new_reaction_',num2str(i)];
    rxnID_file{i} = ['rxnAdd',num2str(i)];
end

model.S = [model.S,S_file];
model.rxns = [model.rxns;rxnNames_file];
model.lb = [model.lb;cell2mat(lb_file)];
model.ub = [model.ub;cell2mat(ub_file)];
model.c  = [model.c;zeros(length(rxnNames_file),1)];
model.rules = [model.rules;cell(length(rxnNames_file),1)];
model.rxnGeneMat = [model.rxnGeneMat ; sparse(length(rxnNames_file), length(model.genes))];
model.rxnNames = [model.rxnNames;rxnNames_file];
model.subSystems = [model.subSystems;cell(length(rxnNames_file),1)];
model.grRules = [model.grRules;cell(length(rxnNames_file),1)];
model.comments = [model.comments;cell(length(rxnNames_file),1)];
model.citations = [model.citations;cell(length(rxnNames_file),1)];
model.rxnConfidenceScores = [model.rxnConfidenceScores;zeros(length(rxnNames_file),1)];
model.rxnECNumbers = [model.rxnECNumbers;cell(length(rxnNames_file),1)];
model.rxnKEGGID = [model.rxnKEGGID;cell(length(rxnNames_file),1)];
model.trRules = [model.trRules; trRules_file];
model.trules = [model.trules; trules_file];
model.rxnTaxMat = [model.rxnTaxMat; rxnTaxMat_file ];
model.correctedRxns = [model.correctedRxns;zeros(length(rxnNames_file),1)];
model.rxnID = [model.rxnID;rxnID_file];
model.rxnCode = [model.rxnCode;cell(length(rxnNames_file),1)];
model.rxnStoichiometry = [model.rxnStoichiometry;Stoichiometry_file];
model.rxnIsTransport = [model.rxnIsTransport;cell(length(rxnNames_file),1)];
model.rxnEquation = [model.rxnEquation;rxnEquation_file];
model.rxnDefinition = [model.rxnDefinition;cell(length(rxnNames_file),1)];
model.rxnReversibility = [model.rxnReversibility;cell(length(rxnNames_file),1)];
model.rxnDirection = [model.rxnDirection;cell(length(rxnNames_file),1)];
model.rxnAbstractReaction = [model.rxnAbstractReaction;cell(length(rxnNames_file),1)];
model.rxnPathways = [model.rxnPathways;cell(length(rxnNames_file),1)];
model.rxnAliases = [model.rxnAliases;cell(length(rxnNames_file),1);];
model.rxnDeltaG = [model.rxnDeltaG;cell(length(rxnNames_file),1)];
model.rxnDeltaGErr = [model.rxnDeltaGErr;cell(length(rxnNames_file),1)];
model.rxnCompoundIDs = [model.rxnCompoundIDs;cell(length(rxnNames_file),1)];
model.rxnStatus = [model.rxnStatus;cell(length(rxnNames_file),1)];
model.rxnIsObsolete = [model.rxnIsObsolete;cell(length(rxnNames_file),1)];
model.rxnLinkedReaction = [model.rxnLinkedReaction;cell(length(rxnNames_file),1)];
model.rxnNotes = [model.rxnNotes;cell(length(rxnNames_file),1)];
model.rxnSource = [model.rxnSource;cell(length(rxnNames_file),1)];
model.rxnOntology = [model.rxnOntology;cell(length(rxnNames_file),1)];

%% Add exchange outside flux for those metabolites of SEED without it
cpd_id_mets = unique(cat(1,cpd_id_mets{:}));

new_cpd_id_mets = table2cell(new_metabolites(:,'metID'));
cpd_oi = cpd_id_mets(~ismember(cpd_id_mets,new_cpd_id_mets));

mets_to_add_exchange = [];
for i = 1:length(cpd_oi)
    %CPD ID with more than one metabolite are supposed to be from AGORA and
    %therefore, they already contain an exchange reaction
    if sum(~cellfun(@isempty,regexp(model.metID,cpd_oi(i)))) > 1
        continue
    else
        idx = find(~cellfun(@isempty,regexp(model.metID,cpd_oi(i))));
        %Checking the exchange reaction of metabolites from SEED
        tmp_S = model.S(:,find(model.S(idx,:)));
        if ~isempty(find(sum(tmp_S~=0, 1)==1))
            continue
        else
            mets_to_add_exchange = [mets_to_add_exchange idx];
        end 
    end
end

%Create the S matrix
tmp_S = zeros(length(model.metID),length(mets_to_add_exchange));
for i = 1:length(mets_to_add_exchange)
    tmp_S(mets_to_add_exchange(i),i) = -1;
end

n_mets = length(mets_to_add_exchange);
id = cellfun(@strcat, repmat({'rxnAddNewEX'},n_mets,1), cellfun(@(x) num2str(x), ...
    num2cell((1:n_mets)'),'UniformOutput',false),'UniformOutput',false);

%Restore the rest of the fields
model.rxnID = [model.rxnID;id];
model.rxns = [model.rxns;strcat('EX_',model.metNames(mets_to_add_exchange))];
model.rxnNames = [model.rxnNames;strcat('EX_',model.metNames(mets_to_add_exchange))];
model.S = [model.S, tmp_S];
model.lb = [model.lb;zeros(n_mets,1)];
idx1 = mets_to_add_exchange(ismember(model.metNames(mets_to_add_exchange),'tetrahydrocurcumin'));
tmp = find(model.S(idx1,:));
model.lb(tmp(find(sum(model.S(:,tmp)~=0,1)==1))) = -1000;
idx1 = mets_to_add_exchange(ismember(model.metNames(mets_to_add_exchange),'Ox-FAD-Flavoproteins'));
tmp = find(model.S(idx1,:));
model.lb(tmp(find(sum(model.S(:,tmp)~=0,1)==1))) = -1000;
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

