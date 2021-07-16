function model = balanceSpecies(model,path)

% balanceSpecies.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

% Perform a FVA for each independent model
nT = length(model.taxonomy);
for i = 1:nT
    fprintf('Analyzing species %d of %d\n',i,nT)
    if ~exist(fullfile(path,[model.taxonomy{i} '.mat']),'file')
        m = mapToAGREDA(model,model.taxonomy(i));
        [a,b] = fastFVA(m,0);
        save(fullfile(path,[model.taxonomy{i} '.mat']),'a','b')
    end
end

% Load flux variability of each independent model
output = zeros(length(model.rxnID),length(model.taxonomy));
nT = length(model.taxonomy);
for i = 1:nT
    fprintf('Analyzing species %d of %d\n',i,nT)
    m = mapToAGREDA(model,model.taxonomy(i));
    load(fullfile(path,[model.taxonomy{i} '.mat']));
    index = find(abs(a)<1e-08 & abs(b)<1e-08);
    id = m.rxnID(index);
    output(ismember(model.rxnID,id),i) = 1;
end

% Remove taxes from the model
model_nuevo = model;
for i = 1:size(output,2)
    fprintf('Analyzing species %d of %d\n',i,nT)
    if sum(output(:,i))==0
        continue
    end
    rxn = find(output(:,i));
    model_nuevo.trRules(rxn) = cellfun(@strrep, model_nuevo.trRules(rxn), repmat(model.taxonomy(i), length(rxn), 1), repmat({'ELIMINAR'}, length(rxn), 1), 'UniformOutput', false); 
end

% Correct trRules and detect rxns to delete
index = find(~cellfun(@isempty, regexp(model_nuevo.trRules,'ELIMINAR')));
model_nuevo2 = model_nuevo;
rxns_to_delete = [];
for i = 1:length(index)
    
    tmp_trRules = split(model_nuevo2.trRules(index(i)), ' or ');
    pos = find(~cellfun(@isempty,regexp(tmp_trRules,'ELIMINAR')));
    
    if length(pos) == length(tmp_trRules)
        rxns_to_delete = [rxns_to_delete, index(i)];
        continue
    end
    
    notpos = find(cellfun(@isempty,regexp(tmp_trRules,'ELIMINAR')));
    model_nuevo2.trRules(index(i)) = join(tmp_trRules(notpos), ' or ');
end

model_nuevo3 = removeRxnsAndMets(model_nuevo2,rxns_to_delete);
model = model_nuevo3;

% Correct trules field
for i = 1:length(model.rxns)
    tmp = split(model.trRules(i),' or ');
    trules = cellfun(@num2str, num2cell(find(ismember(model.taxonomy,tmp))), 'UniformOutput', false);
    model.trules(i) = join(join([repmat({'t('}, length(trules), 1),trules,repmat({')'}, length(trules), 1)],''), ' | ');
end

end