function model = Create_exchange_reaction_agora(model, input_metabolites, so_model_name)

% Create_exchange_reaction_agora.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

clc
disp('Creating the exchange reaction in agora for the mets in our diet');
model.modelName = so_model_name;

present_agora = input_metabolites(:,'AGORA');
present_agora = table2cell(present_agora);

mets_in_agora = input_metabolites(:,'AGORA_SYNONYM');
mets_in_agora = table2cell(mets_in_agora);

mets_in_agora(find(cell2mat(present_agora)==0)) = [];
mets_in_agora = mets_in_agora(~cellfun(@isempty,mets_in_agora));

% Create the exchange reaction information
for j = 1:length(mets_in_agora)
    
    % Find metabolite index
    pos=find(strcmp(model.metNames,mets_in_agora{j}));
    pos_ext = find(~cellfun(@isempty, regexp(model.mets(pos),'\[e\]')));
    if isempty(pos_ext)
        model.S = [model.S, zeros(length(model.mets),1)
            zeros(1,length(model.rxns)), -1];
        
        model.mets = [model.mets;strcat(strtok(model.mets(pos),'['),'[e]')];
        model.metNames = [model.metNames;model.metNames(pos)];
        
        model.S = [model.S, zeros(length(model.mets),1)];
        
        model.S(pos,end) = 1;
        model.S(end,end) = -1;
        
        tmp_taxonomy = split(model.metsTax{pos},' or ');
        tmp_rxnTax = ismember(model.taxonomy,tmp_taxonomy)';
        num_tax = length(tmp_taxonomy);
        
        model.ub = [model.ub;num_tax*1000;num_tax*1000];
        model.lb = [model.lb;num_tax*-1000;num_tax*-1000];
        model.rxns = [model.rxns;strcat('EX_',strtok(model.mets(pos),'['),'(e)');strcat(model.mets(end),'_',model.mets(pos))];
        model.rxnNames = [model.rxnNames;strcat('EX_',strtok(model.mets(pos),'['),'(e)');strcat(model.mets(end),'_',model.mets(pos))];
        
        model.c = [model.c; 0; 0];
        model.b = [model.b; 0];
        model.rules = [model.rules; model.rules(pos);model.rules(pos)];
        model.genes = model.genes;
        model.csense = [model.csense;'E'];
        model.rxnGeneMat = [model.rxnGeneMat; zeros(1,length(model.genes));zeros(1,length(model.genes))];
        model.subSystems = [model.subSystems; cell(2,1)];
        model.grRules = [model.grRules; model.grRules(pos);model.grRules(pos)];
        model.comments = [model.comments; cell(2,1)];
        model.citations = [model.citations; cell(2,1)];
        model.rxnConfidenceScores = [model.rxnConfidenceScores; 0; 0];
        model.rxnECNumbers = [model.rxnECNumbers; cell(2,1)];
        model.rxnKEGGID = [model.rxnKEGGID; cell(2,1)];
        model.metCharges = [model.metCharges;model.metCharges(pos)];
        model.metFormulas = [model.metFormulas;model.metFormulas(pos)];
        model.metChEBIID = [model.metChEBIID;model.metChEBIID(pos)];
        model.metKEGGID = [model.metKEGGID;model.metKEGGID(pos)];
        model.metPubChemID = [model.metPubChemID;model.metPubChemID(pos)];
        model.metInChIString = [model.metInChIString;model.metInChIString(pos)];
        model.metHMDBID = [model.metHMDBID;model.metHMDBID(pos)];
        model.metSmiles = [model.metSmiles;model.metSmiles(pos)];
        model.modelID = model.modelID;
        model.taxonomy = model.taxonomy;
        model.metsTax = [model.metsTax;model.metsTax(pos)];
        model.trRules = [model.trRules;model.metsTax(pos);model.metsTax(pos)];
        
        index_trules = find(tmp_rxnTax~=0);
        tmp = cell(num_tax,1);
        for k =1:num_tax
            tmp{k,1} = strcat('t(',num2str(index_trules(k)),')');
        end
        
        tmp_trules = strjoin(tmp,' | '); 
        model.trules = [model.trules;tmp_trules;tmp_trules];
        model.rxnTaxMat = [model.rxnTaxMat;tmp_rxnTax;tmp_rxnTax];

    end
end

%Add the rxnID field to the AGORA model
model.rxnID = cellfun(@strcat, repmat({'rxnAGORA'},length(model.rxns),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(model.rxns))'),'UniformOutput',false),'UniformOutput',false);

%Add the metID field to the AGORA model
model.metID = cellfun(@strcat, repmat({'cpdAGORA'},length(model.mets),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(model.mets))'),'UniformOutput',false),'UniformOutput',false);

end