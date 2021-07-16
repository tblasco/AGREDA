function [model,count_present, count_not_present] = Create_exchange_reaction_seed(model, input_metabolites, so_model_name)

% Create_exchange_reaction.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

clc
disp('Creating the exchange reaction for the mets in our diet');
model.modelName = so_model_name;

present_agora = input_metabolites(:,'AGORA');
present_agora = table2cell(present_agora);

present_seed = input_metabolites(:,'SEED');
present_seed = table2cell(present_seed);

mets_in_seed = input_metabolites(:,'CPD');
mets_in_seed = table2cell(mets_in_seed);

mets_to_remove_seed = double(cell2mat(present_seed)==0) + double(cell2mat(present_agora)~=0);

mets_in_seed(logical(mets_to_remove_seed)) = [];
mets_in_seed = mets_in_seed(~cellfun(@isempty,mets_in_seed));

count_present = 0;
count_not_present = 0;
cont = 1;
for i = 1:length(mets_in_seed)
    pos=find(strcmp(model.metID,mets_in_seed{i}));
    if ~isempty(pos)
        model.S = [model.S, zeros(length(model.metID),1)];
        model.S(pos,end) = -1;
        model.ub = [model.ub;1000];
        model.lb = [model.lb;-1000];
        model.rxnAbbreviation = [model.rxnAbbreviation;strcat('EX_iDiet_',model.metNames{pos})];
        model.rxnNames = [model.rxnNames;strcat('EX_iDiet_',model.metNames{pos})];
        model.rxnECnumbers = [model.rxnECnumbers; cell(1,1)];
        model.correctedRxns = [model.correctedRxns; 0];
        model.rxnID = [model.rxnID;strcat('rxnDiet',num2str(cont))];
        model.rxnCode = [model.rxnCode; cell(1,1)];
        model.rxnStoichiometry = [model.rxnStoichiometry; cell(1,1)];
        model.rxnIsTransport = [model.rxnIsTransport; cell(1,1)];
        model.rxnEquation = [model.rxnEquation; cell(1,1)];
        model.rxnDefinition = [model.rxnDefinition; cell(1,1)];
        model.rxnReversibility = [model.rxnReversibility; cell(1,1)];
        model.rxnDirection = [model.rxnDirection; cell(1,1)];
        model.rxnAbstractReaction = [model.rxnAbstractReaction; cell(1,1)];
        model.rxnPathways = [model.rxnPathways; cell(1,1)];
        model.rxnAliases = [model.rxnAliases; cell(1,1)];
        model.rxnDeltaG = [model.rxnDeltaG; cell(1,1)];
        model.rxnDeltaGErr = [model.rxnDeltaGErr; cell(1,1)];
        model.rxnCompoundIDs = [model.rxnCompoundIDs; cell(1,1)];
        model.rxnStatus = [model.rxnStatus; cell(1,1)];
        model.rxnIsObsolete = [model.rxnIsObsolete; cell(1,1)];
        model.rxnLinkedReaction = [model.rxnLinkedReaction; cell(1,1)];
        model.rxnNotes = [model.rxnNotes; cell(1,1)];
        model.rxnSource = [model.rxnSource; cell(1,1)];
        model.rxnOntology = [model.rxnOntology; cell(1,1)];
        model.trules = [model.trules; cell(1,1)];
        model.trRules = [model.trRules; cell(1,1)];
        model.rxnTaxMat = [model.rxnTaxMat; sparse(1,size(model.rxnTaxMat,2))];
        
        count_present = count_present+1;
        cont = cont+1;
    else
        count_not_present = count_not_present+1;
    end
end

end