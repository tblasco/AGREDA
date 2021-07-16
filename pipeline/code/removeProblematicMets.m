function model = removeProblematicMets(model,metID)

% removeProblematicMets.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

% Extract reactions to delete
rxns_to_delete = [];
for i = 1:length(metID)
    rxns_to_delete = [rxns_to_delete find(model.S(ismember(model.metID,metID(i)),:))];
end
rxns_to_delete = unique(rxns_to_delete);

model = removeRxnsAndMets(model,rxns_to_delete);

end