function model = mapToAGREDA(model, sp)

% mapToAGREDA.R
%
% Function to contextualize a GSM model based on the
% species present in a sample.
%
% Variables:
%
%   * Inputs:
%       - model: GSM model (cobratoolbox format).
%       - sp: cell containing ASV information of the
%             species present in a sample. Each ASV
%             can be related to multiple species 
%             separated by ";".
%
%   * Outputs:
%       - model: GSM model (cobratoolbox format).
%
% Usage:
%
%       model_mapped <- mapToAGREDA(model, sp)
%
% Author: Telmo Blasco
% Email: tblasco@tecnun.es

% Check if the species are introduced in a cell
if ~iscell(sp)
    error('Introduce the species in a cell array')
end

% Extract number of species
n = length(sp);

% Reactions without taxonomy are not considered
wRxns = cellfun(@isempty, model.trRules);

% Check if ASV have many species related
tmp = cellfun(@strsplit, sp, repmat({';'}, n, 1), 'UniformOutput', false);

% Extract the index of the reactions related to the species
index = false(length(model.rxns),1);
for i = 1:n
    % For ASV with multiple species, extract common reactions
    idx = true(length(model.rxns),1);
    for j = 1:length(tmp{i,1})
        idx = idx & ~cellfun(@isempty, regexp(model.trRules,tmp{i,1}{1,j}));
    end
    index = index | idx;
end

% Remove the unnecessary metabolites and reactions
model = removeRxnsAndMets(model, find(~(index | wRxns)));

end