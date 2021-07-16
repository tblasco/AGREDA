function [finalModel, coreAddedRxns, iDietAddedRxns] = applyFastcoreDiet(model, core, sense, options)

% applyFastcoreDiet.m
%
% Authors: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es, 
% fplanes@tecnun.es

% Set parameters
if nargin == 3
    options.modelName = 'modelFastcoreDiet';
    options.plantRxns = [];
    options.fastFVAName = 'fastFVAFastcoreDiet';
    options.blockedDietRxnsUp = [];
    options.blockedDietRxnsOut = [];
    options.rxnDietID = 'rxnDiet';
else
    if isfield(options, 'modelName')
        modelName = options.modelName;
    else
        modelName = 'modelFastcoreDiet';
    end
    if isfield(options, 'plantRxns')
        plantRxns = options.plantRxns;
    else
        plantRxns = [];
    end
    if isfield(options, 'fastFVAName')
        fastFVAName = options.fastFVAName;
    else
        fastFVAName = 'fastFVAFastcoreDiet';
    end
    if isfield(options, 'blockedDietRxnsUp')
        blockedDietRxnsUp = options.blockedDietRxnsUp;
    else
        blockedDietRxnsUp = [];
    end
    if isfield(options, 'blockedDietRxnsOut')
        blockedDietRxnsOut = options.blockedDietRxnsOut;
    else
        blockedDietRxnsOut = [];
    end
    if isfield(options, 'rxnDietID')
        rxnDietID = options.rxnDietID;
    else
        rxnDietID = 'rxnDiet';
    end
end

if strcmp(sense,'uptake')
    origSense = 'uptake';
    nIter = 1;
elseif strcmp(sense,'outside')
    origSense = 'outside';
    nIter = 1;
elseif strcmp(sense,'both')
    origSense = 'uptake';
    oppositeSense = 'outside';
    nIter = 2;
else
    error('The sense is not correct. Options are uptake, outside or both.')
end
    
% Modify the name of the reactions defined as NaN
nan_name = find(cellfun(@(x) x,cellfun(@length,cellfun(@isnan,model.rxns,'UniformOutput',false),'UniformOutput',false))==1);
model.rxns(nan_name) = cellfun(@strcat, repmat({'Name_'},length(nan_name),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(nan_name))'),'UniformOutput',false),'UniformOutput',false);
nan_name = find(cellfun(@(x) x,cellfun(@length,cellfun(@isnan,model.rxnNames,'UniformOutput',false),'UniformOutput',false))==1);
model.rxnNames(nan_name) = cellfun(@strcat, repmat({'Name_'},length(nan_name),1), cellfun(@(x) num2str(x), ...
    num2cell((1:length(nan_name))'),'UniformOutput',false),'UniformOutput',false);

% Define core reactions
model.core = zeros(length(model.rxns),1);
model.core(core) = 1;

% Identify the reactions from plants
if ~isempty(plantRxns)
    model.plantRxns = zeros(length(model.rxns),1);
    model.plantRxns(ismember(model.rxnID,plantRxns)) = 1;
end

model.modelName = modelName;
model_orig = model;

% Split the model in order to apply the fastFVA
revs_f = find(model.lb < 0 & model.ub > 0);
revRxnsRelation = [revs_f (length(model.rxns)+1:length(model.rxns)+length(revs_f))'];

model.S = [model.S, -model.S(:, revs_f)];
model.ub = [model.ub; -model.lb(revs_f)];
model.lb = zeros(size(model.ub));
model.rxns = [model.rxns; model.rxns(revs_f)];
model.c = [model.c; model.c(revs_f)];
model.rules = [model.rules; model.rules(revs_f)];
model.rxnID = [model.rxnID; model.rxnID(revs_f)];

if exist([fastFVAName '.mat'], 'file')
    load([fastFVAName '.mat'], 'vMin', 'vMax')
    if length(vMin)~=length(model.rxns) || length(vMax)~=length(model.rxns)
        error('The fastFVA fluxes do not correspond to the model reactions')
    end
else
    [vMin, vMax] = fastFVA(model);
    save([fastFVAName '.mat'], 'vMin', 'vMax')
end

% Change reaction bounds and delete blocked reactions from the model
epsilon = 1e-06;
blockedRxns = find(abs(vMin)<=epsilon & vMax<=epsilon);

allRevRxns = [revRxnsRelation(:,1);revRxnsRelation(:,2)];
irrevRxnsToDelete = blockedRxns(~ismember(blockedRxns,allRevRxns));
revRxnsToChange = blockedRxns(ismember(blockedRxns,allRevRxns));

for i = 1:length(revRxnsToChange)
    if ismember(revRxnsToChange(i), revRxnsRelation(:,1))
        model_orig.ub(revRxnsToChange(i)) = 0;
    elseif ismember(revRxnsToChange(i), revRxnsRelation(:,2))
        pos = revRxnsRelation(ismember(revRxnsRelation(:,2),revRxnsToChange(i)),1);
        model_orig.lb(pos) = 0;
    end
end

model_orig.S(:,irrevRxnsToDelete) = [];
model_orig.lb(irrevRxnsToDelete) = [];
model_orig.ub(irrevRxnsToDelete) = [];
model_orig.rxns(irrevRxnsToDelete) = [];
model_orig.rxnNames(irrevRxnsToDelete) = [];
model_orig.c(irrevRxnsToDelete) = [];
model_orig.trules(irrevRxnsToDelete) = [];
model_orig.rules(irrevRxnsToDelete) = [];
model_orig.subSystems(irrevRxnsToDelete) = [];
model_orig.grRules(irrevRxnsToDelete) = [];
model_orig.comments(irrevRxnsToDelete) = [];
model_orig.citations(irrevRxnsToDelete) = [];
model_orig.rxnConfidenceScores(irrevRxnsToDelete) = [];
model_orig.correctedRxns(irrevRxnsToDelete) = [];
model_orig.rxnECNumbers(irrevRxnsToDelete) = [];
model_orig.rxnKEGGID(irrevRxnsToDelete) = [];
model_orig.rxnID(irrevRxnsToDelete) = [];
model_orig.rxnCode(irrevRxnsToDelete) = [];
model_orig.rxnStoichiometry(irrevRxnsToDelete) = [];
model_orig.rxnIsTransport(irrevRxnsToDelete) = [];
model_orig.rxnEquation(irrevRxnsToDelete) = [];
model_orig.rxnDefinition(irrevRxnsToDelete) = [];
model_orig.rxnReversibility(irrevRxnsToDelete) = [];
model_orig.rxnDirection(irrevRxnsToDelete) = [];
model_orig.rxnAbstractReaction(irrevRxnsToDelete) = [];
model_orig.rxnPathways(irrevRxnsToDelete) = [];
model_orig.rxnAliases(irrevRxnsToDelete) = [];
model_orig.rxnDeltaG(irrevRxnsToDelete) = [];
model_orig.rxnDeltaGErr(irrevRxnsToDelete) = [];
model_orig.rxnCompoundIDs(irrevRxnsToDelete) = [];
model_orig.rxnStatus(irrevRxnsToDelete) = [];
model_orig.rxnIsObsolete(irrevRxnsToDelete) = [];
model_orig.rxnLinkedReaction(irrevRxnsToDelete) = [];
model_orig.rxnNotes(irrevRxnsToDelete) = [];
model_orig.rxnSource(irrevRxnsToDelete) = [];
model_orig.rxnOntology(irrevRxnsToDelete) = [];
model_orig.trRules(irrevRxnsToDelete) = [];
model_orig.rxnGeneMat(irrevRxnsToDelete,:) = [];
model_orig.rxnTaxMat(irrevRxnsToDelete,:) = [];
model_orig.core(irrevRxnsToDelete) = [];
if isfield(model_orig, 'plantRxns')
    model_orig.plantRxns(irrevRxnsToDelete,:) = [];
end

model = model_orig;
finalModel = model_orig;

% Calculate the i-Diet blocked metabolites
dietRxnsIdx = find(~cellfun(@isempty,regexp(model.rxnID,rxnDietID)));
dietRxns = model.rxnID(dietRxnsIdx);
for j = 1:nIter
    if j == 1
        sense = origSense;
    else
        sense = oppositeSense;
    end
    if strcmp(sense,'uptake')
        if isempty(blockedDietRxnsUp)
            nDiet = length(dietRxnsIdx);
            blockedDietRxnsUp = zeros(nDiet,1);
            for i = 1:nDiet
                modelTmp = model;
                modelTmp.c = zeros(length(modelTmp.c),1);
                modelTmp.c(dietRxnsIdx(i)) = 1;
                cplex = FBA(modelTmp,'minimize');
                if abs(cplex.Solution.objval) < 1e-08
                    blockedDietRxnsUp(i) = 1;
                end
            end
            save(fullfile('.','input','optionsFastcore',['blockedDietRxnsUp_' modelName '.mat']), 'blockedDietRxnsUp')
            dietRxnsOfIntUp = dietRxns(~logical(blockedDietRxnsUp));
        else
            dietRxnsOfIntUp = dietRxns(~logical(blockedDietRxnsUp));
        end
    elseif strcmp(sense,'outside')
        if isempty(blockedDietRxnsOut)
            nDiet = length(dietRxnsIdx);
            blockedDietRxnsOut = zeros(nDiet,1);
            for i = 1:nDiet
                modelTmp = model;
                modelTmp.c = zeros(length(modelTmp.c),1);
                modelTmp.c(dietRxnsIdx(i)) = 1;
                cplex = FBA(modelTmp,'maximize');
                if abs(cplex.Solution.objval) < 1e-08
                    blockedDietRxnsOut(i) = 1;
                end
            end
            save(fullfile('.','input','optionsFastcore',['blockedDietRxnsOut_' modelName '.mat']), 'blockedDietRxnsOut')
            dietRxnsOfIntOut = dietRxns(~logical(blockedDietRxnsOut));
        else
            dietRxnsOfIntOut = dietRxns(~logical(blockedDietRxnsOut));
        end
    end
end

% Split the model
nOrigRxns = length(model.rxns);
origCore = find(model.core);
revs_f = find(model.lb < 0 & model.ub > 0);
revRxnsRelation = [revs_f (length(model.rxns)+1:length(model.rxns)+length(revs_f))'];

model.S = [model.S, -model.S(:, revs_f)];
model.ub = [model.ub; -model.lb(revs_f)];
model.lb = zeros(size(model.ub));
model.rxns = [model.rxns; model.rxns(revs_f)];
model.c = [model.c; model.c(revs_f)];
model.rules = [model.rules; model.rules(revs_f)];
model.trules = [model.trules; model.trules(revs_f)];
model.rxnID = [model.rxnID; model.rxnID(revs_f)];
model.core = [model.core; model.core(revs_f)];
if isfield(model, 'plantRxns')
    model.plantRxns = [model.plantRxns; model.plantRxns(revs_f)];
end

% Apply fastcore with the actual core
weights = ones(length(model.rxns),1)*100;
weights(~cellfun(@isempty,model.rxnECNumbers)) = 50;
if isfield(model, 'plantRxns')
    weights(logical(model.plantRxns)) = 1000;
end
weights(~cellfun(@isempty,model.trules)) = 0.1;
weights(logical(model.core)) = 0;

fastcoreResult = fastCoreWeighted(find(model.core),model,weights);
coreRxnBoolSplitted = false(length(model.rxns),1);
coreRxnBoolSplitted(fastcoreResult) = true;

firstSenseRxnBool = coreRxnBoolSplitted(1:nOrigRxns);
secondSenseRxnBool = false(nOrigRxns,1);
tmp = fastcoreResult(fastcoreResult>nOrigRxns);
tmp = revRxnsRelation(ismember(revRxnsRelation(:,2),tmp),1);
secondSenseRxnBool(tmp) = true;
coreRxnBoolOrig = firstSenseRxnBool | secondSenseRxnBool;

rxnsFound = find(coreRxnBoolOrig);
coreAddedRxnsIdx = rxnsFound(~ismember(rxnsFound,origCore));
newCore = [origCore; coreAddedRxnsIdx];
model.core(fastcoreResult) = 1;

coreAddedRxns = cell(length(coreAddedRxnsIdx),6);
coreAddedRxns(:,1) = model.rxnID(coreAddedRxnsIdx);
coreAddedRxns(:,2) = model.rxnNames(coreAddedRxnsIdx);
coreAddedRxns(:,3) = model.rxnAliases(coreAddedRxnsIdx);
coreAddedRxns(:,4) = printRxnFormula(model,'rxnAbbrList',model.rxns(coreAddedRxnsIdx),'metNameFlag',true);
clc
coreAddedRxns(:,5) = model.rxnECNumbers(coreAddedRxnsIdx);
coreAddedRxns(:,6) = model.trules(coreAddedRxnsIdx);

% Identify the reactions for i-Diet metabolites
dietRxnsIdx = find(~cellfun(@isempty,regexp(model.rxnID,rxnDietID)));

for j = 1:nIter
    if j == 1
        sense = origSense;
    else
        sense = oppositeSense;
    end

    if strcmp(sense,'uptake')
        dietRxnsUpIdx = dietRxnsIdx(sum(model.S(:,dietRxnsIdx), 1) == 1);
        dietRxnsOutIdx = dietRxnsIdx(sum(model.S(:,dietRxnsIdx), 1) == -1);
        dietRxnsUp = model.rxnID(dietRxnsUpIdx);
        dietRxnsOut = model.rxnID(dietRxnsOutIdx);
        dietRxnsUpIdx = dietRxnsUpIdx(ismember(dietRxnsUp,dietRxnsOfIntUp));
        dietRxnsOutIdx = dietRxnsOutIdx(ismember(dietRxnsOut,dietRxnsOfIntUp));
        nDiet = length(dietRxnsOfIntUp);
    else
        dietRxnsUpIdx = dietRxnsIdx(sum(model.S(:,dietRxnsIdx), 1) == 1);
        dietRxnsOutIdx = dietRxnsIdx(sum(model.S(:,dietRxnsIdx), 1) == -1);
        dietRxnsUp = model.rxnID(dietRxnsUpIdx);
        dietRxnsOut = model.rxnID(dietRxnsOutIdx);
        dietRxnsUpIdx = dietRxnsUpIdx(ismember(dietRxnsUp,dietRxnsOfIntOut));
        dietRxnsOutIdx = dietRxnsOutIdx(ismember(dietRxnsOut,dietRxnsOfIntOut));
        nDiet = length(dietRxnsOfIntOut);
    end
    iDietAddedRxns = cell(nDiet,8);
    for i = 1:nDiet
        clc
        disp([num2str(i) '/' num2str(nDiet) ' sense:' sense]);
        
        % If the iDiet reaction was already in the core, avoid the fastcore
        % step
        if strcmp(sense,'uptake')
            if model.core(dietRxnsUpIdx(i)) == 1
                iDietAddedRxns{i,1} = model.rxns(dietRxnsUpIdx(i));
                iDietAddedRxns{i,8} = {'CONNECTED'};
                continue;
            end
        elseif strcmp(sense,'outside')
            if model.core(dietRxnsOutIdx(i)) == 1
                iDietAddedRxns{i,1} = model.rxns(dietRxnsOutIdx(i));
                iDietAddedRxns{i,8} = {'CONNECTED'};
                continue;
            end
        end

        weights = ones(length(model.rxns),1)*100;
        weights(~cellfun(@isempty,model.rxnECNumbers)) = 30;
        if isfield(model, 'plantRxns')
            weights(logical(model.plantRxns)) = 1000;
        end
        weights(~cellfun(@isempty,model.trules)) = 0.1;
        if strcmp(sense,'uptake')
            model.core(dietRxnsUpIdx(i)) = 1;
            weights(dietRxnsOutIdx(i),1) = 100000;
        elseif strcmp(sense,'outside')
            model.core(dietRxnsOutIdx(i)) = 1;
            weights(dietRxnsUpIdx(i),1) = 100000;
        end
        weights(logical(model.core)) = 0;

        fastcoreResult = fastCoreWeighted(find(model.core),model,weights);
        coreRxnBoolSplitted = false(length(model.rxns),1);
        coreRxnBoolSplitted(fastcoreResult) = true;

        firstSenseRxnBool = coreRxnBoolSplitted(1:nOrigRxns);
        secondSenseRxnBool = false(nOrigRxns,1);
        tmp = fastcoreResult(fastcoreResult>nOrigRxns);
        tmp = revRxnsRelation(ismember(revRxnsRelation(:,2),tmp),1);
        secondSenseRxnBool(tmp) = true;
        coreRxnBoolOrig = firstSenseRxnBool | secondSenseRxnBool;

        rxnsFound = find(coreRxnBoolOrig);
        iDietAddedRxnsIdx = rxnsFound(~ismember(rxnsFound,newCore));

        % Define the reactions added by each iDiet metabolite
        iDietAddedRxns{i,1} = model.rxns(dietRxnsOutIdx(i));
        iDietAddedRxns{i,2} = model.rxnID(iDietAddedRxnsIdx);
        iDietAddedRxns{i,3} = model.rxns(iDietAddedRxnsIdx);
        iDietAddedRxns{i,4} = model.rxnAliases(iDietAddedRxnsIdx);
        iDietAddedRxns{i,5} = printRxnFormula(model,'rxnAbbrList',model.rxns(iDietAddedRxnsIdx),'metNameFlag',true);
        clc
        iDietAddedRxns{i,6} = model.rxnECNumbers(iDietAddedRxnsIdx);
        iDietAddedRxns{i,7} = model.trules(iDietAddedRxnsIdx);
        rxnsWithoutTax = iDietAddedRxnsIdx(cellfun(@isempty,model.trules(iDietAddedRxnsIdx)));
        admittedRxns = (~cellfun(@isempty,regexp(model.rxnID(rxnsWithoutTax),'rxnDiet'))) | (~cellfun(@isempty,regexp(model.rxnID(rxnsWithoutTax),'rxnAdd')));
        if sum(admittedRxns) == length(rxnsWithoutTax)
            iDietAddedRxns{i,8} = {'CONNECTED'};
        else
            iDietAddedRxns{i,8} = {'NOT CONNECTED'};
        end
        
        if strcmp(sense,'uptake')
            model.core(dietRxnsUpIdx(i)) = 0;
        elseif strcmp(sense,'outside')
            model.core(dietRxnsOutIdx(i)) = 0;
        end      
    end
    if j == 1
        iDietAddedRxnsTmp = iDietAddedRxns;
    end
end

% Build the final model
if nIter == 2
    iDietAddedRxns = [iDietAddedRxnsTmp; iDietAddedRxns];
end
newRxns = unique(cat(1,iDietAddedRxns{:,2}));
newRxnsIdx = find(ismember(finalModel.rxnID,newRxns));
rxnsToDelete = true(length(finalModel.rxns),1);
rxnsToDelete([newCore; newRxnsIdx]) = false;

finalModel.S(:,rxnsToDelete) = [];
finalModel.rxns(rxnsToDelete) = [];
finalModel.lb(rxnsToDelete) = [];
finalModel.ub(rxnsToDelete) = [];
finalModel.c(rxnsToDelete) = [];
finalModel.rules(rxnsToDelete) = [];
finalModel.rxnNames(rxnsToDelete) = [];
finalModel.subSystems(rxnsToDelete) = [];
finalModel.grRules(rxnsToDelete) = [];
finalModel.comments(rxnsToDelete) = [];
finalModel.citations(rxnsToDelete) = [];
finalModel.rxnConfidenceScores(rxnsToDelete) = [];
finalModel.correctedRxns(rxnsToDelete) = [];
finalModel.rxnECNumbers(rxnsToDelete) = [];
finalModel.rxnKEGGID(rxnsToDelete) = [];
finalModel.rxnID(rxnsToDelete) = [];
finalModel.rxnCode(rxnsToDelete) = [];
finalModel.rxnStoichiometry(rxnsToDelete) = [];
finalModel.rxnIsTransport(rxnsToDelete) = [];
finalModel.rxnEquation(rxnsToDelete) = [];
finalModel.rxnDefinition(rxnsToDelete) = [];
finalModel.rxnReversibility(rxnsToDelete) = [];
finalModel.rxnDirection(rxnsToDelete) = [];
finalModel.rxnAbstractReaction(rxnsToDelete) = [];
finalModel.rxnPathways(rxnsToDelete) = [];
finalModel.rxnAliases(rxnsToDelete) = [];
finalModel.rxnDeltaG(rxnsToDelete) = [];
finalModel.rxnDeltaGErr(rxnsToDelete) = [];
finalModel.rxnCompoundIDs(rxnsToDelete) = [];
finalModel.rxnStatus(rxnsToDelete) = [];
finalModel.rxnIsObsolete(rxnsToDelete) = [];
finalModel.rxnLinkedReaction(rxnsToDelete) = [];
finalModel.rxnNotes(rxnsToDelete) = [];
finalModel.rxnSource(rxnsToDelete) = [];
finalModel.rxnOntology(rxnsToDelete) = [];
finalModel.trules(rxnsToDelete) = [];
finalModel.trRules(rxnsToDelete) = [];
finalModel.rxnGeneMat(rxnsToDelete,:) = [];
finalModel.rxnTaxMat(rxnsToDelete,:) = [];

metsToDelete = find(sum(finalModel.S~=0, 2)==0);

finalModel.S(metsToDelete,:) = [];
finalModel.mets(metsToDelete) = [];
finalModel.b(metsToDelete) = [];
finalModel.csense(metsToDelete) = [];
finalModel.metNames(metsToDelete) = [];
finalModel.metCharges(metsToDelete) = [];
finalModel.metFormulas(metsToDelete) = [];
finalModel.metChEBIID(metsToDelete) = [];
finalModel.metKEGGID(metsToDelete) = [];
finalModel.metPubChemID(metsToDelete) = [];
finalModel.metInChIString(metsToDelete) = [];
finalModel.metHMDBID(metsToDelete) = [];
finalModel.metSmiles(metsToDelete) = [];
finalModel.metID(metsToDelete) = [];
finalModel.metsTax(metsToDelete) = [];
finalModel.metMass(metsToDelete) = [];
finalModel.metSource(metsToDelete) = [];
finalModel.metInchiKey(metsToDelete) = [];
finalModel.metIsCore(metsToDelete) = [];
finalModel.metIsObsolete(metsToDelete) = [];
finalModel.metLinkedCompound(metsToDelete) = [];
finalModel.metIsCofactor(metsToDelete) = [];
finalModel.metDeltaG(metsToDelete) = [];
finalModel.metDeltaGErr(metsToDelete) = [];
finalModel.metPKA(metsToDelete) = [];
finalModel.metPKB(metsToDelete) = [];
finalModel.metAbstractCompound(metsToDelete) = [];
finalModel.metComprisedOf(metsToDelete) = [];
finalModel.metAliases(metsToDelete) = [];
finalModel.metOntology(metsToDelete) = [];

% Remove unnecessary fields
finalModel = rmfield(finalModel,'core');
if isfield(finalModel,'plantRxns')
    finalModel = rmfield(finalModel,'plantRxns');
end

end