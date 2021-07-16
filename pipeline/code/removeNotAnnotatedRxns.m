function model = removeNotAnnotatedRxns(model,fastFVAName)

% removeNotAnnotatedRxns.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es
    
    if nargin == 1
        fastFVAName = 'removeNotAnnotatedRxnsFastFVA';
    end

    rxns_without_tax_allowed = find(~cellfun(@isempty,regexp(model.rxnID,'rxnDiet')) | ~cellfun(@isempty,regexp(model.rxnID,'rxnAdd')));
    rxns_without_tax = find(cellfun(@isempty,model.trules));

    rxns_to_delete = rxns_without_tax(~ismember(rxns_without_tax,rxns_without_tax_allowed));

    model.S(:,rxns_to_delete) = [];
    model.rxns(rxns_to_delete) = [];
    model.lb(rxns_to_delete) = [];
    model.ub(rxns_to_delete) = [];
    model.c(rxns_to_delete) = [];
    model.rules(rxns_to_delete) = [];
    model.rxnNames(rxns_to_delete) = [];
    model.subSystems(rxns_to_delete) = [];
    model.grRules(rxns_to_delete) = [];
    model.comments(rxns_to_delete) = [];
    model.citations(rxns_to_delete) = [];
    model.rxnConfidenceScores(rxns_to_delete) = [];
    model.correctedRxns(rxns_to_delete) = [];
    model.rxnECNumbers(rxns_to_delete) = [];
    model.rxnKEGGID(rxns_to_delete) = [];
    model.rxnID(rxns_to_delete) = [];
    model.rxnCode(rxns_to_delete) = [];
    model.rxnStoichiometry(rxns_to_delete) = [];
    model.rxnIsTransport(rxns_to_delete) = [];
    model.rxnEquation(rxns_to_delete) = [];
    model.rxnDefinition(rxns_to_delete) = [];
    model.rxnReversibility(rxns_to_delete) = [];
    model.rxnDirection(rxns_to_delete) = [];
    model.rxnAbstractReaction(rxns_to_delete) = [];
    model.rxnPathways(rxns_to_delete) = [];
    model.rxnAliases(rxns_to_delete) = [];
    model.rxnDeltaG(rxns_to_delete) = [];
    model.rxnDeltaGErr(rxns_to_delete) = [];
    model.rxnCompoundIDs(rxns_to_delete) = [];
    model.rxnStatus(rxns_to_delete) = [];
    model.rxnIsObsolete(rxns_to_delete) = [];
    model.rxnLinkedReaction(rxns_to_delete) = [];
    model.rxnNotes(rxns_to_delete) = [];
    model.rxnSource(rxns_to_delete) = [];
    model.rxnOntology(rxns_to_delete) = [];
    model.trules(rxns_to_delete) = [];
    model.trRules(rxns_to_delete) = [];
    model.rxnGeneMat(rxns_to_delete,:) = [];
    model.rxnTaxMat(rxns_to_delete,:) = [];

    metsToDelete = find(sum(model.S~=0, 2)==0);

    model.S(metsToDelete,:) = [];
    model.mets(metsToDelete) = [];
    model.b(metsToDelete) = [];
    model.csense(metsToDelete) = [];
    model.metNames(metsToDelete) = [];
    model.metCharges(metsToDelete) = [];
    model.metFormulas(metsToDelete) = [];
    model.metChEBIID(metsToDelete) = [];
    model.metKEGGID(metsToDelete) = [];
    model.metPubChemID(metsToDelete) = [];
    model.metInChIString(metsToDelete) = [];
    model.metHMDBID(metsToDelete) = [];
    model.metSmiles(metsToDelete) = [];
    model.metID(metsToDelete) = [];
    model.metsTax(metsToDelete) = [];
    model.metMass(metsToDelete) = [];
    model.metSource(metsToDelete) = [];
    model.metInchiKey(metsToDelete) = [];
    model.metIsCore(metsToDelete) = [];
    model.metIsObsolete(metsToDelete) = [];
    model.metLinkedCompound(metsToDelete) = [];
    model.metIsCofactor(metsToDelete) = [];
    model.metDeltaG(metsToDelete) = [];
    model.metDeltaGErr(metsToDelete) = [];
    model.metPKA(metsToDelete) = [];
    model.metPKB(metsToDelete) = [];
    model.metAbstractCompound(metsToDelete) = [];
    model.metComprisedOf(metsToDelete) = [];
    model.metAliases(metsToDelete) = [];
    model.metOntology(metsToDelete) = [];
    
    if exist([fastFVAName '.mat'], 'file')
        load([fastFVAName '.mat'], 'vMin', 'vMax')
        if length(vMin)~=length(model.rxns) || length(vMax)~=length(model.rxns)
            error('The fastFVA fluxes do not correspond to the model reactions')
        end
    else
        [vMin, vMax] = fastFVA(model);
        save([fastFVAName '.mat'], 'vMin', 'vMax')
    end
    
    [vMin, vMax] = fastFVA(model);
    epsilon = 1e-08;
    rxns_to_delete = find(abs(vMax)<epsilon & abs(vMin)<epsilon);

    model.S(:,rxns_to_delete) = [];
    model.rxns(rxns_to_delete) = [];
    model.lb(rxns_to_delete) = [];
    model.ub(rxns_to_delete) = [];
    model.c(rxns_to_delete) = [];
    model.rules(rxns_to_delete) = [];
    model.rxnNames(rxns_to_delete) = [];
    model.subSystems(rxns_to_delete) = [];
    model.grRules(rxns_to_delete) = [];
    model.comments(rxns_to_delete) = [];
    model.citations(rxns_to_delete) = [];
    model.rxnConfidenceScores(rxns_to_delete) = [];
    model.correctedRxns(rxns_to_delete) = [];
    model.rxnECNumbers(rxns_to_delete) = [];
    model.rxnKEGGID(rxns_to_delete) = [];
    model.rxnID(rxns_to_delete) = [];
    model.rxnCode(rxns_to_delete) = [];
    model.rxnStoichiometry(rxns_to_delete) = [];
    model.rxnIsTransport(rxns_to_delete) = [];
    model.rxnEquation(rxns_to_delete) = [];
    model.rxnDefinition(rxns_to_delete) = [];
    model.rxnReversibility(rxns_to_delete) = [];
    model.rxnDirection(rxns_to_delete) = [];
    model.rxnAbstractReaction(rxns_to_delete) = [];
    model.rxnPathways(rxns_to_delete) = [];
    model.rxnAliases(rxns_to_delete) = [];
    model.rxnDeltaG(rxns_to_delete) = [];
    model.rxnDeltaGErr(rxns_to_delete) = [];
    model.rxnCompoundIDs(rxns_to_delete) = [];
    model.rxnStatus(rxns_to_delete) = [];
    model.rxnIsObsolete(rxns_to_delete) = [];
    model.rxnLinkedReaction(rxns_to_delete) = [];
    model.rxnNotes(rxns_to_delete) = [];
    model.rxnSource(rxns_to_delete) = [];
    model.rxnOntology(rxns_to_delete) = [];
    model.trules(rxns_to_delete) = [];
    model.trRules(rxns_to_delete) = [];
    model.rxnGeneMat(rxns_to_delete,:) = [];
    model.rxnTaxMat(rxns_to_delete,:) = [];

    metsToDelete = find(sum(model.S~=0, 2)==0);

    model.S(metsToDelete,:) = [];
    model.mets(metsToDelete) = [];
    model.b(metsToDelete) = [];
    model.csense(metsToDelete) = [];
    model.metNames(metsToDelete) = [];
    model.metCharges(metsToDelete) = [];
    model.metFormulas(metsToDelete) = [];
    model.metChEBIID(metsToDelete) = [];
    model.metKEGGID(metsToDelete) = [];
    model.metPubChemID(metsToDelete) = [];
    model.metInChIString(metsToDelete) = [];
    model.metHMDBID(metsToDelete) = [];
    model.metSmiles(metsToDelete) = [];
    model.metID(metsToDelete) = [];
    model.metsTax(metsToDelete) = [];
    model.metMass(metsToDelete) = [];
    model.metSource(metsToDelete) = [];
    model.metInchiKey(metsToDelete) = [];
    model.metIsCore(metsToDelete) = [];
    model.metIsObsolete(metsToDelete) = [];
    model.metLinkedCompound(metsToDelete) = [];
    model.metIsCofactor(metsToDelete) = [];
    model.metDeltaG(metsToDelete) = [];
    model.metDeltaGErr(metsToDelete) = [];
    model.metPKA(metsToDelete) = [];
    model.metPKB(metsToDelete) = [];
    model.metAbstractCompound(metsToDelete) = [];
    model.metComprisedOf(metsToDelete) = [];
    model.metAliases(metsToDelete) = [];
    model.metOntology(metsToDelete) = [];
end