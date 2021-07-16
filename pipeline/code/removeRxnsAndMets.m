function model = removeRxnsAndMets(model, rxnsToDelete)

% removeRxnsAndMets.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

    if isfield(model,'S')
        model.S(:,rxnsToDelete) = [];
    end
    if isfield(model,'rxns')
        model.rxns(rxnsToDelete) = [];
    end
    if isfield(model,'lb')
        model.lb(rxnsToDelete) = [];
    end
    if isfield(model,'ub')
        model.ub(rxnsToDelete) = [];
    end
    if isfield(model,'c')
        model.c(rxnsToDelete) = [];
    end
    if isfield(model,'rules')
        model.rules(rxnsToDelete) = [];
    end
    if isfield(model,'rxnNames')
        model.rxnNames(rxnsToDelete) = [];
    end
    if isfield(model,'subSystems')
        model.subSystems(rxnsToDelete) = [];
    end
    if isfield(model,'grRules')
        model.grRules(rxnsToDelete) = [];
    end
    if isfield(model,'comments')
        model.comments(rxnsToDelete) = [];
    end
    if isfield(model,'citations')
        model.citations(rxnsToDelete) = [];
    end
    if isfield(model,'rxnConfidenceScores')
        model.rxnConfidenceScores(rxnsToDelete) = [];
    end
    if isfield(model,'correctedRxns')
        model.correctedRxns(rxnsToDelete) = [];
    end
    if isfield(model,'rxnECNumbers')
        model.rxnECNumbers(rxnsToDelete) = [];
    end
    if isfield(model,'rxnKEGGID')
        model.rxnKEGGID(rxnsToDelete) = [];
    end
    if isfield(model,'rxnID')
        model.rxnID(rxnsToDelete) = [];
    end
    if isfield(model,'rxnCode')
        model.rxnCode(rxnsToDelete) = [];
    end
    if isfield(model,'rxnStoichiometry')
        model.rxnStoichiometry(rxnsToDelete) = [];
    end
    if isfield(model,'rxnIsTransport')
        model.rxnIsTransport(rxnsToDelete) = [];
    end
    if isfield(model,'rxnEquation')
        model.rxnEquation(rxnsToDelete) = [];
    end
    if isfield(model,'rxnDefinition')
        model.rxnDefinition(rxnsToDelete) = [];
    end
    if isfield(model,'rxnReversibility')
        model.rxnReversibility(rxnsToDelete) = [];
    end
    if isfield(model,'rxnDirection')
        model.rxnDirection(rxnsToDelete) = [];
    end
    if isfield(model,'rxnAbstractReaction')
        model.rxnAbstractReaction(rxnsToDelete) = [];
    end
    if isfield(model,'rxnPathways')
        model.rxnPathways(rxnsToDelete) = [];
    end
    if isfield(model,'rxnAliases')
        model.rxnAliases(rxnsToDelete) = [];
    end
    if isfield(model,'rxnDeltaG')
        model.rxnDeltaG(rxnsToDelete) = [];
    end
    if isfield(model,'rxnDeltaGErr')
        model.rxnDeltaGErr(rxnsToDelete) = [];
    end
    if isfield(model,'rxnCompoundIDs')
        model.rxnCompoundIDs(rxnsToDelete) = [];
    end
    if isfield(model,'rxnStatus')
        model.rxnStatus(rxnsToDelete) = [];
    end
    if isfield(model,'rxnIsObsolete')
        model.rxnIsObsolete(rxnsToDelete) = [];
    end
    if isfield(model,'rxnLinkedReaction')
        model.rxnLinkedReaction(rxnsToDelete) = [];
    end
    if isfield(model,'rxnNotes')
        model.rxnNotes(rxnsToDelete) = [];
    end
    if isfield(model,'rxnSource')
        model.rxnSource(rxnsToDelete) = [];
    end
    if isfield(model,'rxnOntology')
        model.rxnOntology(rxnsToDelete) = [];
    end
    if isfield(model,'trules')
        model.trules(rxnsToDelete) = [];
    end
    if isfield(model,'trRules')
        model.trRules(rxnsToDelete) = [];
    end
    if isfield(model,'rxnGeneMat')
        model.rxnGeneMat(rxnsToDelete,:) = [];
    end
    if isfield(model,'rxnTaxMat')
        model.rxnTaxMat(rxnsToDelete,:) = [];
    end

    metsToDelete = find(sum(model.S~=0, 2)==0);

    if isfield(model,'S')
        model.S(metsToDelete,:) = [];
    end
    if isfield(model,'mets')
        model.mets(metsToDelete) = [];
    end
    if isfield(model,'b')
        model.b(metsToDelete) = [];
    end
    if isfield(model,'csense')
        model.csense(metsToDelete) = [];
    end
    if isfield(model,'metNames')
        model.metNames(metsToDelete) = [];
    end
    if isfield(model,'metCharges')
        model.metCharges(metsToDelete) = [];
    end
    if isfield(model,'metFormulas')
        model.metFormulas(metsToDelete) = [];
    end
    if isfield(model,'metChEBIID')
        model.metChEBIID(metsToDelete) = [];
    end
    if isfield(model,'metKEGGID')
        model.metKEGGID(metsToDelete) = [];
    end
    if isfield(model,'metPubChemID')
        model.metPubChemID(metsToDelete) = [];
    end
    if isfield(model,'metInChIString')
        model.metInChIString(metsToDelete) = [];
    end
    if isfield(model,'metHMDBID')
        model.metHMDBID(metsToDelete) = [];
    end
    if isfield(model,'metSmiles')
        model.metSmiles(metsToDelete) = [];
    end
    if isfield(model,'metID')
        model.metID(metsToDelete) = [];
    end
    if isfield(model,'metsTax')
        model.metsTax(metsToDelete) = [];
    end
    if isfield(model,'metMass')
        model.metMass(metsToDelete) = [];
    end
    if isfield(model,'metSource')
        model.metSource(metsToDelete) = [];
    end
    if isfield(model,'metInchiKey')
        model.metInchiKey(metsToDelete) = [];
    end
    if isfield(model,'metIsCore')
        model.metIsCore(metsToDelete) = [];
    end
    if isfield(model,'metIsObsolete')
        model.metIsObsolete(metsToDelete) = [];
    end
    if isfield(model,'metLinkedCompound')
        model.metLinkedCompound(metsToDelete) = [];
    end
    if isfield(model,'metIsCofactor')
        model.metIsCofactor(metsToDelete) = [];
    end
    if isfield(model,'metDeltaG')
        model.metDeltaG(metsToDelete) = [];
    end
    if isfield(model,'metDeltaGErr')
        model.metDeltaGErr(metsToDelete) = [];
    end
    if isfield(model,'metPKA')
        model.metPKA(metsToDelete) = [];
    end
    if isfield(model,'metPKB')
        model.metPKB(metsToDelete) = [];
    end
    if isfield(model,'metAbstractCompound')
        model.metAbstractCompound(metsToDelete) = [];
    end
    if isfield(model,'metComprisedOf')
        model.metComprisedOf(metsToDelete) = [];
    end
    if isfield(model,'metAliases')
        model.metAliases(metsToDelete) = [];
    end
    if isfield(model,'metOntology')
        model.metOntology(metsToDelete) = [];
    end 
end