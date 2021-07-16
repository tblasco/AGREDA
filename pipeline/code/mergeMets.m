function model = mergeMets(model,metsMan,metsRem)

% mergeMets.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

metsToDelete = [];
for i = 1:length(metsMan)
    mMan = find(ismember(model.metID,metsMan(i)));
    mRem = find(ismember(model.metID,metsRem(i)));
    relRxns = find(model.S(mRem,:));
    
    for j = 1:length(relRxns)
        model.S(mMan,relRxns(j)) = model.S(mRem,relRxns(j));
        model.S(mRem,relRxns(j)) = 0;
    end
    
    metsToDelete = [metsToDelete, mRem];
    
end
    
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