function cplex = FBA(Network,sense)

% FBA.m
% 
% Author: Telmo Blasco, Francesco Balzerani, Iñigo Apaolaza, Francisco J.
% Planes
% Email: tblasco@tecnun.es, fbalzerani@tecnun.es, iaemparanza@tecnun.es,
% fplanes@tecnun.es

    if nargin == 1
        sense = 'maximize';
    end
    if ~(strcmp(sense,'maximize') || strcmp(sense,'minimize'))
        error('The value given to sense is not correct. Options are maximize or minimize.')
    end
    
    % Define the number of constraints and variables
    nCon = size(Network.S,1);
    nVar = size(Network.S,2);
    
    % Find the biomass reaction
    vbio = find(Network.c==1);

    % Define the variables for Cplex
    A = spalloc(nCon,nVar,sum(sum(Network.S~=0)));
    lhs = zeros(nCon,1);
    rhs = zeros(nCon,1);
    ub = zeros(nVar,1);
    lb = zeros(nVar,1);
    obj = zeros(nVar,1);
    A(:) = Network.S;
    ub(:) = Network.ub;
    lb(:) = Network.lb;
    ctype(1:nVar) = 'C';
    obj(vbio)=1;
    
    % Solve the optimization problem
    cplex = Cplex('FBA');
    cplex.Model.A = A;
    cplex.Model.lhs = lhs;
    cplex.Model.rhs = rhs;
    cplex.Model.ub = ub;
    cplex.Model.lb = lb;
    cplex.Model.sense = sense;
    cplex.Model.ctype = ctype;
    cplex.Model.obj = obj;
    cplex.solve()
    clear ctype
end