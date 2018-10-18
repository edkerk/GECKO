%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,P,C,R,D,L,M,W] = sumBioMass(model)
%
% Eduard Kerkhoven, 2018-10-10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,P,C,R,D,L,M,W] = sumBioMass(model)
%Get main fractions:
[P,X] = getFraction(model,'protein',0);
[C,X] = getFraction(model,'carbohydrate',X);
[R,X] = getFraction(model,'RNA',X);
[D,X] = getFraction(model,'DNA',X);
[L,X] = getFraction(model,'lipid',X);
[M,X] = getFraction(model,'misc',X);
[W,X] = getFraction(model,'cell wall',X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F,X] = getFraction(model,compType,X)
%Define pseudoreaction name:
rxnName = [compType ' pseudoreaction'];

%Add up fraction:
fractionPos = strcmp(model.rxnNames,rxnName);
%fractionPos=getIndexes(model,'BIOMASS_SCO','rxns');
mets = find(model.S(:,fractionPos) ~= 0);
formulas = model.metFormulas(mets);
coeffs = model.S(mets,fractionPos);
F = 0;
[~, ~, ~, MW] = parseFormulas(formulas);

% Pseudometabolites don't have MW, and should be ignored
NaNidx = find(isnan(MW));
if ~isempty(NaNidx)
    NaNname = model.metNames(mets(NaNidx));
    for i = 1:length(NaNidx)
        if contains(NaNname(i),'pseudometabolite')
            MW(NaNidx(i))=[];
            coeffs(NaNidx(i))=[];
        else
            error(['Chemical formulae are not defined for all biomass '...
                'components.']);
        end
    end
end
F = sum(-coeffs .* MW)/1000;
X = X + F;
end