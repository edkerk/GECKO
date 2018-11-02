%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = scaleBioMass(model,Ptot,GAM,scale_comp)
% 
% Benjamin Sanchez. Last update: 2018-10-23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = scaleBioMass(model,Ptot,GAM,scale_comp)

if nargin < 3
    GAM = [];
end

%Option for changing composition & GAM (=true, default) or only GAM (=false):
if nargin < 4
    scale_comp = true;
end

%Compute carbohydrate and lipid new amounts, based on:
%1. Total mass remains constant, i.e. Pbase+Cbase+Lbase+Rbase =
%   Ptot+Ctot+Ltot+Rtot
%2. Difference in mass is distributed proportionally, i.e. Ctot/Ltot = Cbase/Lbase
% [X,P,C,R,D,L,M,W] = sumBioMass(model)
[~,Pbase,Cbase,Rbase,Dtot,Lbase,~,Wtot] = sumBioMass(model);

Ctot = Cbase + (Pbase - Ptot) * (Cbase / (Cbase+Rbase+Lbase));
Rtot = Rbase + (Pbase - Ptot) * (Rbase / (Cbase+Rbase+Lbase));
Ltot = Lbase + (Pbase - Ptot) * (Lbase / (Cbase+Rbase+Lbase));

%Compute rescaling fractions:
fP = Ptot/Pbase;
fC = Ctot/Cbase;
fR = Rtot/Rbase;
fL = Ltot/Lbase;

%Change compositions:
if scale_comp
    model = rescalePseudoReaction(model,'protein',fP);
    model = rescalePseudoReaction(model,'carbohydrate',fC);
    model = rescalePseudoReaction(model,'lipid',fL);
    model = rescalePseudoReaction(model,'RNA',fR);
end

%Fit GAM if not available:
% if isempty(GAM)
%     GAM = fitGAM(model);
% end
% 
% %Change GAM:
% rxnIdx = getIndexes(model,'BM_growth','rxns');
% metIdx = getIndexes(model,{'atp_c','adp_c','h2o_c','h_c','pi_c'},'mets');
% %Polymerization costs from Borodina et al 2005
% GAM = GAM + 40*Ptot + 4.4*Dtot + 1.25*Rtot + 5.026*Wtot;
% model.S(metIdx, rxnIdx) = sign(model.S(metIdx, rxnIdx)) .* GAM;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = rescalePseudoReaction(model,metName,f)

rxnPos  = strcmp(model.rxnNames,[metName ' pseudoreaction']);

mets = find(model.S(:,rxnPos));
pseudoMet = find(strcmp(model.metNames(mets),[metName ' pseudometabolite']));
mets(pseudoMet)=[]; % Don't scale pseudometabolites

model.S(mets,rxnPos) = model.S(mets,rxnPos)*f;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
