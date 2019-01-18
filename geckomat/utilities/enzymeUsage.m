function [enzUsage,protein,reactionIdx] = enzymeUsage(ecModel,fluxes,nonzero)
% enzymeUsage
%
%   Calculate enzyme usage
%
%   Input:
%   ecModel         enzyme-constrained model
%   fluxes          vector of fluxes, for instance sol.x from sol=solveLP
%   nonzero         logical whether nonzero enzyme usages should be
%                   included (opt, default true)
%
%   Output:
%   enzUsage        ratio of enzyme usage
%   protein         array of matching UniProt protein IDs      
%   reactionIdx     vector of indexes of exchange reactions
%
% Usage: [enzUsage,protein,reactionIdx] = enzymeUsage(ecModel,fluxes)
%
% Eduard Kerkhoven  Last edited: 2018-12-19

if nargin<3
    nonzero=true;
end

protIdx     = contains(ecModel.rxnNames,'prot_');
enzUsage    = fluxes(protIdx)./ecModel.ub(protIdx);
protein     = regexprep(ecModel.rxnNames(protIdx),'(draw_)?prot_','');
protein     = regexprep(protein,'_exchange','');
reactionIdx = find(protIdx);

if ~nonzero
    nonzero     = enzUsage>0;
    enzUsage    = enzUsage(nonzero);
    protein     = protein(nonzero);
    reactionIdx = reactionIdx(nonzero);
end
end
