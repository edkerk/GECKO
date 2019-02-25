function [capUsage,absUsage] = enzymeUsage(ecModel,fluxes,zero)
% enzymeUsage
%
%   Calculate enzyme usage
%
%   Input:
%   ecModel         enzyme-constrained model
%   fluxes          vector of fluxes, for instance sol.x from sol=solveLP
%   zero            logical whether also zero enzyme usages should be
%                   included (opt, default true)
%
%   Output:
%   capUsage        ratio of enzyme usage per available enzyme, i.e.
%                   capacity usage (flux divided by ub)
%   absUsage        absolute enzyme usage (flux)
%
% Note: the order of the enzymes in capUsage and absUsage is identical to
% the ecModel.enzGenes and ecModel.enzymes that are used as inputs.
%
% Usage: [capUsage,absUsage] = enzymeUsage(ecModel,fluxes,zero)
%
% Eduard Kerkhoven  Last edited: 2019-01-25

if nargin<3
    zero=true;
end

protIdx     = find(contains(ecModel.rxnNames,ecModel.enzymes));
matchProt   = regexprep(ecModel.rxnNames(protIdx),'(draw_)?prot_','');
matchProt   = regexprep(matchProt,'_exchange','');
[~,b]       = ismember(ecModel.enzymes, matchProt);

absUsage    = fluxes(protIdx);
capUsage    = absUsage./ecModel.ub(protIdx);

absUsage    = absUsage(b);
capUsage    = capUsage(b);

if ~zero
    nonzero     = absUsage>0;
    absUsage    = absUsage(nonzero);
    capUsage    = capUsage(nonzero);
end
end
