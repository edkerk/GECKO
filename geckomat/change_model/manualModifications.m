%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = manualModifications(model)
%
% Benjamin J. Sanchez. Last edited: 2017-10-29
% Ivan Domenzain.      Last edited: 2018-05-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,modifications] = manualModifications(model)

%Read manual data:
fID           = fopen('../../databases/manual_data.txt');
data          = textscan(fID,'%s %s %s %s %f','delimiter','\t');
structure     = data{2};
protGenes     = data{4};
kcats         = data{5}.*3600;
data          = load('../../databases/ProtDatabase.mat');
swissprot     = data.swissprot;
kegg          = data.kegg;
fclose(fID);
modifications{1} = cell(0,1);
modifications{2} = cell(0,1);

%Construct curated complexes:
uniprots = cell(size(kcats));
stoich   = cell(size(kcats));
for i = 1:length(kcats)
    uniprots{i}  = strsplit(structure{i},' + ');
    stoich{i}    = ones(size(uniprots{i}));
    %Separate complex strings in units and amount of each unit:
    for j = 1:length(uniprots{i})
        unit = uniprots{i}{j};
        pos  = strfind(unit,' ');
        if isempty(pos)
            stoich{i}(j)   = 1;
            uniprots{i}{j} = unit;
        else
            stoich{i}(j)   = str2double(unit(1:pos-1));
            uniprots{i}{j} = unit(pos+1:end);
        end
    end
end

for i = 1:length(model.rxns)
    reaction = model.rxnNames{i};
    %Find set of proteins present in rxn:
    S        = model.S;
    subs_pos = find(S(:,i) < 0);
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos);
    prot_set = cell(size(int_pos));
    MW_set   = 0;
    for j = 1:length(int_pos)
        met_name    = model.mets{int_pos(j)};
        prot_set{j} = met_name(6:end);
        MW_set      = MW_set + model.MWs(strcmp(model.enzymes,prot_set{j}));
    end
    %Find intersection with manual curated data:
    for j = 1:length(uniprots)
        int    = intersect(prot_set,uniprots{j});
        if length(int)/max(length(prot_set),length(uniprots{j})) > 0.50 % 50% match
            %Erase previous protein stoich. coeffs from rxn:
            for k = 1:length(prot_set)
                model.S(int_pos(k),i) = 0;
            end
            %If some proteins where not present previously, add them:
            newMets = uniprots{j};
            grRule  = protGenes{j};
            for k = 1:length(uniprots{j})
                if sum(strcmp(model.enzymes,uniprots{j}{k})) == 0
                    model = addProtein(model,uniprots{j}{k},kegg,swissprot);
                end
                newMets{k} = ['prot_' newMets{k}];
            end
            %Add new protein stoich. coeffs to rxn:
            kvalues = kcats(j)./stoich{j};
            rxnID   = model.rxns{i};
            rxnName = model.rxnNames{i};
            model   = addEnzymesToRxn(model,kvalues,rxnID,newMets,{rxnID,rxnName},grRule);
        end
    end
    %Update int_pos:
    S        = model.S;
    subs_pos = find(S(:,i) < 0);
    %Get the proteins that are part of the i-th rxn
    prot_pos = find(~cellfun(@isempty,strfind(model.mets,'prot_')));
    int_pos  = intersect(subs_pos,prot_pos)';
    %%%%%%%%%%%%%%%%%%%%%%%%%%  Individual Changes:  %%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:length(int_pos)
        enzName = model.mets(int_pos(j));
        %%%%%%%%%%%%%%%%%% MANUAL CURATION FOR TOP GROWTH LIMITING ENZYMES:
        [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications);
        if ~isempty(newValue)
            model.S(int_pos(j),i) = newValue;
        end
    end
    disp(['Improving model with curated data: Ready with rxn #' num2str(i)])
end

%%%%%%%%%%%%%%%%%%%%%%%%% Other manual changes: %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove repeated reactions (2017-01-16):
rem_rxn = false(size(model.rxns));
for i = 1:length(model.rxns)-1
    for j = i+1:length(model.rxns)
        if isequal(model.S(:,i),model.S(:,j)) && model.lb(i) == model.lb(j) && ...
                model.ub(i) == model.ub(j)
            rem_rxn(j) = true;
            disp(['Removing repeated rxn: ' model.rxns{i} ' & ' model.rxns{j}])
        end
    end
end
model = removeReactions(model,model.rxns(rem_rxn),true);
% Merge arm reactions to reactions with only one isozyme (2017-01-17):
arm_pos = zeros(size(model.rxns));
p       = 0;
for i = 1:length(model.rxns)
    rxn_id = model.rxns{i};
    if contains(rxn_id,'arm_')
        rxn_code = rxn_id(5:end);
        k        = 0;
        for j = 1:length(model.rxns)
            if ~isempty(strfind(model.rxns{j},[rxn_code 'No']))
                k      = k + 1;
                pos    = j;
                grRule = model.grRules{j};
            end
        end
        if k == 1
            %Condense both reactions in one:
            equations.mets          = model.mets;
            equations.stoichCoeffs  = model.S(:,i) + model.S(:,pos);
            model = changeRxns(model,model.rxns(pos),equations);
            model.grRules{pos} = grRule;
            p          = p + 1;
            arm_pos(p) = i;
            disp(['Merging reactions: ' model.rxns{i} ' & ' model.rxns{pos}])
        end
    end
end
% Remove saved arm reactions:
model = removeReactions(model,model.rxns(arm_pos(1:p)),true);

% Remove unused enzymes after manual curation (2017-01-16):
rem_enz = false(size(model.enzymes));
for i = 1:length(model.enzymes)
    pos_met = strcmp(model.mets,['prot_' model.enzymes{i}]);
    if sum(model.S(pos_met,:)~=0) == 1
        rem_enz(i) = true;
    end
end
rem_enz = model.enzymes(rem_enz);
for i = 1:length(rem_enz)
    model = deleteProtein(model,rem_enz{i});
    disp(['Removing unused protein: ' rem_enz{i}])
end

% Remove incorrect pathways:
model = removeIncorrectPathways(model);

% Map the index of the modified Kcat values to the new model (after rxns removals)
modifications = mapModifiedRxns(modifications,model);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function modified = mapModifiedRxns(modifications,model)
modified = [];
for i=1:length(modifications{1})
    rxnIndex = find(strcmp(model.rxnNames,modifications{2}(i)),1);
    str      = {horzcat(modifications{1}{i},'_',num2str(rxnIndex))};
    modified = [modified; str];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the top growth limiting enzymes that were detected by the
% modifyKcats.m script in a preliminary run.
function [newValue,modifications] = curation_growthLimiting(reaction,enzName,MW_set,modifications)
newValue = [];
reaction = string(reaction);
% chorismate synthase (Q9KXQ4/EC4.2.3.5) - model specified substrate 
% 5-O-(1-Carboxyvinyl)-3-phosphoshikimate is synonymous to the BRENDA
% specified substrate 5-enolpyruvylshikimate 3-phosphate. Changed manually
% to match kcat of wild-type N. crassa enyzme (doi: 10.1111/
% j.1742-4658.2008.06305.x)
if strcmpi('prot_Q9KXQ4',enzName)
    if contains(reaction,'chorismate synthase (')
        newValue         = -(0.87*3600)^-1; % BRENDA: WT N. crassa
        modifications{1} = [modifications{1}; 'Q9KXQ4'];
        modifications{2} = [modifications{2}; reaction];
    end
    
% phosphoribosylformylglycinamidine synthase (Q9RKK5/EC6.3.5.3) - kcat
% automatically suggested was for ammonium as substrate, not glutamine.
% Instead, use specific activity as provided in BRENDA, which was asayed
% with glutamine according to original paper (doi:10.1021/bi00432a017)
elseif strcmpi('prot_Q9RKK5',enzName)
    if contains(reaction,'phosphoribosylformylglycinamidine synthase (')
        newValue         = -(2.15*60*MW_set)^-1; % BRENDA: WT E. coli
        modifications{1} = [modifications{1}; 'Q9RKK5'];
        modifications{2} = [modifications{2}; reaction];
    end
% methylmalonate-semialdehyde dehydrogenase (malonic semialdehyde)
% (Q9L1J1/EC6.3.5.3) - kcat automatically assigned was from archaea,
% instead use value from Bacillus subtilis (doi:10.1074/jbc.M110.213280).
elseif strcmpi('prot_Q9L1J1',enzName)
    if contains(reaction,'methylmalonate-semialdehyde dehydrogenase (malonic semialdehyde) (')
        newValue         = -(2.2*60*MW_set)^-1; % BRENDA: WT B. subtilis
        modifications{1} = [modifications{1}; 'Q9L1J1'];
        modifications{2} = [modifications{2}; reaction];
    end
% phosphoribosyl-ATP pyrophosphatase (Q9EWK0/EC3.6.1.31) - kcat
% automatically assigned was calculated from specific activity in
% Salmonella enterica, but the reported value was measured in cell extract,
% not from purified enzyme. Instead, use specific activity from
% S. cerevisiae (PMID:379004).
elseif strcmpi('prot_Q9EWK0',enzName)
    if contains(reaction,'phosphoribosyl-ATP pyrophosphatase (')
        newValue         = -(332*60*MW_set)^-1; % BRENDA: WT E. coli
        modifications{1} = [modifications{1}; 'Q9EWK0'];
        modifications{2} = [modifications{2}; reaction];
    end
% % glyceraldehyde-3-phosphate dehydrogenase (Q9Z518/EC1.2.1.12) - kcat
% elseif strcmpi('prot_Q9Z518',enzName)
%     if contains(reaction,'glyceraldehyde-3-phosphate dehydrogenase (')
%         newValue         = -(332*60*MW_set)^-1; % BRENDA: WT B. subtilis
%         modifications{1} = [modifications{1}; 'Q9Z518'];
%         modifications{2} = [modifications{2}; reaction];
    end
end

