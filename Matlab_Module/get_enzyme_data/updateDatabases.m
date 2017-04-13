%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% updateDatabases
% Updates all databases for protein matching (KEGG and Swiss-Prot).
%
% Note: Before using this script, one should manually download:
%       *Swissprot: Download a tab delimited file with the following format:
%                   Entry - Protein names - Gene names - EC number - Sequence
%                   http://www.uniprot.org/uniprot/?query=organism:%22yeast%22
%                   OBS: filter with the Swiss-Prot option
% 
% Benjam�n S�nchez. Last edited: 2017-04-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateDatabases

%Retrieve Swissprot Data (Entry - Protein names - Gene names - EC number - Sequence):
cd ../../Databases
fileID_uni        = fopen('uniprot-organism%3Ayeast.tab');
swissprot         = textscan(fileID_uni,'%s %s %s %s %s %s','delimiter','\t');
swissprot         = [swissprot{1} swissprot{2} swissprot{3} swissprot{4} swissprot{5}];
swissprot(1,:)    = [];
fclose(fileID_uni);
cd ../Matlab_Module/get_enzyme_data
for i = 1:length(swissprot)
    %Leave protein name as lower case, remove ';' from ECs & calculate MW:
    prot_name      = lower(swissprot{i,2});
    uni            = swissprot{i,1};
    sequence       = swissprot{i,5};
    MW             = calculateMW(sequence);
    swissprot{i,2} = prot_name;
    swissprot{i,4} = strrep(swissprot{i,4},';','');
    swissprot{i,5} = MW;
    swissprot{i,6} = sequence;
    disp(['Updating Swiss-Prot database: Ready with protein ' uni])
end

%Retrieve KEGG info (gen - uniprot - EC code - name - MW - Pathway):
cd ../../Databases/KEGG
file_names      = dir();
file_names(1:2) = [];
kegg            = cell(100000,7);
n               = 0;
for i = 1:length(file_names)
    file_name = file_names(i).name;
    %1st column: Gene name
    gene_name = file_name(1:end-4);
    %Retrieve all data as a cell with all rows:
    fID  = fopen(file_name);
    text = textscan(fID,'%s','delimiter','\t');
    fclose(fID);
    text = text{1};
    cd ../../Matlab_Module/get_enzyme_data
    
    uni      = '';
    sequence = '';
    MW       = 0;
    pathway  = '';
    for j = 1:length(text)
        line = text{j};
        %2nd column: uniprot number
        if ~isempty(strfind(line,'UniProt:'))
            uni = line(10:end);
            
        %3rd & 4th column: protein name and EC number
        elseif ~isempty(strfind(line,'DEFINITION'))
            pos_EC    = strfind(line,'EC:');
            if isempty(pos_EC)
                prot_name = lower(line(13:end));
                EC_names  = '';
            else
                prot_name = lower(line(13:pos_EC-3));
                EC_names  = line(pos_EC+3:end-1);
            end
            
        %5th column and 7th column: MW & sequence
        elseif ~isempty(strfind(line,'AASEQ'))
            end_seq  = false;
            for k = j+1:length(text)
                if ~isempty(strfind(text{k},'NTSEQ'))
                    end_seq = true;
                elseif ~end_seq
                    sequence = [sequence text{k}];
                end
            end
            MW = calculateMW(sequence);
        
        %6th column: Pathway
        elseif ~isempty(strfind(line,'PATHWAY'))
            start    = strfind(line,'sce');
            pathway  = line(start(1):end);
            end_path = false;
            for k = j+1:length(text)
                nospace = strrep(text{k},'sce01100  Metabolic pathways','');
                nospace = strrep(nospace,' ','');
                if length(nospace) > 10
                    if strcmp(nospace(1:3),'sce') && ~end_path
                        start    = strfind(text{k},'sce');
                        pathway  = [pathway ' ' text{k}(start(1):end)];
                    else
                        end_path = true;
                    end
                end
            end
        end
    end
    %Create one aditional association in kegg:
    n         = n+1;
    kegg{n,1} = uni;
    kegg{n,2} = prot_name;
    kegg{n,3} = gene_name;
    kegg{n,4} = EC_names;
    kegg{n,5} = MW;
    kegg{n,6} = pathway;
    kegg{n,7} = sequence;
    cd ../../Databases/KEGG
    disp(['Updating KEGG database: Ready with gene ' gene_name])
end
kegg(n+1:end,:)         = [];

%Save all databases as .mat files:
cd ..
save('ProtDatabase.mat','kegg','swissprot');
cd ../Matlab_Module/get_enzyme_data

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%