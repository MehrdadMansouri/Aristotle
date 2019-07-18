function [PG,NameG] = Patient2SNP(DataDirectory,n0)

%% Importing Patients
Nuc = {'A','C','G','T','?'};
PG0 = readtable([DataDirectory '\AnthracyclineGWAS_data.csv']);
NameG0 = PG0.Properties.VariableNames;
PG0 = table2cell(PG0);
PG0 = PG0(:,8:end);
NameG0 = NameG0(8:end)';
[np,ng0] = size(PG0);
disp('Patients Import Done')

%% SNP Binary-ization
PG = zeros(np,2*ng0); NameG = cell(2*ng0,1); Keep = ones(2*ng0,1);
fa = zeros(5,1);
for ig = 1:ng0
    Pg = PG0(:,ig);
    cPg = strjoin(Pg);
    for ic = 1:4
        fa(ic) = numel(strfind(cPg,Nuc{ic}));
    end
    [~,aG] = max(fa);
    agn = Nuc{aG};
    for ip = 1:np
        pg = strrep(Pg{ip},Nuc{5},agn);
        if isempty(strfind(pg,agn))
            PG(ip,ig*2-1) = 1;
        end
        if not(strcmp(pg,[agn,'_',agn]))
            PG(ip,ig*2) = 1;
        end
    end
    NameG{ig*2-1} = [NameG0{ig},'_',agn,agn];
    NameG{ig*2} = [NameG0{ig},'_',agn];
    s1 = sum(PG(:,ig*2-1));
    s2 = sum(PG(:,ig*2));
    if s1<n0 || s1>np-n0 || aG == 5
        Keep(ig*2-1) = 0;
    end
    if s2<n0 || s2>np-n0 || aG == 5 || s2<s1+2
        Keep(ig*2) = 0;
    end
end
PG = PG(:,find(Keep));
disp('SNP Binaryization Done')
