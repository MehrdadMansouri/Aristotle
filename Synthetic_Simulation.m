%% Setting
clc; clear variables; close all
dy  = 3;    	% Number of Outcome Subtypes
zy  = 0.1;      % Outcome Noise
ny1 = 150;      % Number of Positive Outcome Samples
ny  = 750;      % Number of Samples
nx  = 100000;	% Number of Features
dx  = 2500;   	% Number of Features Subtypes
dp  = 100;   	% Number of Pathways 
zx  = 0.05;     % Feature Noise
zp  = 0.05;     % Uniformity of Feature Selection for Pathway
%% Feature & Outcome Subgroups' Index  
Ix = Choper(nx,dx); Ix = [Ix;Ix(end)+round(nx/dx)];
Iy = Choper(ny1,dy); Iy = [Iy;Iy(end)+round(ny1/dy)];
%% Outcome Values
Y = [ones(ny1,1);zeros(ny-ny1,1)];
%% Feature Values
X = zeros(ny,nx);
for i = 1:dy
    X(Iy(i):Iy(i+1)-1,Ix(i):Ix(i+1)-1) = 1;
end
for i = dy+1:dx
    X(Picker(ny,ny1/dy),Ix(i):Ix(i+1)-1) = 1;
end
%% Add Noise
Y = Polluter(Y,zy);
X = Polluter(X,zx);
%% Permutate Rows
Perm_Order = randperm(ny);
XX = X(Perm_Order,:);
YY = Y(Perm_Order);
%% Reverse Order Option
[~,Re_Order] = sort(Perm_Order); % Y = YY(Re_Order)
%% Break into Pathways
% Barabasi–Albert
Pathway = zeros(Ix(end)-1,1);
Ip = Choper(nx,dp); Ip = [Ip;Ix(end)];
Unpicked = ones(Ix(end)-1,1);
for i = 1:dp
    P0 = zp*ones(Ix(end)-1,1);
    for j = 1:(Ip(i+1)-Ip(i))
        Px = Unpicked.*P0;
        r = mnrnd(1,Px/sum(Px)); r = find(r);
        Unpicked(r) = 0;
        x0 = find(r>=Ix,1,'last');
        P0(Ix(x0):Ix(x0+1)-1) = P0(Ix(x0):Ix(x0+1)-1) + 1;
        Pathway(r) = i;
    end
end
%% Display
csvwrite('Features9.txt',X)
csvwrite('Outcome9.txt',Y)
csvwrite('Pathway9.txt',Pathway)
% image([X 3*Y],'CDataMapping','scaled')
% A = Pathway(1:Ix(2));
% B = unique(A);
% histc(A,B)
% size(unique(A))
%% Starting Index of b Partitions of total length ~X
function Index = Choper(X,b)
Index = ones(b,1);
for i = 2:b
    Index(i) = Index(i-1) + binornd(X-1,1/b)+1;
end
end
%% Starting Index of Random samples of size ~d up to dmax
function Index = Picker(dmax,d)
Size = binornd(dmax,d/dmax);
Index = randperm(dmax,Size);
end
%% Modified X by floping p portion of it
function Xp = Polluter(X,p)
Changes = randperm(numel(X),round(numel(X)*p));
Xp = X; Xp(Changes) = 1 - Xp(Changes);
end
