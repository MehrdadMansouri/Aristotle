clc; clear variables; close all; warning off
%% Export Patient-SNP Matrix of all Pathways

%% Setting
n0 = 6;
DataDirectory = 'Coclustering\Data';        disp('Initiated ...')
%% Pre-Processing
G2S = Gene2SNP(DataDirectory);              disp('Gene to SNPs Done')
P2S = Pathway2SNP(DataDirectory,G2S);       disp('Pathway to SNPs Done')
% [PG,NameG] = Patient2SNP(DataDirectory,n0); disp('Patients SNPs Done')
PG = load([DataDirectory '\PG1']);

