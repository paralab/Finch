clear all;
clc;
close all;

DIM=4;
% X = sym('x',[1 81],'real');
eleOrder=2;
nrp=eleOrder + 1;
npe=nrp^DIM;
I = eye(DIM);
K1 = sym('A',[1 nrp*nrp],'real');
K2 = sym('B',[1 nrp*nrp],'real');
K3 = sym('C',[1 nrp*nrp],'real');
K4 = sym('D',[1 nrp*nrp],'real');

% M1 = sym('D',[1 nrp*nrp],'real');
% M2 = sym('E',[1 nrp*nrp],'real');
% M3 = sym('F',[1 nrp*nrp],'real');
% 
K1 = reshape(K1,nrp,nrp);
K2 = reshape(K2,nrp,nrp);
K3 = reshape(K3,nrp,nrp);
 K4 = reshape(K4,nrp,nrp);
% 
% M1 = reshape(M1,nrp,nrp);
% M2 = reshape(M2,nrp,nrp);
% M3 = reshape(M3,nrp,nrp);

% K = sym('K', [2 2]);

% IIIAX = kron(I,kron(I,kron(I,A)))*X';
% IIAIX = kron(I,kron(I,kron(A,I)))*X';
% IAIIX = kron(I,kron(A,kron(I,I)))*X';
% AIIIX = kron(A,kron(I,kron(I,I)))*X';

% K = sym('B',[1 4],'real');
mat1 = kron(K1,kron(K2,kron(K3,K4)));
% mat2 = kron(M1,kron(M2,M3));


mat=reshape(mat1,npe*npe,1);
val = ccode(mat);
fid = fopen('4D.txt','w');
fprintf(fid,'%s', val);
fclose(fid);