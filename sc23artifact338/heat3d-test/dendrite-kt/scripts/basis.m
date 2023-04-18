clear all;
clc;
close all;
I = eye(2);
%  X = sym('x',[1 81],'real');
A = sym('A',[1 4],'real');
% 
A = reshape(A,2,2);


K1 = sym('K1', [2 2]);
K2 = sym('K2', [2 2]);
K3 = sym('K2', [2 2]);

val1 = kron(K1,kron(K2,K3));
val2 = kron(K1,kron(I,I))*kron(I,kron(K2,I))*kron(I,kron(I,K3));
% 
% IIIAX = kron(I,kron(I,kron(I,A)))*X';
% IIAIX = kron(I,kron(I,kron(A,I)))*X';
% IAIIX = kron(I,kron(A,kron(I,I)))*X';
% AIIIX = kron(A,kron(I,kron(I,I)))*X';

