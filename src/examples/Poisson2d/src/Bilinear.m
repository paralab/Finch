%This file was generated by Finch.

%{
Bilinear term
%}

dof = size(grid_data.allnodes,2);
ne  = mesh_data.nel;
Np = refel.Np;
I = zeros(ne * Np*Np, 1);
J = zeros(ne * Np*Np, 1);
val = zeros(ne * Np*Np, 1);

% loop over elements
for e=1:ne
    idx = grid_data.loc2glb(:,e)';
    pts = grid_data.allnodes(:,idx);
    [detJ, Jac]  = Utils.geometric_factors(refel, pts);
    
    ind1 = repmat(idx,Np,1);
    ind2 = reshape(repmat(idx',Np,1),Np*Np,1);
    st = (e-1)*Np*Np+1;
    en = e*Np*Np;
    I(st:en) = ind1;
    J(st:en) = ind2;

elMat = ([diag(Jac.rx) diag(Jac.sx)] * [refel.Qr; refel.Qs])' * diag(-refel.wg .* detJ) * [diag(Jac.rx) diag(Jac.sx)] * [refel.Qr; refel.Qs] + ([diag(Jac.ry) diag(Jac.sy)] * [refel.Qr; refel.Qs])' * diag(-refel.wg .* detJ) * [diag(Jac.ry) diag(Jac.sy)] * [refel.Qr; refel.Qs];


    val(st:en) = elMat(:);
end
LHS = sparse(I,J,val,dof,dof);

for i=1:length(grid_data.bdry)
    LHS(grid_data.bdry{i},:) = 0;
    LHS((size(LHS,1)+1)*(grid_data.bdry{i}-1)+1) = 1;
end