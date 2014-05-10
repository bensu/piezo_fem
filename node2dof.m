function v = node2dof(nodeList, nDof)

nNod = length(nodeList);

v = zeros(nNod,nDof);

for iNod = 1:nNod
    v(iNod,:) =  ((nodeList(iNod) - 1)*nDof + 1) : nodeList(iNod)*nDof;
end

v = reshape(v',1,[]);
    