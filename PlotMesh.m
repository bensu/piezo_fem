function PlotMesh(coordinates,nodes,bc,force)
%--------------------------------------------------------------------------
% Code written by : Siva Srinivas Kolukula                                |
%                   Senior Research Fellow                                |
%                   Structural Mechanics Laboratory                       |
%                   Indira Gandhi Center for Atomic Research              |
%                   India                                                 |
% E-mail : allwayzitzme@gmail.com                                         |
%          http://sites.google.com/site/kolukulasivasrinivas/             |
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Purpose:
%         To plot the Finite Element Method Mesh
% Synopsis :
%           PlotMesh(coordinates,nodes)
% Variable Description:
%           coordinates - The nodal coordinates of the mesh
%           -----> coordinates = [node X Y]
%           nodes - The nodal connectivity of the elements
%           -----> nodes = [node1 node2......]
%--------------------------------------------------------------------------
nodes = nodes(:,[1 2 4 3]);
bc = ~bc;
nel = size(nodes,1) ;                  % number of elements
nnode = size(coordinates,1) ;          % total number of nodes in system
nnel = size(nodes,2);                % number of nodes per element

% Initialization of the required matrices
X = zeros(nnel,nel);
Y = zeros(nnel,nel);
Z = zeros(nnel,nel);

for iel = 1:nel
    for i = 1:nnel
        nd(i) = nodes(iel,i);         % extract connected node for (iel)-th element
        X(i,iel) = coordinates(nd(i),1);    % extract x value of the node
        Y(i,iel) = coordinates(nd(i),2);    % extract y value of the node
        Z(i,iel) = coordinates(nd(i),3);
    end
end

% Plotting the FEM mesh, diaplay Node numbers and Element numbers

f1 = figure ;
set(f1,'name','Mesh','numbertitle','off') ;

patch(X,Y,Z,'w'); 
hold on
% fill(X,Y,'w')

quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3), ...
                        force(:,1),force(:,2),force(:,3));
quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3), ...
                        bc(1:nnode,1),bc(1:nnode,2),bc(1:nnode,2));
title('Finite Element Mesh') ;
axis off ;

k = nodes(:,1:end);
nd = k' ;
for i = 1:nel
    text(X(:,i),Y(:,i),Z(:,i),int2str(nd(:,i)),'fontsize',8,'color','k');
    text(sum(X(:,i))/4,sum(Y(:,i))/4,sum(Z(:,i))/4,int2str(i),'fontsize',10,'color','r') ;
end
hold off