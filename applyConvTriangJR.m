function [K,Q]=applyConvTriangJR(indCV,beta,Tinf,K,Q,nodes,elem)
%
%                Only for testing purposes
% 
%         DO      NOT     USE     THIS      FUNCTION
%------------------------------------------------------------------------
% (c) Numerical Factory 2018 Prof. Antonio Susín
%
% Apply a convection BC on a boundary of a triangulated domain. 
%
% Input:
% indCV -> List of nodes on the convection boundary. They can be
%          non-connected but it needs: 
%          i) Each connected component contains at least two nodes 
%          ii) No corners are accepted (two edges of the same triangle). 
%          If an unique triangle is a corner of the CV bondary, we have two call twice 
%          this function (one time with each edge) and sume contributions
%
%          *           *               
%          | \         |               
%          |  \    =   |   +       
%          *---*       *     *---*  
%
% beta -> convection coefficient
% Tinf -> bulk temperature
% K    -> Global stiffness Matrix
% Q    -> Global Secondary variables Vector
% nodes-> all nodes (matrix of size nunNodx2)
% elem -> all elements (matrix of size numElemx3)
%
%------------------------------------------------------------------------
% Output:
% K -> Modified stiffness Matrix with the convection terms Kc on exit
%       K(indCV,indCV)=K(indCV,indCV)+Kc;
% Q -> The same for Q
%       Q(indCV)=Q(indCV)+Qc;
%
%------------------------------------------------------------------------
numElem=size(elem,1); 
numCov=length(indCV);
if numCov==1 
	error('applyConvTriang: Not unic node allow'); 
end

Kaux=(beta/6)*[2 1; 1 2]; 
Faux=0.5*beta*Tinf*[1; 1];
for e=1:numElem
        nodesElem=elem(e,:);
        rows=intersect(indCV,nodesElem);
        numRows=size(rows,2);
        if numRows > 1
            if numRows == 3
                error('Error: no corners allowed!')
            else
                cols = rows;
                h = norm(nodes(rows(1,1),:)-nodes(rows(1,2),:));
                K(rows,cols)=K(rows,cols)+h*Kaux; %stiffness convection matrix
                Q(rows)=Q(rows)+h*Faux;           %convection vector 
            end
        end
end
end
