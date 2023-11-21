function [Q]=applyConstantNaturalBC(nodes,elem,indBC,q0,Q)
%------------------------------------------------------------------------
% (c) Numerical Factory 2020
%
% Apply a constant natural BC on a boundary of a triangulated domain. 
% 
% INPUT
% --------
% indBC: horizontal vector with the boundary list of nodes
% q0: constant value for this natural BC
% Q: global vector of Q's
% nodes, elem: mesh tables
% 
% OUTPUT
% --------
% Q: global vector with modified values according the constant natural BC.
%------------------------------------------------------------------------
[numElem, ndim]=size(elem); 
numBC=length(indBC);
if numBC==1 
	error('applyNaturalBC: Not unic node allow'); 
end
if (ndim==3) % Triangle elements
    for k=1:numElem
      aux=[0,0,0]; %initial values (no node is found)
      for inod=1:3 %loop for the three nodes
          r=find(indBC == elem(k,inod)); %find if one node of element k is in the convection node list 
          if(~isempty(r)) 
              aux(inod)=1; %put 1 on the local position found
          end
      end
        %      
        % Now the aux vector is one of the following possibilities
        % aux=[0,0,0]  -> this element does not contain any BC node
        % aux=[1,0,0], [0,1,0], [0,0,1] -> contains only one BC node 
        % aux=[1,1,0], [1,0,1], [0,1,1] -> contains two BC nodes
        % To classify which is the edge in local numbering, we consider aux as
        % a binary number expressed as power of 2: aux(1)*2^0+aux(2)*2^1+aux(3)*2^2
        % only 3 (first edge), 5 (third edge) or 6 (second edge)
        % are significant.
      number = aux(1)+2*aux(2)+4*aux(3); 
      switch (number) %identify the appropriate edge
          case 3
              ij=[1,2];
          case 5
              ij=[3,1];
          case 6
              ij=[2,3];
          case 7
              error('applyConvTriang: Corners not allowed !!!!\n');           
          otherwise, ij=[0,0];
      end
      if ( ij(1) > 0) %it's an existing edge 
       n1=elem(k,ij(1)); 
       n2=elem(k,ij(2));
       h=norm(nodes(n1,:)-nodes(n2,:));
       fico=[n1,n2];
       Q(fico)=Q(fico)+0.5*q0*h*[1; 1];           
      end
    end
else
    for k=1:numElem
      aux=[0,0,0,0];
      for j=1:4
        r=find(indBC==elem(k,j)); %find if j node of element k is in the convection nodes list 
          if(~isempty(r)), aux(j)=1; end
      end
      %See the Triang version for the explanation of the binary number    
      number=aux*[1;2;4;8];
      % Here we codify both edges and corners on the boundaries
      switch (number) 
        case  3, ij=[1,2,0,0];  % edge 1:   aux=[1,1,0,0];
        case  6, ij=[2,3,0,0];  % edge 2:   aux=[0,1,1,0];
        case  7, ij=[1,2,2,3];  % corner 2: aux=[1,1,1,0];
        case  9, ij=[4,1,0,0];  % edge 4
        case 11, ij=[4,1,1,2];  % corner 1 
        case 12, ij=[3,4,0,0];  % edge 3
        case 13, ij=[3,4,4,1];  % corner 4
        case 14, ij=[2,3,3,4];  % corner 3
        case 15, error('applyConstantNaturalBCQuad: It can''t manage two corners'); % two corners 
        otherwise, ij=[0,0,0,0];
      end
      for j=1:2:3
         if ~ij(j)==0
            n1=elem(k,ij(j)); 
            n2=elem(k,ij(j+1));
            h=norm(nodes(n1,:)-nodes(n2,:)); 
            fico=[n1, n2];   
            Q(fico)=Q(fico)+0.5*q0*h*[1; 1];           
         end
       end  
    end
  end
end