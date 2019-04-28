function [Ke, Fe] = linearTriangElement(coeff,nodes,elem,e)
%For triangular elements, returns the element stiff matrix, 
%Ke, and the Fe vector for each element depending on the 
%coeff vector for the model equation.
%
% coeff: coefficient vector = [a11,a12,a21,a22,a00,f] for 
%        the model equation.
% nodes: matrix with the coordinates of the nodes.
%  elem: connectivity matrix defining the elements.
%     e: number of the present element

%Compute Ke, Fe for each element
%
v1=nodes(elem(e,1),:);
v2=nodes(elem(e,2),:);
v3=nodes(elem(e,3),:);
Area=0.5*det([1,v1;1,v2;1,v3]);
beta(1)=v2(2)-v3(2);
gamma(1)=v3(1)-v2(1);
% cyclic permutation
beta(2)=v3(2)-v1(2);
gamma(2)=v1(1)-v3(1);
beta(3)=v1(2)-v2(2);
gamma(3)=v2(1)-v1(1);
%
% element stiff matrix
%
Ke=zeros(3);
Fe=zeros(3,1);
if (coeff(1) ~= 0) %a11
    for i=1:3
        for j=1:3
            K11(i,j)=beta(i)*beta(j);
        end
    end
        Ke=Ke+coeff(1)*K11/(4*Area);
end
if (coeff(2) ~= 0) %a12
    for i=1:3
        for j=1:3
            K12(i,j)=beta(i)*gamma(j);
        end
    end
    Ke=Ke+coeff(2)*K12/(4*Area);
end
if (coeff(3) ~= 0) %a21
    for i=1:3
        for j=1:3
            K21(i,j)=gamma(i)*beta(j);
        end
    end
    Ke=Ke+coeff(3)*K21/(4*Area);
end
if (coeff(4) ~= 0) %a22
    for i=1:3
        for j=1:3
            K22(i,j)=gamma(i)*gamma(j);
        end
    end
    Ke=Ke+coeff(4)*K22/(4*Area);
end
if (coeff(5) ~= 0) %a00
    K00=Area/12*[2,1,1;1,2,1;1,1,2];
    Ke=Ke+coeff(5)*K00;
end
if (coeff(6) ~= 0) %f
    Fe=(coeff(6)*Area/3)*[1;1;1];
end

end