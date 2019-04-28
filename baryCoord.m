function [alphas, isInside]=baryCoord(V, point)
%BARYCOORD
%Compute the triangular barycentric coordinates of the a point
%INPUT
%      V (array 3x2): array holding the coordinates of the vertexs of
%                     the triangle V=[x1,y1;x2,y2;x3,y3].
%  point (array 1x2): array holding the coordinates of the point.
%
%OUTPUT
% alphas (array 1x3): array holding the barycentric coordinates of the 
%                     point
% isInside (integer): if isInside=1, then the point belongs to the
%                     triangle, if isInside=0, the point lies outside.
%                

isInside=1; %By default, the point belongs the triangle

%Compute the shape functions
A=[ones(3,1),V];
b=[1;0;0];
c1=A\b;
b=[0;1;0];
c2=A\b;
b=[0;0;1];
c3=A\b;

Psi1=@(x,y) c1(1)+c1(2)*x+c1(3)*y;
Psi2=@(x,y) c2(1)+c2(2)*x+c2(3)*y;
Psi3=@(x,y) c3(1)+c3(2)*x+c3(3)*y;

alpha1=Psi1(point(1,1),point(1,2));
alpha2=Psi2(point(1,1),point(1,2));
alpha3=Psi3(point(1,1),point(1,2));

alphas=[alpha1,alpha2,alpha3];

if (min(alphas) < 0) %at least one of the alphas is negative
    isInside=0;
end

end