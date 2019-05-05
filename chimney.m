clearvars
close all

%Define Geometry: node coordinates and elements
nodes=[
    0.1,0;
    0.2,0;
    0.1,0.1;
    0.2,0.1;
    0.2,0.2;
    ];
elem=[4,3,1;
      1,2,4;
      3,4,5;
      ];
numNod= size(nodes,1);
numElem= size(elem,1);
numbering= 1; %=1 shows the numbers for nodes and elements
plotElements(nodes,elem,numbering);

%Define Coefficients vector of the model equation
a11=1;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=0;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
K=zeros(numNod);
F=zeros(numNod,1);
Q=zeros(numNod,1);
for e=1:numElem
    [Ke, Fe] = linearTriangElement(coeff,nodes,elem,e);
    %
    % Assemble the elements
    %
    rows=elem(e,:)';
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end %end for elements
Kini= K; %We save a copy of the initial K and F arrays
Fini= F; %for the post-process step 

%Boundary conditions
fixedNodes= [1,3];
freeNodes= setdiff(1:numNod,fixedNodes);

%------------ Convection BC
convecNodes=[2,4,5];
beta= 20;
Tinf= 30;
[K,Q]= applyConvTriang(convecNodes,beta,Tinf,K,Q,nodes,elem);

%------------ Essential BC
u= zeros(numNod,1); %initialize u vector
u(fixedNodes)=100.0;
Fm= F(freeNodes)-K(freeNodes,fixedNodes)*u(fixedNodes);

%Reduced system
Km=K(freeNodes,freeNodes);
Fm=Fm+Q(freeNodes);

%Compute the solution
format short e
um=Km\Fm;
u(freeNodes)=um;

%Post Process
Q=Kini*u-Fini;
titol='Temp Distribution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale)

%Fancy output
fprintf('%8s%5s%9s%10s%11s\n','Num.Nod','X','Y','T','Q')
fprintf('%5d%10.4f%9.4f%12.4e%12.4e\n',[(1:numNod)',nodes,u,Q]')