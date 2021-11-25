clearvars
close all

fileName='solution.xlsx';

eval('CircleHolemesh01')

numNodes= size(nodes,1);
numElem= size(elem,1);

numbering= 0; %= 1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering)

%Select Boundary points
indT= find(nodes(:,2) > 0.99); %indices of the nodes at the top boundary
indB= find(nodes(:,2) < -0.99);%indices of the nodes at the bottom boundary
indC= find(sqrt(nodes(:,1).^2 + nodes(:,2).^2) < 0.501); %id.on the circle

hold on
plot(nodes(indT,1),nodes(indT,2),'ok','lineWidth',1,'markerFaceColor',...
    'red','markerSize',5)
plot(nodes(indB,1),nodes(indB,2),'ok','lineWidth',1,'markerFaceColor',...
    'blue','markerSize',5)
plot(nodes(indC,1),nodes(indC,2),'ok','lineWidth',1,'markerFaceColor',...
    'green','markerSize',5)
hold off

%Define the coefficients vector of the model equation
a11=1;
a12=0;
a21=a12;
a22=a11;
a00=0;
f=0;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
K=zeros(numNodes);    %global stiff matrix
F=zeros(numNodes,1);  %global internal forces vector
Q=zeros(numNodes,1);  %global secondary variables vector

for e = 1:numElem
    [Ke, Fe] = linearTriangElement(coeff,nodes,elem,e);
    rows= [elem(e,1); elem(e,2); elem(e,3)];
    cols= rows;
    K(rows,cols)= K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)= F(rows) + Fe;
    end
end %end for elements
Kini=K; %We save a copy of the initial K and F arrays
Fini=F; %for the post-process step 

%Booundary Conditions
fixedNodes= [indB', indC'];                %fixed Nodes (global numbering)
freeNodes= setdiff(1:numNodes,fixedNodes); %free Nodes (global numbering)

%------------- Convetion BC
beta=2.0;
Tinf=-5.0;
indCV=indT;
%[K,Q]=applyConvTriangJR(indCV,beta,Tinf,K,Q,nodes,elem); %<--DO NOT USE IT!
[K,Q]=applyConvTriang(indCV,beta,Tinf,K,Q,nodes,elem);

% Essential B.C.
u=zeros(numNodes,1);
u(indB)= 50.0;  %Temperature at the bottom boundary
u(indC)= 15.0;  %Temperature at the circle
Fm = F(freeNodes) - K(freeNodes,fixedNodes)*u(fixedNodes);

%Reduced system
Km = K(freeNodes,freeNodes);
Fm = Fm + Q(freeNodes);

%Compute the solution
um = Km\Fm;
u(freeNodes)= um;

%PostProcess: Compute secondary variables and plot results
Q = Kini*u - Fini;
titol='Temperature Distribution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

%Fancy output
tableSol=[(1:numNodes)',nodes,u,Q];
fprintf('%8s%9s%15s%15s%14s\n','Num.Nod','X','Y','T','Q')
fprintf('%5d%18.7e%15.7e%15.7e%15.7e\n',tableSol')

%write an Excel with the solutions
% format long e
% ts=table(int16(tableSol(:,1)),tableSol(:,2),tableSol(:,3),tableSol(:,4),...
%     tableSol(:,5),'variableNames',{'NumNod','X','Y','T','Q'});
% writetable(ts,fileName);

%Exercise 1:
%Compute the temperature for the point p=[0.5, 0.2].
p= [0.5, 0.2];

for e=1:numElem
    vertexs= nodes(elem(e,:),:);
    [alphas,isInside] = baryCoord(vertexs,p);
    if (isInside > 0)
        pElem = e;
        numNodElem= elem(e,:);
        break;
    end
end

tempP = alphas*u(numNodElem);
fprintf('\nExercise 1:\n')
fprintf('\nCompute the temperature for the point p=[0.5, 0.2].\n')
fprintf('\nSol.:\n\n')
fprintf('Point P = (%.1f,%.1f) belongs to element number: %d\n',p,pElem)
fprintf('Number of nodes of elem %d: %d, %d, %d\n',pElem,numNodElem)
fprintf('Interpolated temperature at point P: %.4e\n',tempP) 
