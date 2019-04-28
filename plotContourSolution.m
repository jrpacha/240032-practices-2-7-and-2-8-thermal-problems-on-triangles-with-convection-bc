function plotContourSolution(nodes,elem,valueToShow,titol,colorScale)
%
% plot color contour values for triangles and quadrilateral meshes
% (c)Numerical Factory
%
 numElem = size(elem,1);     
 numNodesElem = size(elem,2);
for ielem=1:numElem  
 for i=1:numNodesElem
     nd(i)=elem(ielem,i);         % extract nodes for i-th element
     X(i,ielem)=nodes(nd(i),1);    % extract x value of the node
     Y(i,ielem)=nodes(nd(i),2);    % extract y value of the node
 end   
 plotValues(:,ielem) = valueToShow(nd') ; %extract the value for each node 
end
figure()
colormap(colorScale)
 plot(X,Y,'k')
 fill(X,Y,plotValues)
 title(titol) ;         
 axis equal;
 colorbar