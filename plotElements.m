function plotElements(nodes, elem, numbering)
%
%(c)Numerical Factory 2016
%

numNod=size(nodes,1);
numElem=size(elem,1);
xmax=max(nodes(:,1));
xmin=min(nodes(:,1));
ymax=max(nodes(:,2));
ymin=min(nodes(:,2));

shift=2*max(0.02*max(abs([xmax, xmin,ymin,ymax])),0.01);
if (size(elem,2) == 3) %Case Triangles
    A = sparse([elem(:,1); elem(:,1); elem(:,2); elem(:,2); elem(:,3); elem(:,3)], ...
    [elem(:,2); elem(:,3); elem(:,1); elem(:,3); elem(:,1); elem(:,2)],1.0);
    A = min(A, 1);
    gplot(A,nodes)
    hold on;
    h = title('Triangle Elements Visualizer');
    set(h,'FontSize',11,'FontWeight','bold') 
    plot(nodes(:,1),nodes(:,2),'ro','LineWidth',1.5,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',2.5);         
    axis([xmin-2*shift, xmax+2*shift, ymin-2*shift, ymax+2*shift]);
    axis equal;
    if (numbering==1)
        for j=1:numNod   
            text(nodes(j,1),nodes(j,2)-shift,['n' num2str(j,'%2d')]);
        end
        for e=1:numElem   
            xx=sum(nodes(elem(e,:),:),1)/3;
            text(xx(1),xx(2),['E' num2str(e,'%2d')]);
        end
    end
    hold off;
else
   A = sparse(...
        [elem(:,1); elem(:,1); elem(:,2); elem(:,2); elem(:,3); elem(:,3); elem(:,4); elem(:,4)], ...
        [elem(:,4); elem(:,2); elem(:,1); elem(:,3); elem(:,2); elem(:,4); elem(:,3); elem(:,1)], ...
        1.0);
    A = min(A, 1);
    gplot(A,nodes)
    hold on;
    h = title('Quadrilateral Elements Visualizer');
    set(h,'FontSize',11,'FontWeight','bold') 
    plot(nodes(:,1),nodes(:,2),'ro','LineWidth',1.5,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',2.5);         
    axis([xmin-2*shift, xmax+2*shift, ymin-2*shift, ymax+2*shift]);
    axis equal;

    if (numbering==1)
        for j=1:numNod   
            text(nodes(j,1),nodes(j,2)-shift,['n' num2str(j,'%2d')]);
        end
        for e=1:numElem   
            xx=sum(nodes(elem(e,:),:),1)/4;
            text(xx(1),xx(2),['E' num2str(e,'%2d')]);
        end
    end
    hold off;
end



