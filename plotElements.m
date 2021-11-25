function plotElements(nodes, elem, numbering)
%
%(c)Numerical Factory 2021
%
if (numbering==1)
    nE=1;  nN=1;
else 
    nE=0;  nN=0;
end    
[numNodes, dim] = size(nodes);
[numElem, nod4Elem] = size(elem);

if (nod4Elem == 3) %triangle edges
    edges = unique(sort([elem(:,[1,2]);elem(:,[2,3]);elem(:,[3,1])],2),'rows');
elseif (nod4Elem == 4) %quadrilateral edges
    edges = unique(sort([elem(:,[1,2]);elem(:,[2,3]);elem(:,[3,4]);elem(:,[4,1])],2),'rows');
elseif (nod4Elem == 2) %quadrilateral edges
    edges = elem(:,[1,2]);
end    

%Computation of the window size
xmin = min(nodes(:,1));
xmax = max(nodes(:,1));
if (dim==1)
    ymax=2;
    ymin=-2;
else
    ymax = max(nodes(:,2));
    ymin = min(nodes(:,2));
end
lx = xmax-xmin;
ly = ymax-ymin;
xmin = xmin - 0.1*lx;
xmax = xmax + 0.1*lx; 
ymin = ymin - 0.1*ly;
ymax = ymax + 0.1*ly;
lx = xmax-xmin;
ly = ymax-ymin;
figure('units','normalized','outerposition',[0 0 1 1]) % full screen figure
set(groot,'DefaultAxesFontSize',15)
grid on
hold on;
%Line widths and marker/font sizes adapted to the window size and the mesh
minDist1=min(sum(abs(nodes(edges(:,1),:)-nodes(edges(:,2),:)),2));
ratio=minDist1/max(lx,ly);
lwNode=max(0.5,min(3,200*ratio));
if (numElem > 500), lwNode=0.9; end %just to plot big meshes
%lwNodeBd = lwNode;
lwEdge = 2*lwNode;
%lwEdgeBd = 2*lwNode;
if nN
    msNode = 10*lwNode;
%     msNodeBd = 8*lwNode;
else
    msNode = 4*lwNode;
%     msNodeBd = 4*lwNode;
end
fsElem = 8*lwNode;
fsNode = 8*lwNode;
shift=0;

if (dim <= 2)
    if (nod4Elem <= 3 && dim == 1) %1D bar elements
        nodes(:,2)=0;
        shift=0.2;

    elseif(nod4Elem == 2 && dim == 2) % 2D bar elements
        shift=0.03*ly;

    end
        axis equal
        axis([xmin xmax ymin ymax])
        patch('Faces',elem,'Vertices',nodes,'FaceColor','none','Linewidth',lwEdge)
        plot(nodes(:,1),nodes(:,2),'ko',...
            'MarkerFaceColor',[0.99,1,0.8],'Markersize',msNode,...
            'LineWidth',lwNode);
        %Labels for numbering nodes and elements, trying to get a "nice" figure
        if nE
            for e=1:numElem
                xy=sum(nodes(elem(e,:),:),1)/nod4Elem;
                labelNode = ['{\boldmath$',num2str(e),'$}'];
                if (nod4Elem == 2 && dim == 2) %2D bar elements
                    vect=nodes(elem(e,2),:)-nodes(elem(e,1),:);
                    vect=vect/norm(vect);
                    xy(1)=xy(1)-shift*vect(2);
                    xy(2)=xy(2)+shift*vect(1);
                else    
                xy(2)=xy(2)+shift;
                end
                text(xy(1),xy(2),labelNode,'interpreter','LaTeX',...
                    'FontSize',fsElem,'HorizontalAlignment','center')
            end
        end
        if nN
            for n=1:numNodes
                xy=nodes(n,:);
                labelNode = ['{\boldmath$',num2str(n),'$}'];
                text(xy(1),xy(2),labelNode,'interpreter','LaTeX',...
                    'FontSize',fsNode,'HorizontalAlignment','center')
            end
        end

        h=title(['nNodes = ',num2str(numNodes),',  nElem = ',num2str(numElem)]);
        set(h,'FontSize',25,'FontWeight','bold') 

        hold off;
else %3 dimensions
    zmin = min(nodes(:,3));
    zmax = max(nodes(:,3));

    if(nod4Elem == 2 && dim == 3) %3D bar elements
%         plot3(nodes(:,1),nodes(:,2),nodes(:,3),'ro','LineWidth',1.2,...
%                     'MarkerEdgeColor','k',...
%                     'MarkerFaceColor','g',...
%                     'MarkerSize',5);           
%         %axis([xmin-shift, xmax+shift, ymin-shift, ymax+shift, zmin-shift, zmax+shift])
%         view(120,20);
%         hold on;
%         %shiftPoints=0.05*max(max(abs(nodes)));
%       if (numbering==1)
%         for n=1:numNodes
%             labelNode = ['{\boldmath$n_{',num2str(n),'}$}'];
%             text(nodes(n,1)+shift,nodes(n,2)+shift,nodes(n,3)+shift,...
%                 labelNode,'interpreter','LaTeX','fontSize',12)
%         end
%       end
%         %Plot elements
%         shiftEdges=shift;
%         for e=1:numElem
%             labelEdge = ['{\boldmath$E_{',num2str(e),'}$}'];
%             plot3(nodes(elem(e,:),1),nodes(elem(e,:),2),nodes(elem(e,:),3),'b-');
%             p=[sum(nodes(elem(e,:),1)),sum(nodes(elem(e,:),2)),sum(nodes(elem(e,:),3))]/2;
%       if (numbering==1)
%             text(p(1,1)+shiftEdges,p(1,2)+shiftEdges,p(1,3)+shiftEdges,...
%                 labelEdge,'interpreter','LaTeX','fontSize',12);
%         end    
%         end    
% 
%        axis equal vis3d
%         hold off; 
% 
        shift=0.03*ly;
        lz=zmax-zmin;
        dl = 0.03*min([lx,ly,lz]);
        axis equal
        axis([xmin xmax ymin ymax zmin zmax])
        %patch('Faces',elem,'Vertices',nodes,'FaceColor','none','Linewidth',0.8*lwEdge)
        for e=1:numElem
            pa = nodes(elem(e,end),:);
            pb = nodes(elem(e,1),:);
            vect = pb-pa;
            vect = vect/norm(vect);
            qa = pa + dl*vect;
            qb = pb - dl*vect;
            plot3([qa(1),qb(1)],[qa(2),qb(2)],[qa(3),qb(3)],'k-',...
            'LineWidth',lwNode);
        end
        plot3(nodes(:,1),nodes(:,2),nodes(:,3),'ko',...
            'MarkerFaceColor',[0.99,1,0.8],'Markersize',msNode);
        az = 120;
        el = 20;
        view(az,el);
        camera = zeros(3,1);
        camera(1) = cos(pi*el/180)*sin(pi*az/180);
        camera(2) = -cos(pi*el/180)*cos(pi*az/180);
        camera(3) = sin(pi*el/180);
        %Labels for numbering nodes and elements, trying to get a "nice" figure
        if nE
            for e=1:numElem
                xy=sum(nodes(elem(e,:),:),1)/nod4Elem;
                labelNode = ['{\boldmath$',num2str(e),'$}'];
                vect=nodes(elem(e,end),:)-nodes(elem(e,1),:);
                vect=vect/norm(vect);
                vect=cross(vect,camera);
                xy=xy-shift*vect;
                text(xy(1),xy(2),xy(3),labelNode,'interpreter','LaTeX',...
                    'FontSize',fsElem,'HorizontalAlignment','center')
            end
        end
        if nN
            for n=1:numNodes
                xy=nodes(n,:);
                labelNode = ['{\boldmath$',num2str(n),'$}'];
                text(xy(1),xy(2),xy(3),labelNode,'interpreter','LaTeX',...
                    'FontSize',fsNode,'HorizontalAlignment','center')
            end
        end
        
        axis equal vis3d


        h=title(['nNodes = ',num2str(numNodes),',  nElem = ',num2str(numElem)]);
        set(h,'FontSize',25,'FontWeight','bold') 

        hold off;



    elseif (nod4Elem == 4 && dim == 3) % case Tetrahedra
        TRI=[];
        for e=1:numElem   
            TRI=[TRI; elem(e,[1,2,3])]; %face 1
            TRI=[TRI; elem(e,[1,3,4])];
            TRI=[TRI; elem(e,[1,2,4])];
            TRI=[TRI; elem(e,[2,3,4])];
        end
        trimesh(TRI,nodes(:,1),nodes(:,2),nodes(:,3),'EdgeColor', [0.1, 0.8, 0.8]);
        view(120,20);
        axis equal vis3d


        hold off;
    end
end



