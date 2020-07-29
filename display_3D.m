function display_3D(rho,study)
% Decription: visualize the 3D rho input matrix
% ---- Input Arguments ----
% rho: 3d matrix
% study: 
% a) 'unitcell'; Visualization of the unit cells
% b) 'domain'; Visualization of the initial domain
% c) 'topology'; Visualization of the HTOP
%-----------------------------------------------------------------%
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;  % Element size
face = [1 2 3 4; 2 6 7 3; 4 3 7 8; 1 5 8 4; 1 2 6 5; 5 6 7 8];
% Display options
if strcmpi(study,'unitcell')
    Face = reshape(kron([1,.05],ones(ceil(nelz/2),1)),[nelz,1]);
    Edge=reshape(repmat({'k','none'},[ceil(nelz/2),1]),[nelz,1]);
    Face_c1=repmat(.2,[1,3]); 
    Face_c0=[.6974, .4788,.5856];
    Edalpha=ones(nelz,1);
elseif strcmpi(study,'domain')
    Face(logical(rho(:)))=.05; 
    Edge={'k','none'}; 
    Face_c1=repmat(.73,[1,3]);
    Edalpha(logical(rho(:)))=.05;
end
p=0;
for k = 1:nelz
    p=p+1;
    z = (k-1)*hz;
    for i = 1:nelx
        x = (i-1)*hx;
        for j = 1:nely
            y = nely*hy - (j-1)*hy;
            switch study
                case {'unitcell','domain'}  
                    if strcmpi(study,'unitcell')
                            Fc=Face(k); Ed=Edge{k};Eda=Edalpha(k);
                        else
                            Fc=Face(p); Ed='k';Eda=Edalpha(p);
                    end
                    if (rho(j,i,k) == 1)  % Material 1     
                       vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                        vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                        patch('Faces',face,'Vertices',vert,'FaceColor',Face_c1,'FaceAlpha',Fc,'EdgeColor',Ed,'EdgeAlpha',Eda);
                        hold on;
                    elseif(rho(j,i,k)==0)  % Material 2
                        vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                        vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                        patch('Faces',face,'Vertices',vert,'FaceColor',... % [0.6477    0.4509    0.5470]);
                            Face_c0,'FaceAlpha',Fc,'EdgeColor',Ed,'EdgeAlpha',Eda);%char(Edge{Face(p)})
                        hold on;
                    end
                case 'topology'
                    minV = 0.5; % User-defined display density threshold
                    a=-[.9948,.5576,.7712]; b=[1.1948,.7576,.9712];
                    if (rho(j,i,k) >= minV)
                        vert = [x y z; x y-hx z; x+hx y-hx z; x+hx y z; x y z+hx;x y-hx z+hx; x+hx y-hx z+hx;x+hx y z+hx];
                        vert(:,[2 3]) = vert(:,[3 2]); vert(:,2,:) = -vert(:,2,:);
                        patch('Faces',face,'Vertices',vert,'FaceColor',...
                            [a(1)*rho(j,i,k)+b(1),a(2)*rho(j,i,k)+b(2),a(3)*rho(j,i,k)+b(3)]);
                        hold on;
                    end
            end
        end
    end
end
axis equal; axis tight; axis off; box on; view([30,30]); pause(1e-6);
end