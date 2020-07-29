function [voxel,Density] = GenerateVoxel(n,VoxelType,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n: the number of voxel along each axis
% VoxelType: either 'lattice' or 'box'
% p: r or alpha/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
size = 1/n;               % initial size of voxels
voxel = zeros(n,n,n);      % initial grid with zeros
% generate a list of centers of voxel
voxel_c = zeros(n^3,6);
ii = 0;                    % ii count the number of all voxels
for i = 1:n               % i for z axis
    for j = 1:n           % j for y axis
        for k = 1:n       % k for x axis
            ii = ii + 1;
            voxel_c(ii,1:3) = [k,j,i];  % save index along x,y,z axis
            % save coordinate along x,y,z axis
            voxel_c(ii,4:6) = [(k-0.5)*size,(j-0.5)*size,(i-0.5)*size];
        end
    end
end
%% Get the voxel close the  strut witnin a certain distance
% two types: {box, lattice}
switch VoxelType
    case 'box'
        [voxel] = CreateVoxelBox(voxel_c,n,p);
    otherwise
        [node,strut] = ReadStrut(); % get the information of strut
        for i = 1:length(voxel_c)  % for each voxel, deside if it is active
            % for each strut, get the distance to the voxel
            for j = 1:length(strut)
                start_n = node(strut(j,1),:);  % start node coordinate
                end_n = node(strut(j,2),:);    % end node coordinate
                center = voxel_c(i,4:6);        % voxel center position
                % determine alpha and beta are acute angle
                alpha = acosd(((center - start_n)*(end_n - start_n)')...
                    /(norm(center - start_n)*norm(end_n - start_n)));
                beta = acosd(((center - end_n)*(start_n - end_n)')...
                    /(norm(center - end_n)*norm(start_n - end_n)));
                if(alpha<90 && beta<90)% if not acute angle, distance to line
                    distance = norm(cross(end_n - start_n,center - start_n))...
                        /norm(end_n - start_n);
                else                      % if it is acute angle, distance to node
                    distance = min(norm(center - start_n),norm(center - end_n));
                end
                if (distance<=p)     % if distance less than radius, active it
                    voxel(voxel_c(i,1),voxel_c(i,2),voxel_c(i,3)) = 1;
                    continue;             % run for the next voxel
                end
            end
        end
end
Density = sum(sum(sum(voxel)))/n^3;% calculate the relative density
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Import information of strut
function [nodelist,strutlist] = ReadStrut()
    nodelist= [0 0 0; 1 0 0; 0 1 0;...
        0 0 1; 1 0 1; 0 1 1; 1 1 0; 1 1 1]; % equiv. ff2n(3)
    strutlist = [ones(1,3),4,4,3,3,6,2,2,5,7;2:7,6,8,7,5,8,8]'; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [voxel] = CreateVoxelBox(voxel_c,n,p)
% voxel_c: Matrix containing the central coordinates of the unit cell FE
% n: the number of voxel along each axis
% p: r or alpha/2 
    voxel=zeros(n^3,1);
    idx= voxel_c(:,4:6)>.5-p & voxel_c(:,4:6)<.5+p; 
    voxel(sum(idx,2)==3)=1;voxel=reshape(voxel,[n,n,n]);
end