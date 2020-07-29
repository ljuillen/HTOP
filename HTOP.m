function HTOP(nelx,nely,nelz,volfrac,CH,xlimit,rmin,problem)
% Description: Solves the Homogenization-based Topology Optimization Problem
%-----------------------------------------------------------------%
% ---- Input Arguments ----
% [nelx nely nelz]: Number of elements in x, y, z direction
% volfrac: Volume fraction
% CH: Elasticity tensor (symbolic format)
% xlimit: Extreme values of relative density array
% rmin: Filter radius
% problem: Problem type, 1 or 2
%-----------------------------------------------------------------%
clc; close all;
set(0,'defaultTextInterpreter','latex');
warning('off','all'); 
% USER-DEFINED LOOP PARAMETERS
tic; colordef black; % Start counting
maxloop = 40;    % Maximum number of iterations
tolx = 1e-4;     % Terminarion criterion
displayflag = 1; % if 1, display geometry  
[fixeddof,loaddof]=ProblemBoundaries(nelx,nely,nelz,problem); % Select problem

% PREPARE FINITE ELEMENT ANALYSIS
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
% PREPARE FILTER
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely+j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2;
                        k = k+1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2+(k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
% INITIALIZE ITERATION
x=repmat(volfrac,[nely,nelx,nelz]);
xPhys = x; 
loop = 0;
change = 1;
minV = xlimit(1);
maxV = xlimit(2);
% OBTAIN STIFFNESS MATRIX KE
[KE, DKE]= element_stiffness(0.5,0.5,0.5, CH);
stp.Value=0; % display geometry option
if logical(displayflag)
clear stp;
idx=ones(size(x)); % Visualize the initial domain
set(gcf,'numbertitle','off','name','Homogenization based TOP','units','norm','pos',[.15 .15 .75 .75]);
display_3D(idx,'domain'); hold all;
stp=uicontrol(gcf,'Style','toggle','String','Stop HTOP','units','norm','pos',[.1 .1 .2 .12],...
    'FontName','Garamond','FontSize',20,'Backgroundcolor',repmat(.9020,[1,3]));
end
C=zeros([],1);
% START ITERATION
while (change > tolx) && (loop < maxloop) && (stp.Value==0) 
    loop = loop+1;
    % FE-ANALYSIS
    sK=zeros(24*24,nele);
    for i=1:nele
        sK(:,i)=reshape(KE(xPhys(i)),[],1); % Stiffness matrices
    end
    K=sparse(iK(:),jK(:),sK(:),ndof,ndof); K=(K+K')/2;  % [K] total
    try
        L = ichol(K(freedofs,freedofs)); % Incomplete Cholesky for the KE matrix
        U(freedofs,:) = pcg(K(freedofs,freedofs),F(freedofs,:),1e-10,5e3,L,L');
    catch
        U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    end
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    dc=zeros(nele,1); c=zeros(nele,1); % Vector initizalization
    for i=1:nele
        c(i)=(((U(edofMat(i,:))' *KE(xPhys(i))* U(edofMat(i,:))))); %Elastic energy stored at each element
        dc(i)=-(U(edofMat(i,:))' *DKE(xPhys(i))* U(edofMat(i,:)));  %Sensitivity analysis
    end
    % FILTERING AND MODIFICATION OF SENSITIVITIES
    dc(:) = H*(min(dc(:),0)./Hs);
    % OPTIMALITY CRITERIA UPDATE
    l1 = 0; l2 = 1e12;
    move = .1;
    dL= -(volfrac*nele).*dc;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        xnew = max(minV,max(x(:)-move,min(maxV,min(x(:)+move,(x(:).*sqrt(dL./lmid))))));
        xPhys= H*(xnew(:)./Hs);
        if sum(xPhys(:))> volfrac*nele, l1 = lmid; else l2 = lmid; end
    end
    % CURRENT COMPLIANCE VALUE
    C(loop)=sum(c(:));
    change = norm(xPhys(:)-x(:),inf);
    x = xnew;
    % PRINT THE RESULTS TO COMMAND WINDOW
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,sum(c(:)),sum(xPhys(:))/nele,change);
    % DISPLAY THE GEOMETRY AT CURRENT ITERATION
    if logical(displayflag)
    display_3D(reshape(xPhys,[nely,nelx,nelz]),'topology');
    legend({'$Initial\ Domain$',['$No\ Iter:\ ',num2str(loop),'$ ',sprintf('\n'),...
        ' $Max\ Change:.\ ',num2str(change),'$']},'interpreter','latex',...
        'units','norm','Position',[0.6905 0.8437 0.1911 0.0508],'Box','off','fontsize',18);
    drawnow; pause(.5);
    end
end
figure('color',ones(1,3)); colordef white;
title('\textbf{\itshape{{Compliance vs Iterations:. (OC Method)}}}','Interpreter','Latex','FontSize',20);
plot(1:loop,C,'Color',[1 17 181] ./ 255,'LineStyle','--'); hold on; 
legend({['Minimum Compliance Value:. ',num2str(min(C(:)))]},...
    'FontSize',15,'interpreter','latex');
xlabel('No. Iterations','interpreter','latex','FontSize',20);
ylabel('Compliance','interpreter','latex','FontSize',20);
set(gca,'Yscale','log'); grid minor;
toc; % Stop counting
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE,KEx] = element_stiffness(a, b, c, CH)
% Obtain the element's stiffness matrix and its first order derivative
%-----------------------------------------------------------------%
% ---- Output Arguments ----
% KE: Stiffness matrix 
% KEx: 1st order derivative of the stiffness matrix 
%
% ---- Input Arguments ----
% a = half the length of the element at x direction (multiplied by the scale(1) value)
% b = half the length of the element at y direction (multiplied by the scale(2) value)
% c = half the length of the element at z direction (multiplied by the scale(3) value)
% CH
%-----------------------------------------------------------------%
% Three Gauss points 
xx = [-sqrt(3/5), 0, sqrt(3/5)]; yy = xx; zz = xx;
ww = [5/9, 8/9, 5/9];
syms x real;
KE=zeros(24,24);
for ii = 1:length(xx)
    for jj = 1:length(yy)
        for kk = 1:length(zz)
            %integration point
            x = xx(ii); y = yy(jj); z = zz(kk);
            %stress strain displacement matrix
            qx = [ -((y-1)*(z-1))/8, ((y-1)*(z-1))/8, -((y+1)*(z-1))/8,...
                ((y+1)*(z-1))/8, ((y-1)*(z+1))/8, -((y-1)*(z+1))/8,...
                ((y+1)*(z+1))/8, -((y+1)*(z+1))/8];
            qy = [ -((x-1)*(z-1))/8, ((x+1)*(z-1))/8, -((x+1)*(z-1))/8,...
                ((x-1)*(z-1))/8, ((x-1)*(z+1))/8, -((x+1)*(z+1))/8,...
                ((x+1)*(z+1))/8, -((x-1)*(z+1))/8];
            qz = [ -((x-1)*(y-1))/8, ((x+1)*(y-1))/8, -((x+1)*(y+1))/8,...
                ((x-1)*(y+1))/8, ((x-1)*(y-1))/8, -((x+1)*(y-1))/8,...
                ((x+1)*(y+1))/8, -((x-1)*(y+1))/8];
            % Jacobian
            J = [qx; qy; qz]*[-a a a -a -a a a -a; -b -b b b -b -b b b;...
                -c -c -c -c c c c c]';
            qxyz = J\[qx;qy;qz];
            B_e = zeros(6,3,8);
            for i_B = 1:8
                B_e(:,:,i_B) = [qxyz(1,i_B)   0             0;
                    0             qxyz(2,i_B)   0;
                    0             0             qxyz(3,i_B);
                    qxyz(2,i_B)   qxyz(1,i_B)   0;
                    0             qxyz(3,i_B)   qxyz(2,i_B);
                    qxyz(3,i_B)   0             qxyz(1,i_B)];
            end
            B = [B_e(:,:,1) B_e(:,:,2) B_e(:,:,3) B_e(:,:,4) B_e(:,:,5)...
                B_e(:,:,6) B_e(:,:,7) B_e(:,:,8)];
            % Weight factor at this point
            weight = det(J)*ww(ii) * ww(jj) * ww(kk);
            % Element matrices
            KE = KE + weight * B' * CH * B; %symbolic format of KE
        end
    end
end
syms x real; % Symbolic variable
Kex=diff(KE,x); % First order derivative, KEx = dKE/dvuc
KE=matlabFunction(KE); % Symbolic to function, KE
KEx= matlabFunction(Kex); % Symbolic to function, KEx
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fixeddof,loaddof]=ProblemBoundaries(nelx,nely,nelz,problem)
switch problem
    case 1 %  Fixed Beam 
        [il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
        loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
        loaddof = 3*loadnid(:) - 1;                             % DOFs
        % USER-DEFINED SUPPORT FIXED DOFs
        [iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
        fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Node IDs
        fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % DOFs
    case 2 % Simply supported 
        il=nelx/2; jl=0;kl=nelz/2;
        loadnid=kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl);
        loaddof=3*loadnid(:)-1;
        iif=[0 0 nelx nelx]; jf=zeros(1,4);
        kf=[0 nelz 0 nelz];
        fixednid=kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf);
        fixeddof=[3*fixednid(:);3*fixednid(:)-1;fixednid(:)-2];
end
end