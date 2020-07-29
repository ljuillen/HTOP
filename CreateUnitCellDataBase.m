function [CHFunc,dens] = CreateUnitCellDataBase(VoxelNumber,VoxelType,Pol_Order,Lambda,Mu,vis)
% CreateUnitCellDataBase Creates the unit cell data base
% Create unit cell data base with homogenized tensors & relative densities
% of the 1st material
%-----------------------------------------------------------------% 
% --- Input Arguments ---
% VoxelNumber: Unit cell number of voxel along each direction
% VoxelType:   Unit cell type ('lattice' or 'box')
% Pol_Order:   Order of the polynomial for the 'box' unit cell case
%              []  for the 'lattice' unit cell type
% Lambda:      Lame's first parameter for solid materials
% Mu:          Lame's second parameter for solid materials
% vis:         'on', 'off'
% --- Output Arguments ---
% CHFunc:  The sumbolic Homogenized elasticity tensor 
% dens:    The relative density vector       
%-----------------------------------------------------------------%
%% CREATE GEOMETRIC PARAMETER VECTOR 
p=UnitCellType(VoxelNumber,VoxelType);
%% HOMOGENIZATION 
[ch,dens]=Homogenization(VoxelNumber,Lambda,Mu,p,VoxelType,vis);
%% FITTING PROCESS
[CHFunc]=FittingProcess(ch,dens,VoxelType,Pol_Order,vis);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p]=UnitCellType(voxels,VoxelType)
% UnitCellType returns the geometric parameter vector p
% Creates either the vector r or the half length a/2 depending on the 
% VoxelType, 'lattice' or 'box'
%-----------------------------------------------------------------%
% ---- Input Arguments ----
% voxels: Unit cell number of voxel along each axis
% VoxelType:  Unit cell type ('lattice' or 'box')
%
% ---- Output Arguments ----
% p: Unit cell geometric parameter vector
%-----------------------------------------------------------------%
%% GENERATE VECTOR R
sz=1/voxels; % Length of each voxel.
switch VoxelType
    case 'box'
        % Generate the a/2 vector
        p=1e-5:sz:.5-sz/2+1e-3; 
    case 'lattice'
         % Generate the radius vector
        p=[1e-5,sqrt(2)*sz/2:sqrt(2)*sz:.5*sqrt(2)]; 
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ch,dens]=Homogenization(voxels,Lamda,Mu,p,VoxelType,vis)
% Homogenization computes the unit cell elasticity tensors & relative densities
% Compute homogenized elasticity tensors and corresponding
% densities for unit cells of different geometric parameters
%-----------------------------------------------------------------%
% ---- Input Arguments  ----
% voxels:       Unit cell number of voxels along each axis 
% Lamda:        Lame's first parameter for solid materials 
% Mu:           Lame's second parameter for solid materials 
% p:            Design Parameter vector
% VoxelType:  Unit cell type ('lattice' or 'box')
% vis:          {'on', 'off'} visualize or not the unit cells 
% ---- Output Arguments ----
% CH: Homogenized elasticity tensors for each parameter value
% dens: Relative Densities of the 1st material 
%-----------------------------------------------------------------%
set(0,'defaultTextInterpreter','latex'); 
%% INITIALIZE Vector and cell arrays
ch=cell(1,[]); % Store the elasticity tensors 
dens=zeros(1,length(p));  % Store the relative densities 
voxelex=cell(1,[]); % Representation of the voxels as a 3d matrix 
% parfor 
parfor ii=1:length(p)
    %% GENERATE VOXEL 3D MATRIX & DENSITY 
    [voxel,dens(ii)] = GenerateVoxel(voxels,VoxelType,p(ii));
    %% COMPUTE ELASTICITY TENSOR
    ch{1,ii} = homo3D(1,1,1,Lamda,Mu,voxel); % Store the equivalent tensors 
    if strcmpi(vis,'on')
        voxelex{1,ii}=voxel; % Store the voxel matrix
    end
end
% Visualize the Unit Cells
if strcmpi(vis,'on')
    unit_cell_visualization(p,VoxelType,dens,voxelex); % call function
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CHf]=FittingProcess(ch,dens,VoxelType,Pol_Order,vis)
% Description: Fit the obtained CH w.r.t. the Relative Density terms 
%-----------------------------------------------------------------%
% --- Input Arguments ---
% ch:     Elasticity tensor cell array
% dens:   Relative density terms
% VoxelType:  'lattice' or 'box'
% Pol_Order:  Polynomial order for box unit cells 
% vis:    {'on', 'off'} visualize or not the fitting results
% --- Output Arguments ---
% CHf:  Symbolic format of the elasticity tensor w.r.t. the densities
%-----------------------------------------------------------------%
clc; close all;
set(0,'defaultTextInterpreter','latex');% Default Intrepreter
% Color palette
colour_darkblue = [1 17 181]./ 255; colour_peach = [251 111 66]./ 255;
colour_green = [12 195 82]./ 255; colour_teal = [18 150 155]./ 255;
C={colour_teal,colour_darkblue,colour_peach,colour_green};
minV=1e-5; % Threshold, minimum Value of the design variable
ch=double(horzcat(ch{:})); % ch [6 x 6*n] matrix
CH11=ch(1,1:6:end); % CH11 terms
CH12=ch(1,2:6:end); % CH12 terms
CH22=ch(2,2:6:end); % CH22 terms
CH66=ch(6,6:6:end); % CH66 terms
index = {'11','12','66','22'}; % index terms
% Prepare the data;  
switch VoxelType 
    case 'lattice' % sigmoid fit;   
        ft = fittype( 'c/(1+exp(-a*(x-b)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares','Display','Off');
    case 'box' % Polynomial fit of Pol_Order orde
        if isnumeric(Pol_Order)
        options = fitoptions('Method','LinearLeastSquares','lower',zeros(1,3)); 
        eval(['ft=fittype("poly',num2str(Pol_Order),'");']);
        dummy=[1:Pol_Order;flip(1:Pol_Order)];
        for ii =1:length(index)
            eval(['modelfun',char(index{ii}),'= @(b,x)',sprintf(repmat(' b(%d)*x.^%d +',[1,Pol_Order]),dummy),'CH',char(index{ii}),'(1);']);
            eval(['f',char(index{ii}),'="fun',char(index{ii}),'=',sprintf(repmat(' b0(%d)*x^%d +',[1,Pol_Order]),dummy),'CH',char(index{ii}),'(1);";']);
        end
        end
end
% Fit the data
syms x,'real';
for ii=1:length(index)
    eval(['[x',(char(index{ii})),', y',(char(index{ii})),']=prepareCurveData(dens,CH',char(index{ii}),');']);
    if strcmpi(VoxelType,'lattice')
        opts.StartPoint = zeros(1,3);
        eval(['[fitresult',char(index{ii}),',~]=fit( x',char(index{ii}),', y',char(index{ii}),',ft, opts);']);
        eval(['cf',char(index{ii}),'=coeffvalues(fitresult',char(index{ii}),');']);
        eval(['fun',char(index{ii}),'=vpa(cf',char(index{ii}),'(3)/(1+exp(-cf',char(index{ii}),'(1)*(x-cf',char(index{ii}),'(2)))));']);
        eval(['[a,b,c] =deal(cf',char(index{ii}),'(1),cf',char(index{ii}),'(2),cf',char(index{ii}),'(3));']);
        eval(['leg',char(index{ii}),'=["$ CH^{H}_{',char(index{ii}),'}\ =\ \frac{',num2str(c),'}{1\ +\ e^{-',num2str(a),'\ \cdot\ (v_{uc}\ -',num2str(b),')}}$"];']);
        clear a b c 
    elseif strcmpi(VoxelType,'box')
        eval(['[fitresult',char(index{ii}),',~]=fit( x',char(index{ii}),', y',char(index{ii}),',ft, options);']);
        eval(['cf',char(index{ii}),'=coeffvalues(fitresult',char(index{ii}),');']);
        eval(['b0=cf',char(index{ii}),'(1:Pol_Order);']);
        eval(['fnl',char(index{ii}),'=fitnlm(dens,CH',char(index{ii}),'(:),modelfun',char(index{ii}),',b0);']);
        eval(['b0=table2array(fnl',char(index{ii}),'.Coefficients(:,1));']) % for the symbolic part
        eval(eval(['f',char(index{ii})]));
        eval(['leg',char(index{ii}),'=["$ Modelfun\ with\ Rsq:.',num2str(eval(['fnl',char(index{ii}),'.Rsquared.Adjusted'])),' $"];'])
        clear b0
    end
end
% Symbolic format of the CH tensor 
syms x real; digits(4);
CHf=vpa([fun11,fun12,fun12,0,0,0;
         fun12,fun22,fun12,0,0,0;
         fun12,fun12,fun22,0,0,0;
         0,0,0,fun66,0,0;
         0,0,0,0,fun66,0;
         0,0,0,0,0,fun66]);

% Plot the fitting data
if strcmpi(vis,'on')
param='v_{uc}';
for ii=1:length(index)
    subplot(2,2,ii)
    scatter(dens,eval(['CH',char(index{ii}),'(:)']),'ks','MarkerFaceColor','y');
    hold on;
    fplot(matlabFunction(eval(['fun',char(index{ii})])),[min(dens),max(dens)],...
    'Color',C{1},'LineStyle','-.');
    title(['$ CH^{H}_{',char(index{ii}),'}',' - ',param,' $'],'Interpreter','latex','FontSize',20);
    xlabel(['$ ',param,' $'],'Interpreter','latex','FontSize',15); 
    xlim([min(dens),max(dens)]);
    ylabel(['$ CH^{H}_{',char(index{ii}),'( ',param,')} $'],'Interpreter','latex','FontSize',15);
    legend({'Scattered Data',eval(['leg',char(index{ii})])},'Location','best',...
    'fontSize',15,'Interpreter','latex'); grid minor; grid on; drawnow;   
end
uicontrol(gcf,'Style','push','String','HTOP','unit','norm',...
    'pos',[.9 .45 .1 .1],'ForegroundColor','k','Backgroundcolor',...
    [0.7647    0.8706    0.7647],'Callback', 'close all',...
    'FontName','Monotype Corsiva','FontSize',20);
uiwait(gcf);  % Continue to HTOP
end
end