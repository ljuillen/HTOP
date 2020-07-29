function unit_cell_visualization(p,VoxelType,dens,voxelex)
% Decription: Visualize the unit cells
% ---- Input Arguments ----
% p: geometric parameter vector
% VoxelType: 'lattice' or 'box'
% dens: 1st material relative density
% voxelex: Voxel matrix (logical)
close all; colordef white; close(gcf);
step=floor(length(p)/3); % set step to display
% Prepare figure
figure('name','Unit Cell visualization','units','norm','pos',[0.2 0.2 0.6 0.6]);
% Annotations
dim1=[0.052 0.75 0.03 0.022];
annotation('textbox',dim1,'String','$ Material_{2} $','BackgroundColor',[.6974, .4788,.5856],...
    'FitBoxToText','on','FaceAlpha',0.8,'interpreter','latex');
p_disp=p(1:step:end);
dim2=[0.052 0.86 0.03 0.022];
annotation('textbox',dim2,'String','$ Material_{1} $','BackgroundColor',[.2 .2 .2],...
    'FitBoxToText','on','FaceAlpha',0.8,'interpreter','latex');
if strcmpi(VoxelType,'lattice')
    param='radius';
else
    param='\frac{\alpha}{2}';
end
set(gca,'pos',[0.1300 0.1100 0.7750 0.8150]); %Prepare axis
% Plot unit cells
j=0;
for i=1:step:length(p)
    j=j+1;
    subplot(ceil(sqrt(length(p_disp))),ceil(sqrt(length(p_disp))),j);
    vlx=voxelex{1,i}; % Current voxel matrix
    display_3D(vlx,'unitcell'); % display voxel
    title(['$Relative Density_{1} = ',...
        num2str(dens(1,i)),'$',sprintf('\n'),'$',char(param),'= ',num2str(double(p(i))),'\ [Half\ section]$'],...
        'interpreter','latex','FontSize',15); drawnow; pause(.5);
end
hold all;
% Continue to the FittingProcess function.
uicontrol(gcf,'Style','push','String','Continue','ForegroundColor',...
    'k','Backgroundcolor',[0.9137    0.9490    0.9137],'Callback','close(gcf)','units',...
    'norm','pos',[0.052 0.5 0.12 0.1],...
    'FontName','Monotype Corsiva','FontSize',20);
uiwait(gcf);
end