clc; clear; close all;
%% run CreateUnitCellDataBase
[CHFunc,dens] = CreateUnitCellDataBase(20,'lattice',3,[0.5769,1e-2*.5769],...
     [0.3846,1e-2*.3846],'on');
%% run HTOP
HTOP(20,10,20,0.2,CHFunc,[dens(1),dens(end)],1.5,2);