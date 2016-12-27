%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MMCLAB - Mesh-based Monte Carlo for MATLAB/Octave by Qianqina Fang
%
% In this example, we show the most basic usage of MMCLAB.
%
% This file is part of Mesh-based Monte Carlo (MMC) URL:http://mcx.sf.net/mmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prepare simulation input
tic;
load('fourLayer_2mm.mat');
cfg.nphoton=1e6;
% [cfg.node face cfg.elem]=meshabox([0 0 0],[60 60 30],6);

cfg.node = node(:,1:3);
face = face(:,1:3);
cfg.elem =  elem;

cfg.elemprop=ones(size(cfg.elem,1),1);
cfg.srcpos=[1 1 0.15];
cfg.srcdir=[0 0 1];

% [mua mus g n]
cfg.prop=[0.01 10 1.0 1.38;  % Up layer
0.1 11.5 0.87 1.36;          % layer1
0.1 45 0.78 1.38;            % layer2
0.2 23 0.87 1.36;            % layer3
0.2 90 0.78 1.38;            % layer4
0.01 10 1.0 1.38             % Bottom layer
];

cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.debuglevel='TP';

% run the simulation

flux=mmclab(cfg);

% plotting the result

% if you have the SVN version of iso2mesh, use the next line to plot:
% qmeshcut(cfg.elem(:,1:4),cfg.node(:,1:3),log10(abs(flux.data(:))),'y=30','linestyle','none');
intensity = abs(flux.data(:));
cutElem = cfg.elem(elem(:,5) ~= 5 & elem(:,5) ~= 6,1:4);

% cfg.elem(:,1:4)
% [cutpos,cutvalue,facedata]=qmeshcut(cutElem,cfg.node(:,1:3),intensity,'y=1');
cfg.node(:,3) = cfg.node(:,3) - 0.2;
cfg.node(:,1) = cfg.node(:,1) - 1;
hcut = figure(1)
[cutpos,cutvalue,facedata]=qmeshcut(cutElem,cfg.node(:,1:3),intensity,'y=1');
cutvalue = cutvalue/max(cutvalue);
patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',log10(cutvalue),'facecolor','flat','linestyle','none');
view([0.0001 1 0]);
set(gca,'dataaspectratio',[1 1 1]);
set(hcut, 'Position', [100, 100, 400, 300]);
caxis([-11, 0])
zlim([0 0.9])
xlim([-0.8 0.8])
colorbar;
axis equal;
set(gca,'xdir','reverse')
set(gca,'zdir','reverse')
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
