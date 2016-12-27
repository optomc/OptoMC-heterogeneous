% load 'Simulated data\1e6 num of photons\M2_IntensityMap_1e6.mat'
% load 'Simulated data\1e6 num of photons\M2_IntensityMap_1e5.mat'
% load 'Simulated data\1e6 num of photons\M1_IntensityMap_1e6.mat'
% load 'Simulated data\1e6 num of photons\M1_IntensityMap_1e5.mat'

load 'Simulated data\1e6 num of photons\M2_AbsorptionMap_1e6.mat'
% load 'Simulated data\1e6 num of photons\M2_AbsorptionMap_1e5.mat'
% load 'Simulated data\1e6 num of photons\M1_AbsorptionMap_1e6.mat'
% load 'Simulated data\1e6 num of photons\M1_AbsorptionMap_1e5.mat'

%% Coronal (AP)
hc = 1:4;
for figNum = 1:4
    if ishandle(hc(figNum))
        close(hc(figNum))
    end
end
xAdd = 0.5;
xCutPos = MC.source(1) + xAdd;
plane = [xCutPos 0 0; xCutPos 0 1; xCutPos 1 0];

[cutPos, faceValue, faceData, faceElemId]=irrMeshCut(elem(:,1:4),node(:,1:3),elem(:,5),plane);

% subplot(121)
hc(1) = figure(21); axis on; grid on; hold on;
colormap hsv;
plotmesh(node(:,1:3), elem(elem(:,5) == 3,:), 'edgecolor', 'none', 'facecolor', [0.4980 0 0], 'facealpha',0.1);
plotmesh(node(:,1:3), elem(elem(:,5) == 2,:), 'edgecolor', 'none', 'facecolor', [0.5608 1.0000 0.4353], 'facealpha',0.1);
plotmesh(node(:,1:3), elem(elem(:,5) == 1,:), 'edgecolor', 'none', 'facecolor', [0 0 0.5608], 'facealpha',0.1);
plot3(MC.source(1),MC.source(2),MC.source(3),'+r');
patch('Faces',faceData,'Vertices',cutPos,'CData',faceValue,'facecolor','flat', 'facealpha', 0.5, 'edgecolor', 'none');
set(hc(1), 'Position', [660,200,xMax/xyRes*2.5740,yMax/xyRes*1.225]);
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
set(gca,'zdir','reverse');

% making mask
opengl software;
[cutPos,faceValue,faceData, faceElemId]=irrMeshCut(elem(:,1:4),node(:,1:3),elem(:,5),plane);
hc(2) = figure(22);
patch('Faces',faceData,'Vertices',cutPos,'CData',faceValue,'facecolor','flat', 'edgecolor', 'none'); axis off; hold on;
view(90, 0);
ylim([-yMax yMax]);
zlim([0 zMax]); axis off; grid off;
set(hc(2), 'Position', [700,200,yMax/xyRes*2.5740,zMax/zRes*1.225]);
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
set(gca,'zdir','reverse');

% making outline
F = getframe(gca);
opengl hardware;

mask = F.cdata;
hc(3) = figure(23); imagesc(mask); axis off; grid off; hold on;
set(hc(3), 'Position', [700,200,yMax/xyRes*2.5740,zMax/zRes*1.225]);
set(hc(2), 'Position', [960,220,yMax/xyRes*1.2870,zMax/zRes*0.6125]);
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
close figure 23;

hc(4) = figure(24);
set(hc(4), 'Position', [700,550,yMax/xyRes*2.974,zMax/zRes*1.225]); 
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
ylim([0 round(yMax/xyRes)]);
zlim([0 round(zMax/zRes)]);
gMask=rgb2gray(mask);
bEdges = (edge(gMask,'canny'));
bEdges = abs(bEdges-1);

sIrr = squeeze(nIrr(round((xCutPos-min(node(:,1)))/xyRes)+1,:,:))'; 
% padding data
if size(bEdges,2) < size(sIrr,2)
    pad = zeros(1, size(bEdges,1))';
    bEdges = [bEdges pad];
elseif size(bEdges,2) > size(sIrr,2)
    pad = zeros(1, size(bEdges,1))';
    sIrr = [sIrr pad];
end
if size(bEdges,1) < size(sIrr,1)
    pad = zeros(1, size(bEdges,2));
    bEdges = [bEdges; pad];
elseif size(bEdges,1) > size(sIrr,1)
    pad = zeros(1, size(bEdges,2));
    sIrr = [sIrr; pad];
end

sy = [round(yMin/xyRes):round(yMax/xyRes)]/20;
sz = [0:round(zMax/zRes)]/20;

if scale == 1
    scaleIrr = fiberTipPower*sIrr;
else
    scaleIrr = log(fiberTipPower*sIrr);    
end

h = imagesc(sy, sz, scaleIrr);  %log/linear
set(h, 'AlphaData', bEdges(1:size(sIrr,1),1:size(sIrr,2))); colorbar;
set(gca,'ydir','reverse')
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
xlim([round(yMin/xyRes)/20 round(yMax/xyRes)/20]);
zlim([0 round(zMax/zRes)/20]);

titleStr = ['x = ' num2str(MC.source(3) + xAdd) ' mm (AP)'];
title(titleStr)

figure(22); grid on; axis on;
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
plot3(MC.source(1),1.5,MC.source(3),'+r');

%% Sagital (ML)
hs = 1:4;
for figNum = 1:4
    if ishandle(hs(figNum))
        close(hs(figNum))
    end
end
yAdd = 0;
yCutPos = MC.source(2) - yAdd;
plane = [0 yCutPos 0; 0 yCutPos 1; 1 yCutPos 0;];
[cutPos,faceValue,faceData, faceElemId]=irrMeshCut(elem(:,1:4),node(:,1:3),elem(:,5),plane);
 
% subplot(121)
hs(1) = figure(31); axis on; grid on; hold on;
colormap hsv;
plotmesh(node(:,1:3), elem(elem(:,5) == 3,:), 'edgecolor', 'none', 'facecolor', [0.4980 0 0], 'facealpha',0.1);
plotmesh(node(:,1:3), elem(elem(:,5) == 2,:), 'edgecolor', 'none', 'facecolor', [0.5608 1.0000 0.4353], 'facealpha',0.1);
plotmesh(node(:,1:3), elem(elem(:,5) == 1,:), 'edgecolor', 'none', 'facecolor', [0 0 0.5608], 'facealpha',0.1);
plot3(MC.source(1),MC.source(2),MC.source(3),'+r');
patch('Faces',faceData,'Vertices',cutPos,'CData',faceValue,'facecolor','flat', 'facealpha', 0.5, 'edgecolor', 'none');
set(hs(1), 'Position', [1520,200,(xMax-xMin)/xyRes*0.6435,yMax/xyRes*1.225]);
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
set(gca,'zdir','reverse');

% making mask
opengl software;
[cutPos,faceValue,faceData, faceElemId]=irrMeshCut(elem(:,1:4),node(:,1:3),elem(:,5),plane);
hs(2) = figure(32);
patch('Faces',faceData,'Vertices',cutPos,'CData',faceValue,'facecolor','flat', 'edgecolor', 'none'); axis off; hold on;
view(0, 360);
xlim([xMin xMax]);
zlim([0 zMax]); axis off; grid on;
set(hs(2), 'Position', [1000,200,(xMax-xMin)/xyRes*1.287,zMax/zRes*1.225]);
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
set(gca,'zdir','reverse');

% making outline
F = getframe(gca);
opengl hardware;
mask = F.cdata;
hs(3) = figure(33); imagesc(mask); axis off; grid off; hold on;
set(hs(3), 'Position', [1200,210,(xMax-xMin)/xyRes*1.287,zMax/zRes*1.225]);
set(hs(2), 'Position', [1200,210,(xMax-xMin)/xyRes*0.6435,zMax/zRes*0.6125]);
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
close figure 33;

hs(4) = figure(34);
set(hs(4), 'Position', [1200,550,(xMax-xMin)/xyRes*1.487,zMax/zRes*1.225]); 
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
set(gcf,'color','w');
set(gca,'zdir','reverse');
xlim([0 round((xMax-xMin)/xyRes)]);
zlim([0 round(zMax/zRes)]);
gMask=rgb2gray(mask);
bEdges = (edge(gMask,'canny'));
bEdges = abs(bEdges-1);

sIrr = squeeze(nIrr(:,round((yCutPos-min(node(:,2)))/xyRes)+1,:))'; 
% padding data
if size(bEdges,2) < size(sIrr,2)
    pad = zeros(1, size(bEdges,1))';
    bEdges = [bEdges pad];
elseif size(bEdges,2) > size(sIrr,2)
    pad = zeros(1, size(bEdges,1))';
    sIrr = [sIrr pad];
end
if size(bEdges,1) < size(sIrr,1)
    pad = zeros(1, size(bEdges,2));
    bEdges = [bEdges; pad];
elseif size(bEdges,1) > size(sIrr,1)
    pad = zeros(1, size(bEdges,2));
    sIrr = [sIrr; pad];
end

tx = [round(xMin/xyRes):round(xMax/xyRes)]/20;
tz = [0:round(zMax/zRes)]/20;

if scale == 1
    scaleIrr = fiberTipPower*sIrr;
else
    scaleIrr = log(fiberTipPower*sIrr);    
end

h = imagesc(tx, flipud(tz), scaleIrr);  %log/linear
set(h, 'AlphaData', bEdges(1:size(sIrr,1),1:size(sIrr,2))); colorbar;
set(gca,'ydir','reverse')
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
xlim([round(xMin/xyRes)/20 round(xMax/xyRes)/20]);
zlim([0 round(zMax/zRes)/20]);

titleStr = ['y = ' num2str(MC.source(2) - yAdd) ' mm (ML)'];
title(titleStr);

figure(32); grid on; axis on;
set(gcf,'color','w');
plot3(MC.source(1),MC.source(2),MC.source(3),'+r');

%% Transverse (DV)
ht = 1:4;
for figNum = 1:4
    if ishandle(ht(figNum))
        close(ht(figNum))
    end
end
zAdd = 0.5;
zCutPos = MC.source(3) + zAdd; % target plane: 500 um away from the stimulation site
plane=[0 0 zCutPos; 1 0 zCutPos; 0 1 zCutPos];

[cutPos,faceValue,faceData, faceElemId]=irrMeshCut(elem(:,1:4),node(:,1:3),elem(:,5),plane);

% subplot(121)
ht(1) = figure(11); axis on; grid on; hold on;
colormap hsv;
plotmesh(node(:,1:3), elem(elem(:,5) == 3,:), 'edgecolor', 'none', 'facecolor', [0.4980 0 0], 'facealpha',0.1);
plotmesh(node(:,1:3), elem(elem(:,5) == 2,:), 'edgecolor', 'none', 'facecolor', [0.5608 1.0000 0.4353], 'facealpha',0.1);
plotmesh(node(:,1:3), elem(elem(:,5) == 1,:), 'edgecolor', 'none', 'facecolor', [0 0 0.5608], 'facealpha',0.1);
patch('Faces',faceData,'Vertices',cutPos,'CData',faceValue,'facecolor','flat', 'facealpha', 0.5, 'edgecolor', 'none');
set(ht(1), 'Position', [10,200,(xMax-xMin)/xyRes*0.6435,yMax/xyRes*1.225]);
set(gcf,'color','w');
set(gca,'zdir','reverse');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')

% making mask
opengl software;
[cutPos,faceValue,faceData, faceElemId]=irrMeshCut(elem(:,1:4),node(:,1:3),elem(:,5),plane);
ht(2) = figure(12);
patch('Faces',faceData,'Vertices',cutPos,'CData',faceValue,'facecolor','flat', 'edgecolor', 'none'); axis off; hold on;
xlim([xMin xMax]);
ylim([yMin yMax]); axis off; grid off;
set(ht(2), 'Position', [10,200,(xMax-xMin)/xyRes*1.287,yMax/xyRes*2.45]);
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
% making outline
F = getframe(gca);
opengl hardware;

mask = F.cdata;
ht(3) = figure(13); imagesc(mask); axis off; hold on;
set(ht(3), 'Position', [330,200,(xMax-xMin)/xyRes*1.287,yMax/xyRes*2.45]);
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
set(ht(2), 'Position', [330,200,(xMax-xMin)/xyRes*0.6435,yMax/xyRes*1.225]);
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
close figure 13;

ht(4) = figure(14);
set(ht(4), 'Position', [10,500,(xMax-xMin)/xyRes*1.487,yMax/xyRes*2.45]); 
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
xlim([round(xMin/xyRes) round(xMax/xyRes)]);
ylim([round(yMin/xyRes) round(yMax/xyRes)]);
gMask=rgb2gray(mask);
bEdges = (edge(gMask,'canny'));
bEdges = abs(bEdges-1);

sIrr = fliplr(rot90(nIrr(:,:,zCutPos/zRes+1),2)');

% padding data
if size(bEdges,2) < size(sIrr,2)
    pad = zeros(1, size(bEdges,1))';
    bEdges = [bEdges pad];
elseif size(bEdges,2) > size(sIrr,2)
    pad = zeros(1, size(bEdges,1))';
    sIrr = [sIrr pad];  
end
if size(bEdges,1) < size(sIrr,1)
    pad = zeros(1, size(bEdges,2));
    bEdges = [bEdges; pad];
elseif size(bEdges,1) > size(sIrr,1)
    pad = zeros(1, size(bEdges,2));
    sIrr = [sIrr; pad];   
end

rx = [round(xMin/xyRes):round(xMax/xyRes)-1]/20;
ry = [round(yMax/xyRes)-1:-1:round(yMin/xyRes)]/20;

if scale == 1
    scaleIrr = fiberTipPower*sIrr;
else
    scaleIrr = log(fiberTipPower*sIrr);    
end

h = imagesc(rx, ry, scaleIrr);  %log/linear
set(h, 'AlphaData', bEdges); colorbar;
set(gca,'ydir','normal')
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')
xlim([round(xMin/xyRes)/20 round(xMax/xyRes)/20]);
ylim([round(yMin/xyRes)/20 round(yMax/xyRes)/20]);

titleStr = ['z = ' num2str(-MC.source(3) - zAdd) ' mm (DV)'];
title(titleStr);

figure(12); grid on; axis on;
plot3(MC.source(1),MC.source(2),MC.source(3),'+r');
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')

% %%
% ex = figure(1);
% plotmesh(node(:,1:3), elem(faceElemId,:), 'edgecolor', [0 0 0.4980], 'edgealpha',0.1, 'facecolor', [0 0 0.4980], 'facealpha',0.1); hold on;
% plot3(MC.source(1),MC.source(2),MC.source(3),'+r');
% set(gcf,'color','w'); grid on; axis on;
% set(gca,'FontName', 'Times New Roman', 'FontSize', 12, 'FontWeight', 'bold')
% set(ex, 'Position', [10,200,400,560]);