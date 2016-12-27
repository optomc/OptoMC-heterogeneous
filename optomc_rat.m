clc; clear; close all; cos0_eps  = 1.0-1.0E-14;  cos90_eps = 1.0E-7; eps=1.0E-10; warning off MATLAB:colon:operandsNotRealScalar

load 'ratMesh';

% generating source index to find initial element: if we already has a determinent list, just skip this procedure.
detcoord = zeros(length(elem),1); % determinent of all elements
for iel = 1:length(elem)
    nd = elem(iel,:); % connected node for (iel)-th element, (1x4)
    coord = [node(nd(1),1:3); node(nd(2),1:3); node(nd(3),1:3); node(nd(4),1:3)]; % (4x3), column represents x,y,z respectively.
    detcoord(iel) = det([[1 1 1 1]' coord]); % negative values
end

MC.source = [1 -2.5 1.5];   % [x (AP), y (ML), z(-DV)]center of the optical fiber: define position of initial photon entry

plane = [0 0 MC.source(3); 0 1 MC.source(3); 1 0 MC.source(3)]; 
[faceElemId]=myMeshCutSourceElem(elem(:,1:4),node(:,1:3),plane);
idx = findSourceElem(elem(faceElemId,1:4),node(:,1:3), MC.source', detcoord);
MC.sourceElem = faceElemId(idx);         % initial souce element related with center of the optical fiber

MC.numPhotons = 1e5;    % # of photon for simulation
fiberTipPower = 8;      % unit: mW
scale = 0;              % 0: log scale, 1: linear scale

% Define fiber-optic properties (unit: mm)
MC.fiberWaist = 0.1; 
MC.fiberNA = 0.39;

zRes  = 0.05;           % default resolution of imagenary plane - 50 um (unit: mm)
xyRes = 0.05;
xMax = max(node(:,1)); yMax = max(node(:,2)); zMax = max(node(:,3));
xMin = min(node(:,1)); yMin = min(node(:,2)); zMin = min(node(:,3));
nIrr = zeros(ceil((xMax-min(node(:,1)))/xyRes)+1, ceil((yMax-min(node(:,2)))/xyRes)+1, ceil(zMax/zRes)+1);
maxIdx = size(nIrr);
area  = xyRes*xyRes;

% Gray matter (473nm)
Med(1).mus = 11.5;
Med(1).mua = 0.1;
Med(1).n = 1.36;
Med(1).g = 0.87;
% CSF (473nm)
Med(2).mus = 0.01;
Med(2).mua = 0.004;
Med(2).n = 1.37;
Med(2).g = 0.9;
% White matter (473nm)
Med(3).mus = 45;
Med(3).mua = 0.2;
Med(3).n = 1.38;
Med(3).g = 0.78;

hitBoundary = 0;
photon.Dead = 0;

tic;
fprintf('# of photon : %d\n', MC.numPhotons);

for iPhoton=1:MC.numPhotons
    photon.Dead = 0;
    radius = sqrt(rand()*MC.fiberWaist^2);          %uniform transverse distribution (top-hat profile)
    theta = rand()*2*pi;
    x = radius*cos(theta) + MC.source(1);
    y = radius*sin(theta) + MC.source(2);
    
    photon.Pos=[x, y, MC.source(3)];
     
    % find current sourceElem within initial plane
    idx = findSourceElem(elem(faceElemId,1:4),node(:,1:3), photon.Pos', detcoord);
    MC.sourceElem = faceElemId(idx);
    photon.CurElem = MC.sourceElem;
    photon.CurMed  = elem(photon.CurElem, 5); 
    
    photon.Dir=launchEffNA(MC.fiberNA, Med(photon.CurMed).n);       %uniform angular distribution                
     
    maxDist = exprnd(1/Med(photon.CurMed).mua);
    totalDist = 0;
        
    if mod(iPhoton, 100) == 0
        disp(num2str(iPhoton));
    end
    
    while (photon.Dead ~= 1)            
        photon.PrevPos=photon.Pos;
        if (hitBoundary ~= 1)  % fist photon movement
            s=stepSize(Med(photon.CurMed).mus);		
        end        
        if (s+totalDist > maxDist)
            s = maxDist-totalDist;
        end
         
        [minDist, minNodes, minCos, minNormal] = rayTetrahedronIntersects(photon.CurElem, elem, node, photon.Pos, photon.Dir);
        if minDist == -1
            break;
        end
        faceNormal = minNormal;
       
        if (s > minDist)
            hitBoundary = 1;
            
            % Check photon termination due to absorption
            if (totalDist+minDist >= maxDist-eps)
                minDist = maxDist-totalDist;
                photon.Dead = 1;
            end                 
            
            % Move photon to boundary of the element
            photon.Pos=photon.Pos+minDist.*photon.Dir; %%%%%            
            
            % Find out imaginary plane between each step
            if (sign(photon.Dir(3))>0)
                % find closest and furthermost imaginary plane (photon direction: upward)
                zClosestPlane = ceil(photon.PrevPos(3)/zRes)*zRes;
                zFurtherPlane = floor(photon.Pos(3)/zRes)*zRes;
                zPlanePos = [zClosestPlane:zRes:zFurtherPlane];                
            else
                % find closest and furthermost imaginary plane (photon direction: downward)
                zClosestPlane = floor(photon.PrevPos(3)/zRes)*zRes;
                zFurtherPlane = ceil(photon.Pos(3)/zRes)*zRes;
                zPlanePos = [zClosestPlane:-zRes:zFurtherPlane];
            end
            
            for i=1:numel(zPlanePos)
                % ray/plane intersect for recording
                zPlane.Pos = [0 0 zPlanePos(i)];
                zPlane.Normal = [0 0 -1];
                t = rayPlaneIntersects(zPlane, photon.PrevPos, photon.Dir);   
                photonPlane = photon.PrevPos+t*photon.Dir;
                    
                xImgId = ceil(real(photonPlane(1)-min(node(:,1)))/xyRes);
                yImgId = ceil(real(photonPlane(2)-min(node(:,2)))/xyRes);
                zImgId = round(real(zPlanePos(i)/zRes));                     
                               
                if (zImgId < 1 || zImgId > maxIdx(3) || xImgId < 1 || xImgId > maxIdx(1) || yImgId < 1 || yImgId > maxIdx(2))
                    break; % out of matrix index
                end

                nIrr(xImgId, yImgId, zImgId) = nIrr(xImgId, yImgId, zImgId) + 1;
%                 nIrr(xImgId, yImgId, zImgId) = nIrr(xImgId, yImgId, zImgId) + Med(photon.CurMed).mua;                

            end

            totalDist = totalDist + minDist;
            s = s-minDist;
            sRawRest = s*Med(photon.CurMed).mus; %sRawRest = RAW rest step remained.
           
            % find next elements (transmission)
            tmp1 = ismember(adjElemNode(photon.CurElem,:), minNodes);
            tmp1 = reshape(tmp1,4,4);
            tmp2 = sum(tmp1(2:4,:));
            tmp3 = reshape(adjElemNode(photon.CurElem,:), 4, 4);
            tmp4 = tmp3(:, tmp2 == 3);
            if ~isempty(tmp4)
                nextElem  = tmp4(1);        %nextNodes = tmp4(2:4);
            else
                break;
            end

            % Reflection and transmission
            photon.PrevMed  = photon.CurMed;
            photon.CurMed   = elem(nextElem,5);  
            photon.PrevElem = photon.CurElem;
            photon.CurElem  = nextElem;
            
            if photon.CurMed ~= photon.PrevMed
                cosa = -minCos;
                s = sRawRest/Med(photon.CurMed).mus;        % update rest step.

                % Calculate rFresnel
                if cosa > cos0_eps        % normal incident
                    rFresnel = (Med(photon.PrevMed).n-Med(photon.CurMed).n)/(Med(photon.PrevMed).n+Med(photon.CurMed).n);
                    rFresnel = rFresnel.*rFresnel;
                    cosb = 1;
                elseif cosa < cos90_eps    % very slant, it doesn't need cosb (always rand() < rFresnel)
                    rFresnel = 1;
                else
                    sina = sqrt(1-cosa^2);       % general incident  
                    sinb = Med(photon.PrevMed).n * sina / Med(photon.CurMed).n;                       
                    cosb = sqrt(1-sinb^2);
                    cap = cosa*cosb - sina*sinb;
                    cam = cosa*cosb + sina*sinb;
                    sap = sina*cosb + cosa*sinb;
                    sam = sina*cosb - cosa*sinb;
                    rFresnel = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
                end

                if rand() <= rFresnel
                    % Calculate Reflection vector V = 2*cos(a)*N + U
                    photon.CurMed  = photon.PrevMed;
                    photon.CurElem = photon.PrevElem;
                    photon.Dir = 2*cosa*faceNormal + photon.Dir;
                else
                    % Calculate Transmit vector   V = -cos(b)*N + (sinb/sina)*(cosa*N+U), sinb/sina = n1/n2       
                    photon.Dir = -cosb*faceNormal + (Med(photon.PrevMed).n/Med(photon.CurMed).n)*(cosa*faceNormal+photon.Dir);   
                                     
                    % update maximum travel distance
                    restDistRaw = (maxDist - totalDist)*Med(photon.PrevMed).mua;
                    restDist = restDistRaw/Med(photon.CurMed).mua;
                    maxDist  = restDist + totalDist; % Updated maxDist as it passes through the medium
                end                    
            end
            continue;                    
        else
            % Move photon within the same element, and scatter photon.
            hitBoundary = 0;
            
            % Check photon termination due to absorption
            if (totalDist+s >= maxDist-eps)
                photon.Dead = 1;
            end                
            
            photon.Pos=photon.Pos+s.*photon.Dir; 
            
            % Find out imaginary plane between each step
            if (sign(photon.Dir(3))>0)
                % find closest and furthermost imaginary plane (photon direction: upward)
                zClosestPlane = ceil(photon.PrevPos(3)/zRes)*zRes;
                zFurtherPlane = floor(photon.Pos(3)/zRes)*zRes;
                zPlanePos = [zClosestPlane:zRes:zFurtherPlane];                
            else
                % find closest and furthermost imaginary plane (photon direction: downward)
                zClosestPlane = floor(photon.PrevPos(3)/zRes)*zRes;
                zFurtherPlane = ceil(photon.Pos(3)/zRes)*zRes;
                zPlanePos = [zClosestPlane:-zRes:zFurtherPlane];
            end
            
            for i=1:numel(zPlanePos)
                % ray/plane intersect for recording
                zPlane.Pos = [0 0 zPlanePos(i)];
                zPlane.Normal = [0 0 -1];
                t = rayPlaneIntersects(zPlane, photon.PrevPos, photon.Dir);   
                photonPlane = photon.PrevPos+t*photon.Dir;
                
                xImgId = ceil(real(photonPlane(1)-min(node(:,1)))/xyRes);
                yImgId = ceil(real(photonPlane(2)-min(node(:,2)))/xyRes);
                zImgId = round(real(zPlanePos(i)/zRes));                   

                if (zImgId < 1 || zImgId > maxIdx(3) || xImgId < 1 || xImgId > maxIdx(1) || yImgId < 1 || yImgId > maxIdx(2))
                    break; % out of matrix index
                end

                nIrr(xImgId, yImgId, zImgId) = nIrr(xImgId, yImgId, zImgId) + 1;
%                 nIrr(xImgId, yImgId, zImgId) = nIrr(xImgId, yImgId, zImgId) + Med(photon.CurMed).mua;      
                
            end
                        
            photon.Dir=spin(Med(photon.CurMed).g, photon.Dir);
            totalDist = totalDist + s; 
        end          
    end        
end

time=toc;
fprintf('Memory space for nIrr: %f [MB]\n', ((size(nIrr(:),1)*8)/1024)/1024);
fprintf('Episode time: %f [sec]\n',time);
fprintf('Average time per one photon: %f [sec]\n',time/MC.numPhotons);

% Calculate normalized irradiance
nIrr = (nIrr./MC.numPhotons)./area;
% nIrrHigh = (nIrrHigh./MC.numPhotons)./areaHigh;
%% Coronal (AP)
hc = 1:4;
for figNum = 1:4
    if ishandle(hc(figNum))
        close(hc(figNum))
    end
end
xAdd = 0;
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
plot3(MC.source(1),MC.source(2),MC.source(3),'+r');

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

titleStr = ['z = ' num2str(MC.source(3) + zAdd) ' mm (-DV)'];
title(titleStr);

figure(12); grid on; axis on;
plot3(MC.source(1),MC.source(2),MC.source(3),'+r');
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 15, 'FontWeight', 'bold')