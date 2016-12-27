clc; close all; cos0_eps  = 1.0-1.0E-14;  cos90_eps = 1.0E-7; eps=1.0E-10; warning off MATLAB:colon:operandsNotRealScalar

load 'fourLayerMesh';

% generating source index to find initial element: if we already has a determinent list, just skip this procedure.
detcoord = zeros(length(elem),1); % determinent of all elements
for iel = 1:length(elem)
    nd = elem(iel,:); % connected node for (iel)-th element, (1x4)
    coord = [node(nd(1),1:3); node(nd(2),1:3); node(nd(3),1:3); node(nd(4),1:3)]; % (4x3), column represents x,y,z respectively.
    detcoord(iel) = det([[1 1 1 1]' coord]); % negative values
end

MC.source = [1 1 0.15]; % z = 0.15 is superficial layer (due to small fluctuation of the mesh coordinate), x = y = 1 is the center of the cube.
idx = findSourceElem(elem(:,1:4),node(:,1:3), MC.source', detcoord);        % serach initial
MC.sourceElem = idx;                                            % initial souce element for the pencil beam

MC.numPhotons = 1e5;    % # of photon for simulation
scale = 0;              % 0: log scale, 1: linear scale

zRes  = 0.01;           % default resolution of imagenary plane - 50 um (unit: mm)
xyRes = 0.01;
rRes = 0.01;
xMax = max(node(:,1)); yMax = max(node(:,2)); zMax = max(node(:,3));
xMin = min(node(:,1)); yMin = min(node(:,2)); zMin = min(node(:,3));
nIrr = zeros(ceil((xMax-min(node(:,1)))/xyRes)+1, ceil((yMax-min(node(:,2)))/xyRes)+1, ceil(zMax/zRes)+1);
maxIdx = size(nIrr);
area  = xyRes*xyRes;
onAxis = zeros(1, ceil(zMax/zRes)+1);

% Layer1
Med(1).mus = 11.5;
Med(1).mua = 0.1;
Med(1).n = 1.36;
Med(1).g = 0.87;
% Layer2
Med(2).mus = 45;
Med(2).mua = 0.1;
Med(2).n = 1.38;
Med(2).g = 0.78;
% Layer3
Med(3).mus = 23;
Med(3).mua = 0.2;
Med(3).n = 1.36;
Med(3).g = 0.87;
% Layer4
Med(4).mus = 90;
Med(4).mua = 0.2;
Med(4).n = 1.38;
Med(4).g = 0.78;
% Layer 5 (up layer)
Med(5).mus = 10;
Med(5).mua = 0.01;
Med(5).n = 1.36;
Med(5).g = 1.0;
% Layer 6 (bottom layer)
Med(6).mus = 10;
Med(6).mua = 0.01;
Med(6).n = 1.38;
Med(6).g = 1.0;

hitBoundary = 0;
photon.Dead = 0;

tic;
fprintf('# of photon : %d\n', MC.numPhotons);

for iPhoton=1:MC.numPhotons
    
    photon.Dead = 0;    
    photon.CurElem = MC.sourceElem;
    photon.CurMed  = elem(photon.CurElem, 5);     
    
    photon.Pos=MC.source;
    photon.Dir=[0 0 1]; % Simulate Pencil Beam (No angular spreading)
     
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
                    
                xImgId = ceil(real(photonPlane(1))/xyRes);
                yImgId = ceil(real(photonPlane(2))/xyRes);
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
                
                xImgId = ceil(real(photonPlane(1))/xyRes);
                yImgId = ceil(real(photonPlane(2))/xyRes);
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
% nIrr = (nIrr./MC.numPhotons)./area;
Red = [237 28 36]./255; Yellow = [253 185 19]./255; Sky = [0 173 220]./255; Blue = [82 79 161]./255; Green = [203 219 42]./255; Gray = [180 180 180]./255;


map0 = squeeze(nIrr(100,:,:))';
map1 = squeeze(nIrr(101,:,:))';
map = map0+map1;

offset = round(0/xyRes);
map_offset = (map(:,100+offset)+map(:,101+offset));

drectAvg = map_offset(16:end);
drectAvg = drectAvg(1:100);
direct = (drectAvg)/drectAvg(1);

zp = [0:99]*0.01;

h = figure(1); 
semilogy(zp, direct, 'Color', Red, 'LineWidth', 2.5); hold on;
hold off;

set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 21, 'FontWeight', 'bold')
xlabel('Depth {\itz} (mm)');
ylabel('Relative Intensity (A.U.)');
title('Off-axis {\mum}');
% set(gca, 'YTick', [10^-4 10^-2 10^-0]);
set(gca, 'XTick', [0.2 0.5 0.7 1.0]);
grid on;  box off;
set(h, 'Position', [100, 100, 450, 350]);

ryAxis = [10:190]*0.01-1;
zAxis = [1:91]*0.01;
mapCut = map(16:90+16,10:190);
h2 = figure(2); 
imagesc(ryAxis, zAxis, log10(mapCut./drectAvg(1))); 
% colorbar;
set(gcf,'color','w');
set(gca,'FontName', 'Times New Roman', 'FontSize', 21, 'FontWeight', 'bold')
ylabel('Depth {\itz} (mm)');
xlabel('Transverse Distance {\ity} (mm)');
title('Direct photon flux (3-D)');
grid off;  box off;
set(h2, 'Position', [100, 100, 450, 350]);
caxis([-11, 0])
