function [sourceidx]=findSourceElem(nodes, gcoord, Rq, detcoord)
% calculate element index that includes initial photon.

nel = size(nodes,1); % # of elements
sourceidx = 0;
dummy = [1 1 1 1]';
denominator = detcoord;
nominator = zeros(4,1);

for iel = 1 : nel
    nd = nodes(iel,:); % connected node for (iel)-th element, (1x4)
    nominator(1) = det([dummy [Rq(:)'; gcoord(nd(2),:); gcoord(nd(3),:); gcoord(nd(4),:)]]);
    nominator(2) = det([dummy [gcoord(nd(1),:); Rq(:)'; gcoord(nd(3),:); gcoord(nd(4),:)]]);
    nominator(3) = det([dummy [gcoord(nd(1),:); gcoord(nd(2),:); Rq(:)'; gcoord(nd(4),:)]]);
    nominator(4) = det([dummy [gcoord(nd(1),:); gcoord(nd(2),:); gcoord(nd(3),:); Rq(:)']]);
    % determinent of nominator has minus sign when source is inside due to dummy is located left

    judge = nominator./denominator(iel);

    if(sum(judge(:) >= 0) == 4 || sum(judge(:) == 0) >= 1)
        sourceidx = iel;
        break;
    end
end

