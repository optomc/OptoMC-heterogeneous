function [adjElem adjElemNode] = getAdjElemNode(elem)
%% initiation
% This function automatically calculate adjacent elements of each element.
% It takes quite long time and user of our mesh-based MC needn't run this code again
% since we provide mesh binary adjacent elements were precomputed (ratMesh_adjElem.mat).
%
% Usage: 
%  load rat_rawMesh.mat
%  getAdjElemNode(elem)
%
elem=elem(:,1:4);
f1=[1 2 3];
f2=[1 2 4];
f3=[2 3 4];
f4=[1 3 4];
adjElem = ones(length(elem),4);
adjElemNode = ones(length(elem),16);

for i = 1:length(elem)
    % getting surface
    s=[elem(i,f1);elem(i,f2);elem(i,f3);elem(i,f4)];
    % find common surface
    tmp0 = ismember(elem,s);
    tmp1=sum(ismember(elem,s),2);
    tmp2=find(tmp1==3);

    mulElem = elem(tmp2,:) .* tmp0(tmp2,:);
    xl=length(mulElem(:,1));
    yl=length(mulElem);
    stack = zeros(1,xl*yl-xl);
    k = 1;
    for x = 1:xl
        for y = 1:yl
            if mulElem(x,y) ~= 0
                stack(k) = mulElem(x,y);
                k = k+1;
            end
        end
    end
    stack2 = reshape(stack, yl-1,xl)';
    adjElemNodeTrans = [tmp2, stack2]';
    adjElemNodeList  = adjElemNodeTrans(:);

    if(~isempty(tmp2))
        adjElem(i,1:length(tmp2))=tmp2;
        adjElemNode(i,1:length(adjElemNodeList))=adjElemNodeList;
    end
end
end