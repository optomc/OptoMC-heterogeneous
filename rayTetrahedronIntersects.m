function [minDist, minNodes, minCos, minNormal] = rayTetrahedronIntersects(elemIdx, elem, node, photonPos, photonDir)  
    % compute ray/tetrahedron intersection test
    minDist = 1e10;    
    minNodes = [-1 -1 -1];
    
    % find distance to collision point, t      
    nodeIdx = [1 2 3 4 1 2];
    
    for i = 1:4             
        nodeIdx1 = elem(elemIdx,nodeIdx(i));
        nodeIdx2 = elem(elemIdx,nodeIdx(i+1));
        nodeIdx3 = elem(elemIdx,nodeIdx(i+2));
        
        n1 = node(nodeIdx1, 1:3);
        n2 = node(nodeIdx2, 1:3);
        n3 = node(nodeIdx3, 1:3);

        faceNormal = cross(n2-n1, n3-n1);
        faceNormal = faceNormal/norm(faceNormal);
              
        t = dot(faceNormal, (n1 - photonPos)) / dot(faceNormal,photonDir);          % (n1 - photonPos) is vector -P in the paper
        if (t < minDist && t > 1e-12) 
%             % make sure the normal vector points to inside of the tetrahedron.
            tmp = -(faceNormal(1)*photonPos(1) + faceNormal(2)*photonPos(2) + faceNormal(3)*photonPos(3) - dot(faceNormal,n1));
            if tmp > 0             
                minNormal = -faceNormal;    % points to outside, therefore change its direction.
            else
                minNormal = faceNormal;
            end                     
            minCos = dot(minNormal,photonDir);
            minDist  = t;   
            minNodes = [nodeIdx1 nodeIdx2 nodeIdx3];
        end
    end 

    if ~exist('minCos','var')
        minCos = -1; minDist = -1; minNormal = [-1 -1 -1]; 
    end
end