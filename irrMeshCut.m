function [cutPos,faceValue,faceData, faceElemId]=irrMeshCut(elem,node,medium,plane)

% this function is based on qmeshcut which is built-in function of iso2mesh
% in http://iso2mesh.sourceforge.net/cgi-bin/index.cgi

[a,b,c,d]=getplanefrom3pt(plane);
% elem = elem(:,1:4); node = node(:,1:3);
% compute which side of the plane for all nodes in the mesh

co=repmat([a b c],size(node,1),1);
dist=sum( (co.*node)' )+d;
asign=dist;
asign(find(asign>=0))=1;
asign(find(asign<0))=-1;

% get all the edges of the mesh
esize=size(elem,2);
edges=[elem(:,[1,2]);elem(:,[1,3]);elem(:,[1,4]); elem(:,[2,3]);elem(:,[2,4]);elem(:,[3,4])];

% find all edges with two ends at the both sides of the plane
edgemask=sum(asign(edges(:,1:2)),2);
cutedges=find(edgemask==0);

% calculate the distances of the two nodes, and use them as interpolation weight 
cutweight=dist(edges(cutedges,1:2));
totalweight=diff(cutweight');
cutweight=abs(cutweight./repmat(totalweight(:),1,2));
cutPos=node(edges(cutedges,1),:) .* repmat(cutweight(:,2),[1 3]) + node(edges(cutedges,2),:) .* repmat(cutweight(:,1),[1 3]);
   
% organize all cross-cuts into patch facedata format
emap=zeros(size(edges,1),1);
emap(cutedges)=1:length(cutedges);
emap=reshape(emap,[size(elem,1),esize*(esize-1)/2]); % C^n_2

etag=sum(emap>0,2);

tricut=find(etag==3);
quadcut=find(etag==4);

% fast way (vector-form) to get all triangles
tripatch=emap(tricut,:)';
tripatch=reshape(tripatch(find(tripatch)),[3,length(tricut)])';

quadpatch=emap(quadcut,:)';
quadpatch=reshape(quadpatch(find(quadpatch)),[4,length(quadpatch)])';

% combine the two sets to create the final facedata
% using the matching-tetrahedra algorithm as shown in 
% https://visualization.hpc.mil/wiki/Marching_Tetrahedra

faceData=[tripatch(:,[1 2 3 3]); quadpatch(:,[1 2 4 3])];
faceElemId=[tricut; quadcut];
faceValue=medium(faceElemId);