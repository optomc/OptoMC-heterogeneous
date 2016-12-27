function [faceElemId]=myMeshCutSourceElem(elem,node,plane)
%
% compute tetrahedral mesh cross-section 
% this code is based on qmeshcut which is function of iso2mesh.

[a,b,c,d]=getplanefrom3pt(plane); % function of iso2mesh toolbox
% elem = elem(:,1:4); node = node(:,1:3);
% compute which side of the plane for all nodes in the mesh

co=repmat([a b c],size(node,1),1);
dist=sum( (co.*node)' )+d;
asign=dist;
asign(asign>=0)=1;
asign(asign<0)=-1;

% get all the edges of the mesh
esize=size(elem,2);
edges=[elem(:,[1,2]);elem(:,[1,3]);elem(:,[1,4]); elem(:,[2,3]);elem(:,[2,4]);elem(:,[3,4])];

% find all edges with two ends at the both sides of the plane
edgemask=sum(asign(edges(:,1:2)),2);
cutedges=find(edgemask==0);

% organize all cross-cuts into patch facedata format
emap=zeros(size(edges,1),1);
emap(cutedges)=1:length(cutedges);
emap=reshape(emap,[size(elem,1),esize*(esize-1)/2]); % C^n_2

etag=sum(emap>0,2);

tricut=find(etag==3);
quadcut=find(etag==4);

faceElemId=[tricut; quadcut];