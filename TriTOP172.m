function TriTOP172(nx,ny,ER,tau,Vmax,maxedge,minedge,E,nu)
%% PARAMETER CHECK
if (tau/E>1e-12) || (ER>0.1) error('Improper input'); end
%% INITIALIZATION
BDY = [-.5*(nx) -.5*(ny); .5*(nx), .5*(ny)]/100;
[xn,yn] = meshgrid(BDY(1,1):(BDY(2,1)-BDY(1,1))/(nx):BDY(2,1), BDY(1,2):(BDY(2,2)-BDY(1,2))/(ny):BDY(2,2));
dN = sin(xn/BDY(2,1)*6*pi).*cos(yn/BDY(2,1)*6*pi)+0.5;
%% MAIN LOOP
for iterNum = 1:200
%% GENERATE MESH AND NONLINEAR DIFFUSION TERM
[p,t,t1,t2,Ve] = GenerateMesh(xn,yn,dN,maxedge,minedge,BDY,80);
[nonlidiff] = NonlinearDiffsion(p,t,t1,Ve);
%% FINITE ELEMENT ANALYSIS AND SENSITIVITY
[Ce,J,vol] = FEA(t,t1,t2,p,Ve,BDY,E,nu);
V = max(Vmax,vol-ER);
Cn = E2N(t,p,Ce,Ve)+tau*nonlidiff';
if iterNum > 1  Cnlast = griddata(pold(:,1),pold(:,2),Cnold,p(:,1),p(:,2),'cubic'); end
Cnold = Cn; pold = p;
if iterNum >1   Cn = 0.5*(Cn+Cnlast); end
%% PRINT RESULTS
clf; patch('Faces',t1,'Vertices',p,'EdgeColor',[250 250 250]/255,'FaceColor',[192 192 192]/255); hold on;
patch('Faces',t2,'Vertices',p,'EdgeColor','k','FaceColor',[255 255/2 0]/255); axis off equal; pause(1e-3);
disp(['It.: ' num2str(iterNum) '  Obj.: ' sprintf('%6.4f',J) '  Vol.: ' sprintf('%6.4f' ,vol) ]);
%% CONVERGENCE CHECK
if iterNum>10 && abs(vol-Vmax)<0.005 && all(abs(Obj(end)-Obj(end-5:end-1)) <0.005*abs(Obj(end)))
return;   end
%% DESIGN UPDATE USING BESO METHOD
l1 = 0; l2 = 1e-3;
while ((l2-l1)/l2 > 1e-3)
CnNew = Cn-(l1+l2)/2;
dN = imgaussfilt(griddata(p(:,1),p(:,2),CnNew,xn,yn,'cubic'),0.8);
if isempty(contour(xn,yn,dN,[0 0]))==0
[~,~,t1,~,Ve] = GenerateMesh(xn,yn,dN,maxedge,minedge,BDY,1);
Vtemp = 1-sum(Ve(1:length(t1)))/sum(Ve);
else  Vtemp = dot(max(0,sign(mean(CnNew(t),2))),Ve)/sum(Ve);
end
if Vtemp > V  l1 = (l1+l2)/2; else  l2 = (l1+l2)/2;  end
end
Obj(iterNum) = J;
end
%% MESH GENERATOR
function [p,t,t1,t2,Ve] = GenerateMesh(xn,yn,dN,maxedge,minedge, BDY,maxiter)
x = xn; x(2:2:end,:) = x(2:2:end,:)+0.005; pi = [x(:),yn(:)];
C = contour(xn,yn,dN,[0 0]); C1 = [];
dump = find((C(1,:)==0) & (C(2,:)>=2));
for i = 1:length(dump)-1
C1 = [C1 (RedisPoints (C(:,dump(i)+1:dump(i+1)-1),0.004,0.01))'];
end
C = [C1 (RedisPoints (C(:,dump(end)+1:end),0.004,0.01))'];
d = zeros(size(pi,1),1);
for i = 1:size(pi,1)
d(i) = sqrt(min((pi(i,1)-C(1,:)).^2+(pi(i,2)-C(2,:)).^2));
end
r0 = 1./min(max(minedge,d),maxedge).^2;
pfix=[C'; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),0];
p = [pfix; pi(r0./max(r0)>0.5,:)];
p1 = 0; warning off;
for i = 1:maxiter
if max(sum((p-p1).^2,2))>1e-6
t = delaunayn(p);
edges = unique(sort([t(:,[1,2]);t(:,[1,3]); t(:,[2,3])],2),'rows');
p1 = p;
end 
midpoint = (p(edges(:,1),:)+p(edges(:,2),:))/2;
d = zeros(size(midpoint,1),1);
for j = 1:size(midpoint,1)
d(j) = sqrt(min((midpoint(j,1)-C(1,:)).^2+(midpoint(j,2)-C(2,:)).^2));
end
L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2));
L1 = min(max(minedge,d),maxedge);
L0 = 1.2*L1*sqrt(sum(L.^2)/sum(L1.^2));
Fe = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:));
Fp = full(sparse(edges(:,[1,1,2,2]),ones(size(d)) *[1,2,1,2],[Fe,-Fe],size(p,1),2));
Fp(1:size(pfix,1),:) = 0;
p = p+0.2*Fp;
p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1)));
p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2)));
end
[p, t] = UniqueNode(p,t);
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');
t1 = t(dE<0,:); t2 = t(dE>=0,:); t = [t1;t2];
for kk = 1:length(t)
Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
end
%% REMOVE REPEATED NODES
function [P,T] = UniqueNode(p,t)
for kk = 1:length(t)
Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
end
t((Ve==0),:) = [];
P=unique(p(unique(sort(t(:))),:),'rows');
for i = 1:length(t)
for j = 1:3  T(i,j) = find(P(:,1)==p(t(i,j),1)&P(:,2)==p(t(i,j),2)); end
end
%% CONTOUR POINTS DISTANCE ADJUST
function [C1] = RedisPoints(C,d1,d2)
C = C'; C1 = C;
CL = sqrt(sum(diff(C1,1,1).^2,2));
for i = 1:(size(C,1)-1)
if  CL(i) < d1  C1(i,:) = [0;0]; C1((i+1),:) = 0.5*(C(i,:)+C((i+1),:)); end
end
C1(all(C1==0,2),:) = [];
CL2 = sqrt(sum(diff(C1,1,1).^2,2));
Cmid = [];
for i = 1:(size(C1,1)-1)
if  CL2(i) > d2  Cmid = [Cmid; 0.5*(C1(i,:)+C1((i+1),:))]; end
end
if isempty(Cmid)==0
C1 = union(C1,Cmid,'rows');
end
%% CALCULATE NODE SENSITIVITY
function [dN] = E2N(t,p,x,Ve)
dN = zeros(length(p),1);
for i = 1:length(p)
[row,~] = find(t==i);
dN(i) = dot(Ve(row),x(row))/sum(Ve(row));
end
%% FINITE ELEMENT ANALYSIS
function [Ce,J1,vol] = FEA(t,t1,t2,p,Ve,BDY,E,nu)
NT = length(t);
KK = zeros(6,6*NT);
for i = length(t1)+1:NT
KK(:,6*i-5:6*i) = GetKe(p(t(i,:),1),p(t(i,:),2),E,nu);
end
elemDof = zeros(NT,6);
elemDof(:,[1 3 5]) = 2*t-1;
elemDof(:,[2 4 6]) = 2*t;
iK = reshape(kron(elemDof,ones(6,1))',36*NT,1);
jK = reshape(kron(elemDof,ones(1,6))',36*NT,1);
sK = reshape(KK,36*NT,1);
NK = sparse(iK,jK,sK,2*length(p),2*length(p));
NK = (NK+NK')/2;
fixedNodes = find(p(:,1)==BDY(1,1));
forceNodes = find(p(:,1)==BDY(2,1) & p(:,2)==0);
fixedDof = [2*fixedNodes-1; 2*fixedNodes];
SolidDof = [2*unique(sort(t2(:)))-1; 2*unique(sort(t2(:)))];
freeDofs = setdiff(SolidDof,fixedDof);
U = zeros(2*length(p),1);
F = sparse(2*forceNodes,1,-100,2*length(p),1);
U(freeDofs,:) = NK(freeDofs,freeDofs) \ F(freeDofs,1);
for i = length(t1)+1:NT
Ce(i) = 0.5 .* sum((U(elemDof(i,:))'*KK(:,6*i-5:6*i)).*U(elemDof(i,:))',2);
end
J1 = sum(Ce);
vol = 1-sum(Ve(1:length(t1)))/sum(Ve);
%% ELEMENT STIFFNESS MATRIX
function [Ke] = GetKe(X,Y,E0,nu)
D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2];
J = [X(1)-X(3) Y(1)-Y(3);X(2)-X(3) Y(2)-Y(3)];
Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
Ke = 1/2*det(J)*Be'*D*Be;
%% NONLINEAR DIFFUSION TERM DERIVATION
function [nonlidiff] = NonlinearDiffsion(p,t,t1,Ve)
f = E2N(t,p,sparse(length(t1)+1:length(t),1,1,length(t),1),Ve);
for A0 = 1:length(p)
Lf = []; [row , ~] = find(t==A0); ts = t(row,:);
PT = setdiff(ts,A0);
for i = 1:length(PT)
angl = [];
[secp , ~] = find(ts==PT(i));
for k = 1:length(secp)
A1 = PT(i); A2 = setdiff(ts(secp(k),:),[A0 A1]);
angl(k) = atan2(2*Ve(secp(k)),dot(p(A0,:)-p(A2,:),p(A1,:)-p(A2,:)));
end
sigma = .8; g = 1/(1+((f(A1)-f(A0))^2/sigma^2));
if length(angl)==2
Lf(i) = (cot(angl(1))+cot(angl(2)))*(f(A1)-f(A0))*g;
else  Lf(i) = (2*cot(angl(1)))*(f(A1)-f(A0))*g; end
end
nonlidiff(A0) = 1/(2*sum(Ve(row)))*sum(Lf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by ZC. Zhuang, YM. Xie, Q. Li and SW. Zhou %
% Centre for Innovative Structures and Materials, RMIT University         %
% Please send your comments to: zhuanginmelbourne@gmail.com               %
%                                                                         %
% The program is proposed for educational purposes and is introduced in   %
% the paper - A 172-line Matlab code for structural topology              %
% optimization in the body-fitted mesh, SMO, 2022                         %
%                                                                         %
% Disclaimer:                                                             %
% The authors reserve all rights but do not guarantee that the code is    %
% free from errors. Furthermore, we shall not be liable in any event.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CANTILEVER BEAM OPTIMIZATION
% TriTOP172(80,50,0.05,2e-8,0.5,0.2,0.025,1e5,0.3)

%% MBB BEAM OPTIMIZATION
% TriTOP172(120,40,0.05,2e-8,0.5,0.2,0.025,1e5,0.3)
%% INITIAL PATTERN (REPLACE LINE 7)
% dN = sin(xn/BDY(2,1)*6*pi).*cos(yn/BDY(2,1)*6*pi)+0.5; dN (1:5,:)= 1;
%% LOAD AND BOUNDARY CONDITIONS (REPLACE LINES 134-140)
% fixedNodes1 = find(p(:,1)==BDY(1,1));
% fixedDof1 = 2*fixedNodes1-1;
% fixedNodes2 = find(p(:,1)==BDY(2,1) & p(:,2)==BDY(1,2));
% fixedDof2 = 2*fixedNodes2;
% fixedDof = [fixedDof1;fixedDof2];
% forceNodes = find(p(:,1)==BDY(1,1) & p(:,2)==BDY(1,2));
% SolidDof = [2*unique(sort(t2(:)))-1; 2*unique(sort(t2(:)))];
% freeDofs = setdiff(SolidDof,fixedDof);
% U = zeros(2*length(p),1);
% F = sparse(2*forceNodes,1,-50,2*length(p),1);

%% L-BRACKET OPTIMIZATION
% TriTOP172(100,100,0.05,2e-13,0.4,0.2,0.025,2,0.3)
%% INITIAL PATTERN (REPLACE LINE 7)
% dN = sin(xn/BDY(2,1)*4*pi).*cos(yn/BDY(2,1)*4*pi)+0.5;
%% PASSIVE ELEMENTS DEFINITION (ADD AFTER LINE 43)
% px = 0.2 * BDY(1,1); py = 0.2 * BDY(1,2);
% [xp,yp] = meshgrid(px,py:0.01:BDY(2,2));
% [xp2,yp2] = meshgrid(px:0.01:BDY(2,1),py);
% Forcepts = [BDY(2,1) py; BDY(2,1)-0.01 py ;BDY(2,1)-0.02 py ;BDY(2,1)-0.03 py ;BDY(2,1)-0.04 py; BDY(2,1)-0.05 py];
% P = [xp(:),yp(:); xp2(:),yp2(:)];
%% PASSIVE ELEMENTS DEFINITION (REPLACE LINE 55)
% pfix=[C'; P; Forcepts; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2)];
%% PASSIVE ELEMENTS DEFINITION (ADD AFTER LINE 69)
% L = min(max(0.001,L),0.15);
%% PASSIVE ELEMENTS DEFINITION (REPLACE LINE 82)
% tnp = t(~(pmid(:,1)>px & pmid(:,2)>py),:);
% tp = t((pmid(:,1)>px & pmid(:,2)>py),:);
% dEnp = dE(~(pmid(:,1)>px & pmid(:,2)>py),:);
% tv = tnp(dEnp<0,:);
% t1 = [tp; tnp(dEnp<0,:)];
% t2 = tnp(dEnp>0,:);
% t = [t1;t2];
%% LOAD AND BOUNDARY CONDITIONS (REPLACE LINES 123-124)
% for i = 1:NT
%   if i<=length(t1) x=1e-5; else x=1; end
%   KK(:,6*i-5:6*i) = x*GetKe(p(t(i,:),1),p(t(i,:),2),E,nu);
%% LOAD AND BOUNDARY CONDITIONS (REPLACE LINES 134-140)
% fixedNodes=find(p(:,2)==BDY(2,1));
% fixedDof = [2*fixedNodes-1; 2*fixedNodes];
% SolidDof  = 1:2*length(p);
% forceNodes = find(p(:,1)<=BDY(2,1) & p(:,1)>=BDY(2,1)-0.05 & p(:,2)==0.2*BDY(1,2));
% freeDofs = setdiff(SolidDof,fixedDof);
% U = zeros(2*length(p),1);
% F = sparse(2*forceNodes,1,-0.1,2*length(p),1);