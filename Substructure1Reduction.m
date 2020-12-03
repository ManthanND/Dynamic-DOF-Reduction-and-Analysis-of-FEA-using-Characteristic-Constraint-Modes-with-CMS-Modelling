%% CMS Modelling using Craig Bampton method for Sub-structure 1
K1=full(K1);
% Partitioning K1 matrix into Interior and Boundary nodes
A1=K1([629:639],[629:639]);
A2=K1([640:660],[640:660]);
A3=K1([661:671],[661:671]);
A=blkdiag(A1,A2,A3);
A=K1([629:671],[629:671]);
A4=K1([1:628,672:1342],[629:671]);
A5=K1([1:628,672:1342],[1:628,672:1342]);
% Forming Mii(Mass submatrix with interior modes
MA5=M1([1:628,672:1342],[1:628,672:1342]);
psic=-inv(A5)*A4;
% Primary modal analysis of interior nodes
[X,e,s]=eig(A5,MA5);
%Creating the transformation matrix Rcb
RCB1=horzcat(X(:,[1:10]),psic);
I=eye(43,43);
Z=zeros([43,10]);
IZ=horzcat(Z,I);
RCB=vertcat(RCB1,IZ);
% Rearranging the Mass and Stiffness matrices according to the interior and boundary nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%
A=K1([629:671],[629:671]);
A4=K1([1:628,672:1342],[629:671]);
A5=K1([1:628,672:1342],[1:628,672:1342]);
A6=K1([629:671],[1:628,672:1342]);
H=horzcat(A5,A4);
H1=horzcat(A6,A);
Kf=vertcat(H,H1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
KA=M1([629:671],[629:671]);%Mbb
KA4=M1([1:628,672:1342],[629:671]);%Mib
KA5=M1([1:628,672:1342],[1:628,672:1342]);%Mii
KA6=M1([629:671],[1:628,672:1342]);%Mbi
KH=horzcat(KA5,KA4);
KH1=horzcat(KA6,KA);
Mf=vertcat(KH,KH1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transformation matrix for Mass and Stiffness reduction
MT=transpose(RCB)*Mf*RCB;
KT=transpose(RCB)*Kf*RCB;
[X1,e1,s1]=eigs(KT,MT,6,'sm');

%%%%%%%%Extracting individual submatrices from reduced mass and stiffness matrix%%%%%%%%%
X=X(:,[1:10])
Mnn=eye(10,10);
Mnb=(transpose(X)*KA5*psic)+(transpose(X)*KA4);
Mbn=transpose(Mnb);
Mbb=transpose(psic)*KA5*psic+transpose(psic)*KA4+KA6*psic+KA;
Mr1=horzcat(Mnn,Mnb);
Mr2=horzcat(Mbn,Mbb);
Mr=vertcat(Mr1,Mr2);
MnnA=Mnn
MnbA=Mnb
MbnA=Mbn
MbbA=Mbb

Knn=transpose(X)*A5*X;
Knb=(transpose(X)*A5*psic)+(transpose(X)*A4);
Kbn=transpose(Knb);
Kbb=transpose(psic)*A5*psic+transpose(psic)*A4+A6*psic+A;
KnnA=Knn
KnbA=Knb
KbnA=Kbn
KbbA=Kbb

Kr1=horzcat(Knn,Knb);
Kr2=horzcat(Kbn,Kbb);
Kr=vertcat(Kr1,Kr2);