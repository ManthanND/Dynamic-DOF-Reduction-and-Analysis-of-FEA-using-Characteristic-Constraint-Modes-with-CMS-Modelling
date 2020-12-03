%% CMS Modelling using Craig Bampton method for Sub-structure 1
K3=full(K3);
% Partitioning K1 matrix into Interior and Boundary nodes
A22=K3([421:441],[421:441]);
A42=K3([1:420,442:882],[421:441]);
A52=K3([1:420,442:882],[1:420,442:882]);
% Forming Mii(Mass submatrix with interior modes
MA52=M3([1:420,442:882],[1:420,442:882]);
psic=-inv(A52)*A42;
% Primary modal analysis of interior nodes
[X,e2,s2]=eig(A52,MA52);
%Creating the transformation matrix Rcb
RCB12=horzcat(X(:,[1:10]),psic);
I=eye(21,21);
Z=zeros([21,10]);
IZ=horzcat(Z,I);
RCB=vertcat(RCB12,IZ);
% Rearranging the Mass and Stiffness matrices according to the interior and boundary nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%
A=K3([421:441],[421:441]);
A4=K3([1:420,442:882],[421:441]);
A5=K3([1:420,442:882],[1:420,442:882]);
A6=K3([421:441],[1:420,442:882]);
H=horzcat(A5,A4);
H1=horzcat(A6,A);
Kf=vertcat(H,H1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
KA=M3([421:441],[421:441]);
KA4=M3([1:420,442:882],[421:441]);
KA5=M3([1:420,442:882],[1:420,442:882]);
KA6=M3([421:441],[1:420,442:882]);
KH=horzcat(A5,A4);
KH1=horzcat(A6,A);
Mf=vertcat(KH,KH1);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Transformation matrix for Mass and Stiffness reduction
MT=transpose(RCB)*Mf*RCB;
KT=transpose(RCB)*Kf*RCB;
[X12,e12,s12]=eigs(KT,MT,6,'sm')

%%%Extracting individual submatrices from reduced mass and stiffness matrix%%%%%%%%%
X=X(:,[1:10])
Mnn=eye(10,10);
Mnb=(transpose(X)*KA5*psic)+(transpose(X)*KA4);
Mbn=transpose(Mnb);
Mbb=transpose(psic)*KA5*psic+transpose(psic)*KA4+KA6*psic+KA;
Mr1=horzcat(Mnn,Mnb);
Mr2=horzcat(Mbn,Mbb);
Mr=vertcat(Mr1,Mr2);
MnnC=Mnn
MnbC=Mnb
MbnC=Mbn
MbbC=Mbb

Knn=transpose(X)*A5*X;
Knb=(transpose(X)*A5*psic)+(transpose(X)*A4);
Kbn=transpose(Knb);
Kbb=transpose(psic)*A5*psic+transpose(psic)*A4+A6*psic+A;
KnnC=Knn
KnbC=Knb
KbnC=Kbn
KbbC=Kbb

Kr1=horzcat(Knn,Knb);
Kr2=horzcat(Kbn,Kbb);
Kr=vertcat(Kr1,Kr2);