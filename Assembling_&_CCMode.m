% Assembling 2 3 4 substructures
KbbAB=blkdiag(KbbB,KbbC,KbbD)
MbbAB=blkdiag(MbbB,MbbC,MbbD)

KbnAB=blkdiag(KbnB,KbnC,KbnD)
MbnAB=blkdiag(MbnB,MbnC,MbnD)

KnbAB=blkdiag(KnbB,KnbC,KnbD)
MnbAB=blkdiag(MnbB,MnbC,MnbD)

MnnAB=eye(30,30)
KnnAB=eye(30,30)

% Assembling with Substructure 1
Mnassemb1=horzcat(MnnA,zeros(10,30),MnbA)
Mnassemb2=horzcat(zeros(30,10),MnnAB,MnbAB)
Mnassemb3=horzcat(MbnA,MbnAB,MbbA+MbbAB)
Mnassemb=vertcat(Mnassemb1,Mnassemb2,Mnassemb3)

Knassemb1=horzcat(KnnA,zeros(10,30),zeros(10,43))
Knassemb2=horzcat(zeros(30,10),KnnAB,zeros(30,43))
Knassemb3=horzcat(zeros(43,10),zeros(43,30),KbbA+KbbAB)
Knassemb=vertcat(Knassemb1,Knassemb2,Knassemb3)
% Characteristic modes from assembled matrix 
[Xassemb,es,ss]=eig(KbbA+KbbAB,MbbA+MbbAB)
%[XCB,eCB,sCB]=eigs(Knassemb,Mnassemb,50,'sm')
% Eigen modal analysis for Craig Bampton Reduction 
[XCB,eCB,sCB]=eig(Knassemb,Mnassemb)

Xassemb=Xassemb
%% Selecting i as number of CC modes to form CMS model reduction
i=10% No. of CC modes
MnbAT=transpose(Xassemb([1:10],[1:i]))*MnbA
MnbABT=transpose(Xassemb([1:30],[1:i]))*MnbAB
MnbAT=MnbA*Xassemb([1:43],[1:i])
MnbABT=MnbAB*Xassemb([1:43],[1:i])
MbnAT=transpose(Xassemb([1:43],[1:i]))*MbnA
MbnABT=transpose(Xassemb([1:43],[1:i]))*MbnAB
Mlast=transpose(Xassemb([1:43],[1:i]))*(MbbA+MbbAB)*Xassemb([1:43],[1:i])
Klast=transpose(Xassemb([1:43],[1:i]))*(KbbA+KbbAB)*Xassemb([1:43],[1:i])

Mrassemb1=horzcat(eye(10,10),zeros(10,30),MnbAT)
Mrassemb2=horzcat(zeros(30,10),eye(30,30),MnbABT)
Mrassemb3=horzcat(MbnAT,MbnABT,Mlast)
MR=vertcat(Mrassemb1,Mrassemb2,Mrassemb3)

Krassemb1=horzcat(KnnA,zeros(10,30),zeros(10,i))
Krassemb2=horzcat(zeros(30,10),KnnAB,zeros(30,i))
Krassemb3=horzcat(zeros(i,10),zeros(i,30),Klast)
KR=vertcat(Krassemb1,Krassemb2,Krassemb3)

%[XR,eR,sR]=eigs(KR,MR,20,'sm')
% Eien values for CC mode reduced assembly matrix
[eR]=eig(KR,MR)
