%input
clc;
clear;
n=9; % number of members
I=[1 1 1 1 1 1 1 1 1]; %Moment of inertia in m4
L=[3 3 3 3 3 3 4.242 3 4.242]; %length in m
A=[0.005 0.004 0.004 0.005 0.004 0.004 0.006 0.005 0.006]; %Area in m2
theta=[90 0 0 -90 0 0 -45 -90 -135];  %angle in degrees
uu=8; %Number of unrestrained degrees of freedom
ur=4; %Number of restrained degrees of freedom 
uul=[1 2 3 4];
url=[5 6 7 8 9 10 11 12];
l1=[10 2 9 1];
l2=[2 4 1 3];
l3=[4 6 3 5];
l4=[6 12 5 11];
l5=[8 12 7 11];
l6=[10 8 9 7];
l7=[2 8 1 7];
l8=[4 8 3 7];
l9=[6 8 5 7];
l=[l1; l2; l3; l4; l5; l6; l7; l8; l9];
dof= uu + ur;
Ktotal=zeros (dof);
Tt1=zeros (4);
Tt2=zeros(4);
fem1=[0;0;0;0];
fem2=[0;0;0;0];
fem3=[0;0;0;0];
fem4=[0;0;0;0];
fem5=[0;0;0;0];
fem6=[0;0;0;0];
fem7=[0;0;0;0];
fem8=[0;0;0;0];
fem9=[0;0;0;0];
%% rotation coefficient for each member
rc1 =4.*I./L;
rc2 =2.*I./L;
rc3=A./L;
cx=cosd(theta);
cy=sind(theta);
%%stiffness matrix
for i=1:n
    Knew =zeros(dof);
    k1=[0; 0; 0; 0];
    k2=[0; 0; 0; 0];
    k3=[0;0;rc3(i);-rc3(i)];
    k4=-k3;
    K=[k1 k2 k3 k4];
    fprintf('Member Number=');
    disp(i);
    fprintf('Local stiffness matrix of member, [k]=\n');
    disp(K);
    T1=[cx(i);0;cy(i); 0];
    T2=[0;cx(i);0;cy(i)];
    T3=[-cy(i);0;cx(i);0];
    T4=[0;-cy(i);0;cx(i)];
    T=[T1 T2 T3 T4];
    fprintf('Transformation of member, [T]=n\n');
    disp(T);
    Ttr=T';
    fprintf('Transformation of Transpose, [T]=\n');
    disp(Ttr);
    Kg=Ttr*K*T;
    fprintf('Global Matrix, [K global]=\n');
    disp(Kg);
    for p=1:4
        for q=1:4
            Knew((l(i,p)),(l(i,q)))=Kg(p,q);
        end
    end
    Ktotal=Ktotal + Knew;
        if i==1
            Tr1=T;
            Kg1=Kg;
            fembar1=Tt1'*fem1;
        elseif i==2
            Tt2=T;
            Kg2=Kg;
            fembar2=Tt2'*fem2;
        elseif i==3
            Tt3=T;
            Kg3=Kg;
            fembar3=Tt3'*fem3;
        elseif i==4
            Tt4=T;
            Kg4=Kg;
            fembar4=Tt4'*fem4;
        elseif i==5
            Tt5=T;
            Kg5=Kg;
            fembar5=Tt5'*fem5;
          elseif i==6
            Tt6=T;
            Kg6=Kg;
            fembar6=Tt6'*fem6;
        elseif i==7
            Tt7=T;
            Kg7=Kg;
            fembar7=Tt7'*fem7;
        elseif i==8
            Tt8=T;
            Kg8=Kg;
            fembar8=Tt8'*fem8;
        elseif i==9
            Tt9=T;
            Kg9=Kg;
            fembar9=Tt9'*fem9;
        end
end
fprintf('Stiffness Matrix of complete structure, Ktotal =\n');
disp(Ktotal);
Kunr=zeros(4);
for x=1:uu
    for y=1:uu
        Kunr(x,y)=Ktotal (x,y);
    end
end
fprintf('unrestraied Stiffness of sub-matrix, [Kuu]=\n');
disp(Kunr);
KuuInv=inv(Kunr);
fprintf('Inverse of unrestrained sub-matrix, Kuu inverse=\n');
disp(KuuInv);
%% Creation of Joint vector
jl=[20;0;0;-70;0;0;0;0;0;0;0;0];
jlu=[20;0;0;-70;0;0;0;0];
delu=KuuInv*jlu;
fprintf('Joint Load vector, Jl =\n');
disp(jl);
fprintf('displacements are Delu =\n');
disp(delu);
delr=[0;0;0;0];
del=zeros(dof,1);
del=[delu;delr];
deli=zeros(4,1);
for i=1:n
    for p=1:4
        deli(p,1)=del((l(i,p)),1);
    end
    if i==1 
        delbar1=deli;
        mbar1=(Kg1*delbar1)+fembar1;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar1);
        fprintf('Global End moment matrix=\n');
        disp(mbar1);
    elseif i==2
        delbar2=deli;
        mbar2=(Kg2*delbar2)+fembar2;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar2);
        fprintf('Global End moment matrix=\n');
        disp(mbar2);
    elseif i==3
        delbar3=deli;
        mbar3=(Kg3*delbar3)+fembar3;
        fprintf('Member Number=');
        disp(i); 
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar3);
        fprintf('Global End moment matrix=\n');
        disp(mbar3);
   elseif i==4
        delbar4=deli;
        mbar4=(Kg4*delbar4)+fembar4;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar4);
        fprintf('Global End moment matrix=\n');
        disp(mbar4);
   elseif i==5
        delbar5=deli;
        mbar5=(Kg5*delbar5)+fembar5;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar5);
        fprintf('Global End moment matrix=\n');
        disp(mbar5);
     elseif i==6
        delbar6=deli;
        mbar6=(Kg6*delbar6)+fembar6;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar6);
        fprintf('Global End moment matrix=\n');
        disp(mbar6);
    elseif i==7
        delbar7=deli;
        mbar7=(Kg7*delbar7)+fembar7;
        fprintf('Member Number='); 
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar7);
        fprintf('Global End moment matrix=\n');
        disp(mbar7);
   elseif i==8
        delbar8=deli;
        mbar8=(Kg8*delbar8)+fembar8;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar8);
        fprintf('Global End moment matrix=\n');
        disp(mbar8);
   elseif i==9
        delbar9=deli;
        mbar9=(Kg9*delbar9)+fembar9;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacemant matrix [DeltaBar]=\n');
        disp(delbar9);
        fprintf('Global End moment matrix=\n');
        disp(mbar9);
    end
end
        
        


