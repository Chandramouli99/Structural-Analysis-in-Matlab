%input
clc;
clear;
n=3; % number of members
I=[1 1 1]; %Moment of inertia in m4
L=[4.242 3 3]; %length in m
A=[5 5 5]; %Area in m2
theta=[45 -90 0];  %angle in degrees
uu=2; %Number of unrestrained degrees of freedom
ur=4; %Number of restrained degrees of freedom 
uul=[1 2];
url=[3 4 5 6];
l1=[3 1 4 2];
l2=[1 5 2 6];
l3=[3 5 4 6 ];
l=[l1; l2; l3];
dof= uu + ur;
Ktotal=zeros (dof);
Tt1=zeros (4);
Tt2=zeros(4);
fem1=[0;0;0;0];
fem2=[0;0;0;0];
fem3=[0;0;0;0];
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
        end
end
fprintf('Stiffness Matrix of complete structure, Ktotal =\n');
disp(Ktotal);
Kunr=zeros(2);
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
jl=[-10;40;0;0;0;0];
jlu=[-10;40];
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
    end
end
        
        


