% Q4 element patch test
X=[3 5 5 3 4 5 4 3 4];
Y=[2 2 4 4 2 3 4 3 3.2];
conn=[1 5 9 8;5 2 6 9;9 6 3 7;8 9 7 4];
% get destination array
for ielem=1:size(conn,1) % loop over elements starts
    n1=conn(ielem,1);
    n2=conn(ielem,2);
    n3=conn(ielem,3);
    n4=conn(ielem,4);
    dest(ielem,:)=[2*n1-1 2*n1 2*n2-1 2*n2 2*n3-1 2*n3 2*n4-1 2*n4];
end
E=30e06;
nu=0.25
delta=1e-06;
fixed=[1 2 4 6 7 8 10 14 15];
fixed_vals=[0 0 0 delta 0 delta 0 delta 0];
nelem=4;
nnode=9;
ngauss=4; % use 2X2 gauss quadrature
%initialize K and F
Kg=zeros(2*nnode);
Fg=zeros(2*nnode,1);
Klocal=[]; %to store Ke for 4th element

% form element stiffness and assemble
for ielem=1:size(conn,1)
Ke=zeros(8); %local stiffness matrix
Fe=zeros(8,1); % local force vector

weights=[1 1 1 1];

r=[-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
s=[-1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];

for l= 1:length(weights) % loop over gauss points start, weights are the GP weights you have chosen
N1 =((r(l) + 1)*(s(l) + 1))/4;
N2 =-((r(l) - 1)*(s(l) + 1))/4;
N3 =((r(l) - 1)*(s(l) - 1))/4;
N4 =-((r(l) + 1)*(s(l) - 1))/4;

N1r =(s(l)/4) + (1/4);
N1s =(r(l)/4) + (1/4);
N2r =-(s(l)/4) - (1/4);
N2s =(1/4) - (r(l)/4);
N3r =s(l)/4 - 1/4;
N3s =r(l)/4 - 1/4;
N4r =1/4 - s(l)/4;
N4s =-(r(l)/4) + (-1/4);

J=[N1r*X(conn(ielem,1))+N2r*X(conn(ielem,2))+N3r*X(conn(ielem,3))+N4r*X(conn(ielem,4)) N1r*Y(conn(ielem,1))+N2r*Y(conn(ielem,2))+N3r*Y(conn(ielem,3))+N4r*Y(conn(ielem,4));N1s*X(conn(ielem,1))+N2s*X(conn(ielem,2))+N3s*X(conn(ielem,3))+N4s*X(conn(ielem,4)),N1s*Y(conn(ielem,1))+N2s*Y(conn(ielem,2))+N3s*Y(conn(ielem,3))+N4s*Y(conn(ielem,4))];

L=inv(J);

B=[L(1,1)*N1r+L(1,2)*N1s 0 L(1,1)*N2r+L(1,2)*N2s 0 L(1,1)*N3r+L(1,2)*N3s 0 L(1,1)*N4r+L(1,2)*N4s 0;
    0 L(2,1)*N1r+L(2,2)*N1s 0 L(2,1)*N2r+L(2,2)*N2s 0 L(2,1)*N3r+L(2,2)*N3s 0 L(2,1)*N4r+L(2,2)*N4s;
    L(2,1)*N1r+L(2,2)*N1s L(1,1)*N1r+L(1,2)*N1s L(2,1)*N2r+L(2,2)*N2s  L(1,1)*N2r+L(1,2)*N2s L(2,1)*N3r+L(2,2)*N3s L(1,1)*N3r+L(1,2)*N3s L(2,1)*N4r+L(2,2)*N4s L(1,1)*N4r+L(1,2)*N4s ];

D = (E/(1-nu*nu))*[1 nu 0 ; nu 1 0 ;0 0 (1-nu)/2];   %Elasticity Matrix

Ke=B'*D*B*det(J)*weights(l)+Ke;

end  % gauss loop ends
Klocal=[Klocal,Ke];
% assemble
for idof=1:8
for jdof=1:8
Kg(dest(ielem,idof),dest(ielem,jdof)) = Kg(dest(ielem,idof),dest(ielem,jdof)) + Ke(idof,jdof);
end
end
end % element loop ends
Kg
Ke=Klocal(1:8,25:32)   % store element stiffness for element 4

% insert boundary conditions
fixed=[1 2 4 6 7 8 10 14 15];
Fg=Fg - (delta)*Kg(:,6)-(delta)*Kg(:,8)-(delta)*Kg(:,14);
Fg(fixed)=0;
Fg(6)=delta;Fg(8)=delta;Fg(14)=delta;
Kg(fixed,:)=0; Kg(:,fixed)=0;
for jj=1:length(fixed)
    Kg(fixed(jj),fixed(jj))=1;
end
Kg
Fg
% compute displacement vector 
U = inv(Kg)*Fg
% find stresses
stress=zeros(3,1);
U3=[U(17);U(18);U(11);U(12);U(5);U(6);U(13);U(14)];
weights=[1 1 1 1];
r=[-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
s=[-1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
for l = 1:length(weights)
N1 =((r(l) + 1)*(s(l) + 1))/4;
N2 =-((r(l) - 1)*(s(l) + 1))/4;
N3 =((r(l) - 1)*(s(l) - 1))/4;
N4 =-((r(l) + 1)*(s(l) - 1))/4;

N1r =(s(l)/4) + (1/4);
N1s =(r(l)/4) + (1/4);
N2r =-(s(l)/4) - (1/4);
N2s =(1/4) - (r(l)/4);
N3r =s(l)/4 - 1/4;
N3s =r(l)/4 - 1/4;
N4r =1/4 - s(l)/4;
N4s =-(r(l)/4) + (-1/4);

J=[N1r*X(9)+N2r*X(6)+N3r*X(3)+N4r*X(7) N1r*Y(9)+N2r*Y(6)+N3r*Y(3)+N4r*Y(7);N1s*X(9)+N2s*X(6)+N3s*X(3)+N4s*X(7),N1s*Y(9)+N2s*Y(6)+N3s*Y(3)+N4s*Y(7)];

L=inv(J);

B=[L(1,1)*N1r+L(1,2)*N1s 0 L(1,1)*N2r+L(1,2)*N2s 0 L(1,1)*N3r+L(1,2)*N3s 0 L(1,1)*N4r+L(1,2)*N4s 0;
    0 L(2,1)*N1r+L(2,2)*N1s 0 L(2,1)*N2r+L(2,2)*N2s 0 L(2,1)*N3r+L(2,2)*N3s 0 L(2,1)*N4r+L(2,2)*N4s;
    L(2,1)*N1r+L(2,2)*N1s L(1,1)*N1r+L(1,2)*N1s L(2,1)*N2r+L(2,2)*N2s  L(1,1)*N2r+L(1,2)*N2s L(2,1)*N3r+L(2,2)*N3s L(1,1)*N3r+L(1,2)*N3s L(2,1)*N4r+L(2,2)*N4s L(1,1)*N4r+L(1,2)*N4s ];

D = (E/(1-nu*nu))*[1 nu 0 ; nu 1 0 ;0 0 (1-nu)/2]; %Elasticity Matrix

stress=D*B*U3+stress;
end
% average yy stress in element 3
stressyy3 = stress(2,1)/4