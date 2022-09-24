% Q4 element patch test
X=[3 5 5 3]; % global x coords
Y=[2 2 4 4]; % global y coords
conn=[1 2 3 4]; % element connectivity
E=30e06; % E
nu=0.25; % nu
delta=1e-06; % delta_0
% global stiffness and global force initialised
Kg=zeros(8);
Fg=zeros(8,1);
for ielem=1:size(conn,1) % loop over elements starts
    n1=conn(ielem,1);
    n2=conn(ielem,2);
    n3=conn(ielem,3);
    n4=conn(ielem,4);
    dest(ielem,:)=[2*n1-1 2*n1 2*n2-1 2*n2 2*n3-1 2*n3 2*n4-1 2*n4];

Ke=zeros(8); %local stiffness matrix
Fe=zeros(8,1); % local force vector
weights=[1 1 1 1];
r=[-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
s=[-1/sqrt(3) 1/sqrt(3) -1/sqrt(3) 1/sqrt(3)];
for lint = 1:length(weights) % loop over gauss points start, weights are the GP weights you have chosen
N1 =((r(lint) + 1)*(s(lint) + 1))/4;
N2 =-((r(lint) - 1)*(s(lint) + 1))/4;
N3 =((r(lint) - 1)*(s(lint) - 1))/4;
N4 =-((r(lint) + 1)*(s(lint) - 1))/4;

N1r =(s(lint)/4) + (1/4);
N1s =(r(lint)/4) + (1/4);
N2r =-(s(lint)/4) - (1/4);
N2s =(1/4) - (r(lint)/4);
N3r =s(lint)/4 - 1/4;
N3s =r(lint)/4 - 1/4;
N4r =1/4 - s(lint)/4;
N4s =-(r(lint)/4) + (-1/4);

J=[N1r*X(1)+N2r*X(2)+N3r*X(3)+N4r*X(4) N1r*Y(1)+N2r*Y(2)+N3r*Y(3)+N4r*Y(4);N1s*X(1)+N2s*X(2)+N3s*X(3)+N4s*X(4),N1s*Y(1)+N2s*Y(2)+N3s*Y(3)+N4s*Y(4)];

L=inv(J);

B=[L(1,1)*N1r+L(1,2)*N1s 0 L(1,1)*N2r+L(1,2)*N2s 0 L(1,1)*N3r+L(1,2)*N3s 0 L(1,1)*N4r+L(1,2)*N4s 0;
    0 L(2,1)*N1r+L(2,2)*N1s 0 L(2,1)*N2r+L(2,2)*N2s 0 L(2,1)*N3r+L(2,2)*N3s 0 L(2,1)*N4r+L(2,2)*N4s;
    L(2,1)*N1r+L(2,2)*N1s L(1,1)*N1r+L(1,2)*N1s L(2,1)*N2r+L(2,2)*N2s  L(1,1)*N2r+L(1,2)*N2s L(2,1)*N3r+L(2,2)*N3s L(1,1)*N3r+L(1,2)*N3s L(2,1)*N4r+L(2,2)*N4s L(1,1)*N4r+L(1,2)*N4s ];


D = (E/(1-nu*nu))*[1 nu 0 ; nu 1 0 ;0 0 (1-nu)/2] ;  %Elasticity Matrix

Ke=B'*D*B*det(J)*weights(lint)+Ke;
end  % gauss loop ends
% assemble
for idof=1:8
for jdof=1:8
Kg(dest(ielem,idof),dest(ielem,jdof)) = Kg(dest(ielem,idof),dest(ielem,jdof)) + Ke(idof,jdof);
end
end
end % element loop ends

% insert boundary conditions
fixed=[1 2 4 6 7 8];
Fg = Fg - (delta)*Kg(:,6)-(delta)*Kg(:,8);
Fg(fixed)=0;
Fg(6)=delta;Fg(8)=delta;
Kg(fixed,:)=0; Kg(:,fixed)=0;
for jj=1:length(fixed)
    Kg(fixed(jj),fixed(jj))=1;
end
% solve for U, the nodal displacements
U=inv(Kg)*Fg
