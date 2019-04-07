function f=orbitfun(x,mu);

% Function Routine for Orbital Equations Runge Kutta Function

f=zeros(6,1);
r=sqrt(x(1)^2+x(2)^2+x(3)^2);
r3=r^3;
f(1)=x(4);
f(2)=x(5);
f(3)=x(6);
f(4)=-mu/r3*x(1); % force per unit mass or acceleration in x direction
f(5)=-mu/r3*x(2); % force per unit mass or acceleration in y direction
f(6)=-mu/r3*x(3); % force per unit mass or acceleration in z direction