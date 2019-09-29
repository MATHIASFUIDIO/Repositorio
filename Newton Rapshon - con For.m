F=@(x)[x(1) - (1/3)*x(2) + 1 + (e^( x(1) + x(2) ) ) ; -(1/3)*x(1) + 3*x(2) - 2 + (e^( x(1) + x(2) ) )];
J=@(x)[1 + (e^(x(1) + x(2)) ) , -(1/3) + (e^(x(1) + x(2)) ) ; -(1/3) + (e^(x(1) + x(2)) ) ,3 + (e^(x(1) + x(2)) )];
x=[0;0]; %punto inicial
norma=[norm(x)]
i=0;
for i=1:10
    y=-J(x)\F(x);
    x=x+y;
    norma=[norma;norm(x)];
    i=i+1; 
end
x
norma
plot(norma)