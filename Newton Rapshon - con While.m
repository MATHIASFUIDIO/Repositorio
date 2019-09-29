##  Codigo de NR con un while.
# Variables de entrada: F, J,xo , MaxIteraciones y  MaxTolerancia
# 


function [x, Fx ,i] = NewtonRapshon (F,J,x , MaxIter , MaxTol )
  i=0;
  while i<MaxIter

    y=-J(x)\F(x);
    if sqrt(norm(y)/norm(x))< MaxTol
        break;
    endif
    x=x+y;
    i=i+1; 
    if i>=MaxIter
        disp('Se ha soprepasado el numero de iteracciones');
    endif
  endwhile
   x
  Fx = F(x)
   i
endfunction


## Ejemplo de NR

F=@(x)[x(1) - (1/3)*x(2) + 1 + (e^( x(1) + x(2) ) ) ; -(1/3)*x(1) + 3*x(2) - 2 + (e^( x(1) + x(2) ) )];
J=@(x)[1 + (e^(x(1) + x(2)) ) , -(1/3) + (e^(x(1) + x(2)) ) ; -(1/3) + (e^(x(1) + x(2)) ) ,3 + (e^(x(1) + x(2)) )];
x=[0;0]; %punto inicial

NewtonRapshon (F,J,x , 10 , 0.001 )