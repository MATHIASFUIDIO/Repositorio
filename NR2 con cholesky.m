##
% Funcion NR2, realiza una modificacion de NR, con Cholesky para matrices tridiay y actualizacion de rango 1
% Entrada: F - Ecuacion, x- punto inicial, Q- Matriz tridiagonal , a - Vector de actualizacion
% MaxITER maxima cantidad de iteraciones, Tol - Tolerancia
% Salida: XSol - Vector Solucion, Fxsol - F evaluada en solucion ,IterTOT -  CantIteraciones, Tiempo - Tiempo de ejecucion
##
function [ XSol, Fxsol ,IterTOT ,Tiempo ] = NR2 (F,x,Q,a , MaxIter , MaxTol )
  tic;
  i=0;
  n= size(Q);
  uno = ones(1, n);

  #Calculo Descomposicion sobre Q
  [L , T] = CholeskyTriDiag (Q);
 

  while i<MaxIter
      #Aplico actualizacion de rango uno para obtener un B*B' = Jf
      [B]   = CholeskyModifyB(T,a);  
      Bt    = B';
      Faux  = -F(x);
      
      yaux  = linsolve(Bt , Faux, struct('LT', true));
      y     = linsolve(B  , yaux, struct('UT', true));
  
     if sqrt(norm(y)/norm(x))< MaxTol;
        break;
     endif
    x=x+y';
    i=i+1;
    a = uno'.*sqrt(e.^sum(x));
    if i>=MaxIter
        disp('Se ha soprepasado el numero de iteracciones');
    endif
  endwhile
  XSol      = x;
  Fxsol     = F(x);
  IterTOT   = i;
  Tiempo    = toc;
endfunction

