# Descomposicion de Cholesky para matrices tridiagonales.
% Entrada: Matriz A
% Salida: Descomposicion  con L triangular inferior, T triangular superior
# Este codigo es de orden n.

function [L , T] = CholeskyTriDiag (A)  
n=size(A);
L=zeros(n,n);  
T=zeros(n,n);
  for i = 1:n
    if (i == 1)
      L(i,i) = sqrt( A(i,i) );
      T(i,i) = sqrt( A(i,i) );
    elseif  (i <n)
      L(i,i-1) =  A(i,i-1) / L(i-1,i-1);
      T(i-1,i) =  A(i,i-1) / L(i-1,i-1);
      
      L(i,i) =  sqrt(A(i,i) - (L(i,i-1)^2)  );
      T(i,i) =  sqrt(A(i,i)- (L(i,i-1)^2)  );
    else
      L(i,i-1) =  A(i,i-1) / L(i-1,i-1);
      T(i-1,i) =  A(i,i-1) / L(i-1,i-1);
      
      L(i,i) =  sqrt(A(i,i)- (L(i,i-1)^2)  );
      T(i,i) =  sqrt(A(i,i)- (L(i,i-1)^2)  );
    
  
    endif
  endfor  
L;
T;
endfunction


