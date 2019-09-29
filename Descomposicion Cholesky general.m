#Funcion para hacer Descomposiocion de CHOLESKY.
% Entrada: A - MAtriz
# Devuelve la matriz B que es triangular inferior tal que B*B'= A
function  [B] = CholeskyGral (A)
  n=length(A);
  B=zeros(n,n);
  for j=1:n
    for i=j:n
      if i==j && j==1
        B(i,j)=sqrt(A(i,j));
      elseif j==1 && j!=i
        B(i,j)=A(i,j)/B(j,j);
      elseif j>1 && i==j
        B(i,j)=sqrt(A(i,j)-(B(i,[1:j-1])*B(i,[1:j-1])'));
      else
        B(i,j)=(A(i,j)-(B(j,[1:j-1])*B(i,[1:j-1])'))/B(j,j);
      endif
    endfor
  endfor
  B
endfunction

