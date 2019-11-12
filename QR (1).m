function [Q, R] = gsC(A)
clc;
[m, n] = size(A);
R = zeros(m, n);
V = A;
Q=zeros(m, n);

k=1; % inicializo eta variable, sólo la uso para matrices de rango deficiente

for i =1:n
    R(i,i)= norm(V(:,i));  
   
    if abs(R(i,i)) > 1.0e-10
        Q(:,i)= V(:,i)/R(i,i);
    else  % salió cero! La matriz A es de rango deficiente
        r= zeros(1, m); % este va a ser el vector que reemplazará a la fila problemática de A
              
        Q(:,i)= r;
        Q(k, i)=1;
        k=k+1;
        Q2 = Q; % Q2 tendrá los valores de Q, en Q irá la nueva base ortonormal

        % ortogonalizo otra vez aplicando Gram-Schmidt a la matriz Q para
        % hallar un vector con el qué reemplazar a la fila problemática y
        % seguir teniendo una base ortonormal
        for j = 1 : i-1
            r(j) = (Q2(:,j)')*Q2(:,i);
            Q(:,i)=Q(:,i) - r(j)*Q2(:,j);
            Q(:,i)=Q(:,i)/norm(Q(:,i));
        end
    end
   
    for j=i+1:n
       R(i,j)= (Q(:,i)')*V(:,j);
       V(:,j)=V(:,j) - R(i,j)*Q(:,i);
    end
end

% hasta aquí Q y R son la descomposición reducida
% acá hallo la completa, Gram Schmidt otra vez
if m > n   
    I = eye(m, m - n);
    Q=[Q I];
    V=Q;
    r= zeros(1, m);
   
    k=0;
    for i=n+1:m
        for j = 1 : n + k
            r(j) = (V(:,j)')*V(:,i);
            Q(:,i)=Q(:,i) - r(j)*V(:,j);
            Q(:,i)=Q(:,i)/norm(Q(:,i));
        end
        k=k+1;
    end
end
end