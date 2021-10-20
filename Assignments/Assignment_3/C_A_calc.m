syms Xu Yv Yr Nv Nr
syms u v r

A11 = [-Xu*u; -Yv*v; 0];
A12 = [0; -Yr*r; 0];
A21 = [0; 0; -Nv*v];
A22 = [0; 0; -Nr*r];
C_A = [zeros(3)          -Smtrx(A11+A12);
      -Smtrx(A11 + A12) -Smtrx(A21 + A22)]
  
C_A126 = [C_A(1:2, 1:2) C_A(1:2, 6);
          C_A(6, 1:2)   C_A(6, 6);]