%K = 0.0075;
%T = 169.5493;
%h = 0.1;

syms K T h
Ac = [0 1 0; 0 -1/T -K/T; 0 0 0];
Bc = [0; K/T; 0];
Ec = [0 0; 1 0; 0 1];
Cc = [1 0 0];
Dc = 0;

Ad = eye(3) + h*Ac;
Bd = h*Bc;
Ed = h*Ec;
Cd = Cc;
Dd = Dc;
