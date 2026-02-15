function dYdt = ode45StateRhs3DFrame(t, y, M, C, K, groundLoad, agFun)
n = size(M,1);
u = y(1:n);
v = y(n+1:end);
ag = agFun(t);
rhs = -C*v - K*u + groundLoad*ag;
a = M \ rhs;
dYdt = [v; a];
end
