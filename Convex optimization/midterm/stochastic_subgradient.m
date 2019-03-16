function x = stochastic_subgradient(G, prox, C, N, x0)

x = [x0];
for i = 1:1000
    i
    eta  = C/sqrt(i);
    xi = randi([1, N]);
    x1 = prox(x0 - eta*G(x0, xi), eta);
    x = [x x1];
    x0 = x1;
end