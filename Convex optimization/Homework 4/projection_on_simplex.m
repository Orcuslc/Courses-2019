function w = projection_on_simplex(v)
% Projection onto Probability simplex
% copied from `ProjectOntoL1Ball`;
% Reference: https://eng.ucmerced.edu/people/wwang5/papers/SimplexProj.pdf

% % PROJECTONTOL1BALL Projects point onto L1 ball of specified radius. 
% w = ProjectOntoL1Ball(v, b) returns the vector w which is the solution 
% to the following constrained minimization problem: 
% % min ||w - v||_2 
% s.t. ||w||_1 <= b. 
% % That is, performs Euclidean projection of v to the 1-norm ball of radius 
% b. 
% % Author: John Duchi (jduchi@cs.berkeley.edu) 

u = sort(abs(v),'descend'); 
sv = cumsum(u);
rho = find(u+(1-sv)./(1:length(u))' > 0, 1, 'last'); 
lambda = 1/rho*(1-sv(rho));
% theta = max(0, (sv(rho) - b) / rho); 
% w = sign(v) .* max(abs(v) - theta, 0);
w = max(v+lambda, 0);


