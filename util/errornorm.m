function e = errornorm(level, u)
%ERROR_NORM Error norm at a certain coarsening level.
%   E = ERROR_NORM(LEVEL,U) computes the grid-scale L2 residual norm
%   |F-L(U)|_2, where F and L are stored in the LEVEL structure.

r = level.residual(u);
e = norm(r(:))/sqrt(numel(r));

end