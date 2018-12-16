
function [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)
% cgsolve.m
%
% Solve a symmetric positive definite system Ax = b via conjugate gradients.
%
% Usage: [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose)
%
% A - function handle which computes A_matrix*x.
%
% b - N vector
%
% tol - Desired precision.  Algorithm terminates when 
%    norm(Ax-b)/norm(b) < tol .
%
% maxiter - Maximum number of iterations.
%
% verbose - If 0, do not print out progress messages.
%    If and integer greater than 0, print out progress every 'verbose' iters.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

  x = zeros(length(b),1);
  r = b; % = b - A*x0, x0 = 0.
  d = b;
  delta = r'*r;
  delta_init = b'*b;
  
  numiter = 0;
  bestx = x;
  bestres = sqrt(delta/delta_init);
  while ((numiter < maxiter) && (delta > tol^2*delta_init))
    q = A(d);
    alpha = delta/(d'*q);
    x = x + alpha*d;
    
    if (mod(numiter+1,50) == 0)
      r = b - A(x);
    else
      r = r - alpha*q;
    end
    
    deltaold = delta;
    delta = r'*r;
    beta = delta/deltaold;
    d = r + beta*d;
    numiter = numiter + 1;
    if (sqrt(delta/delta_init) < bestres)
      bestx = x;
      bestres = sqrt(delta/delta_init);
    end
    
    if ((verbose) && (mod(numiter,verbose)==0))
      fprintf('cg: Iter = %d, Best residual = %8.3e, Current residual = %8.3e\n', ...
        numiter, bestres, sqrt(delta/delta_init));
    end
    
  end
  
  if (verbose)
    fprintf('cg: Iterations = %d, best residual = %14.8e\n', numiter, bestres);
  end
  x = bestx;
  res = bestres;
  iter = numiter;
end