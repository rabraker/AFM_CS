function [xk,niter] =NESTA_mine(A,At,b,muf,delta,opts)
% [xk,niter,residuals,outputData] =NESTA(A,At,b,muf,delta,opts)
%
% Solves a L1 minimization problem under a quadratic constraint using the
% Nesterov algorithm, with continuation:
%
%     min_x || U x ||_1 s.t. ||y - Ax||_2 <= delta
% 
% Continuation is performed by sequentially applying Nesterov's algorithm
% with a decreasing sequence of values of  mu0 >= mu >= muf
%
% The primal prox-function is also adapted by accounting for a first guess
% xplug that also tends towards x_muf 
%
% The observation matrix A is a projector
%
% Inputs:   A and At - measurement matrix and adjoint (either a matrix, in which
%               case At is unused, or function handles).  m x n dimensions.
%           b   - Observed data, a m x 1 array
%           muf - The desired value of mu at the last continuation step.
%                 A smaller mu leads to higher accuracy.
%           delta - l2 error bound.  This enforces how close the variable
%               must fit the observations b, i.e. || y - Ax ||_2 <= delta
%               If delta = 0, enforces y = Ax
%               Common heuristic: delta = sqrt(m + 2*sqrt(2*m))*sigma;
%               where sigma=std(noise).
%           opts -
%               This is a structure that contains additional options,
%               some of which are optional.
%               The fieldnames are case insensitive.  Below
%               are the possible fieldnames:
%               
%               opts.U and opts.Ut - Analysis/Synthesis operators as function
%               handles.
%               opts.normU - if opts.U is provided, this should be norm(U)
%                   otherwise it will have to be calculated via power iteration on 
%                   on U'*U (to estimate the largest singular value, which is 
%                            potentially expensive).
%               opts.MaxIntIter - number of continuation steps.
%                 default is 5
%               opts.maxiter - max number of iterations in an inner loop.
%                 default is 10,000
%               opts.TolVar - tolerance for the stopping criteria, i.e.,
%                             when the relative change in the objective 
%                             function is less than TolVar.
%               opts.TypeMin - if this is 'L1' (default), then
%                   minimizes a smoothed version of the l_1 norm.
%                   If this is 'tv', then minimizes a smoothed
%                   version of the total-variation norm.
%                   The string is case insensitive.
%               opts.Verbose - if this is 0 or false, then very
%                   little output is displayed.  If this is 1 or true,
%                   then output every iteration is displayed.
%                   If this is a number p greater than 1, then
%                   output is displayed every pth iteration.
%               opts.errFcn - if this is a function handle,
%                   then the program will evaluate opts.errFcn(xk)
%                   at every iteration and display the result.
%                   ex.  opts.errFcn = @(x) norm( x - x_true )
%
%  Outputs:
%           xk  - estimate of the solution x
%           niter - number of iterations
%
% This code is a heavily modified version of the published NESTA code
% availible at https://statweb.stanford.edu/~candes/nesta/
% 
% The original code was written by Jerome Bobin and Stephen Becker.
% 


  if nargin < 6 || (isempty(opts) && isnumeric(opts))
    opts = NESTA_opts();
  end
  
  % -- We can handle non-projections IF a (fast) routine for computing
  %    the psuedo-inverse is available.
  %    We can handle a nonzero delta, but we need the full SVD
  % Check if A is a partial isometry, i.e. if AA' = I
  z = randn(size(b));
  AAtz = A(At(z));
  if norm( AAtz - z )/norm(z) > 1e-8
    error('Measurement matrix A must be a partial isometry: AA''=I');
  end

  % -- Find a initial guess.
  % Use min-energy solution: x_ref = A'*inv(A*A')*b. Since we assume that
  % that AA'=I, then x_ref = A'*b.
  x_ref=At(b);
  opts.xplug = x_ref;

  % x_ref itself is used to calculate mu_0
  %   in the case that xplug has very small norm

  % use x_ref, not xplug, to find mu_0
  Ux_ref = opts.U(x_ref);
  switch lower(opts.TypeMin)
    case 'l1'
      mu0 = 0.9*max(abs(Ux_ref));
    case 'tv'
      mu0 = ValMUTv(Ux_ref);
  end

  opts = set_normU(opts);

  niter = 0;
  Gamma = (muf/mu0)^(1/opts.MaxIntIter);
  mu = mu0;
  Gammat= (opts.TolVar/0.1)^(1/opts.MaxIntIter);
  opts.TolVar = 0.1;
  
  for nl=1:opts.MaxIntIter
    
    mu = mu*Gamma;
    opts.TolVar=opts.TolVar*Gammat;    
    
    if opts.Verbose
      fprintf('   \nBeginning %s Minimization: mu = %g\n\n',opts.TypeMin,mu);
    end
    [xk, niter_int] = Core_Nesterov_mine(A, At, b, mu, delta, opts);
    
    opts.xplug = xk;
    niter = niter_int + niter;
    
  end


  %---- internal routine for setting mu0 in the tv minimization case
  function th=ValMUTv(x)

    N = length(x);
    n = floor(sqrt(N));
    Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
                  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
    Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
                  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
    D = sparse([Dh;Dv]);


    Dhx = Dh*x;
    Dvx = Dv*x;
    
    sk = sqrt(abs(Dhx).^2 + abs(Dvx).^2);
    th = max(sk);

  end

end %-- end of NESTA function

function opts = set_normU(opts)
  % -- If U was set by the user and normU not supplied, then calcuate norm(U)
  if opts.U_userSet && isempty(opts.normU)
    % simple case: U*U' = I or U'*U = I, in which case norm(U) = 1
    z = randn(size(opts.xplug));
    
    UtUz = opts.Ut(opts.U(z));
    
    if norm( UtUz - z )/norm(z) < 1e-8
      opts.normU = 1;
    else
      z = randn(size(Ux_ref));
    
      UUtz = opts.U(opts.Ut(z)); 
    
      if norm( UUtz - z )/norm(z) < 1e-8
        opts.normU = 1;
      end
    end
    
    if isempty(opts.normU)
      % have to actually calculate the norm
      [opts.normU,cnt] = my_normest(opts.U, opts.Ut, length(opts.xplug), 1e-3,30);
      if cnt == 30
        fprintf('Warning: norm(U) may be inaccurate\n'); 
      end
    end
    % opts.normU = normU;
  end
end
  
  
  
%%%%%%%%%%%% POWER METHOD TO ESTIMATE NORM %%%%%%%%%%%%%%%
% Copied from MATLAB's "normest" function, but allows function handles,
% not just sparse matrices.
function [e,cnt] = my_normest(S,St,n,tol, maxiter)
%MY_NORMEST Estimate the matrix 2-norm via power method.
  if nargin < 4, tol = 1.e-6; end
  if nargin < 5, maxiter = 20; end
  if isempty(St)
    St = S;  % we assume the matrix is symmetric;
  end
  x = ones(n,1);
  cnt = 0;
  e = norm(x);
  if e == 0, return, end
  x = x/e;
  e0 = 0;
  while abs(e-e0) > tol*e && cnt < maxiter
    e0 = e;
    Sx = S(x);
    if nnz(Sx) == 0
      Sx = rand(size(Sx));
    end
    e = norm(Sx);
    x = St(Sx);
    x = x/norm(x);
    cnt = cnt+1;
  end
end