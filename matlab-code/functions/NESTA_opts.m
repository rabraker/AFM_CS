%               This is a structure that contains additional options,
%               some of which are optional.
%               The fieldnames are case insensitive.  Below
%               are the possible fieldnames:
%               
%               opts.xplug - the first guess for the primal prox-function, and
%                 also the initial point for xk.  By default, xplug = At(b)
%               opts.U and opts.Ut - Analysis/Synthesis operators
%                 (either matrices of function handles).
%               opts.normU - if opts.U is provided, this should be norm(U)
%               opts.maxiter - max number of iterations in an inner loop.
%                 default is 10,000
%               opts.TolVar - tolerance for the stopping criteria
%               opts.stopTest - which stopping criteria to apply
%                   opts.stopTest == 1 : stop when the relative
%                       change in the objective function is less than
%                       TolVar
%                   opts.stopTest == 2 : stop with the l_infinity norm
%                       of difference in the xk variable is less
%                       than TolVar
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
%               opts.fid - if this is 1 (default), the display is
%                   the usual Matlab screen.  If this is the file-id
%                   of a file opened with fopen, then the display
%                   will be redirected to this file.
%               opts.errFcn - if this is a function handle,
%                   then the program will evaluate opts.errFcn(xk)
%                   at every iteration and display the result.
%                   ex.  opts.errFcn = @(x) norm( x - x_true )

function opts = NESTA_opts(varargin)
    
   p = inputParser();
   U_default_fun = @(x) x;
   p.addParameter('fid', 1);
   p.addParameter('Verbose', true);
  
  
   p.addParameter('MaxIntIter',5);
   p.addParameter('TypeMin','L1');
   p.addParameter('TolVar',1e-5);
   
   p.addParameter('xplug',[]);
   p.addParameter('normU',[]);  % so we can tell if it's been set
   
   p.addParameter('errFcn',[]);

   p.addParameter('stopTest', 1);
   p.addParameter('U', U_default_fun);
   p.addParameter('Ut', U_default_fun);
  
   p.parse(varargin{:});
   
   opts = p.Results;
   
   if ~isequal(opts.U, U_default_fun)
       opts.U_userSet = true;
   else
       opts.U_userSet = false;
       opts.normU = 1;
   end
end
