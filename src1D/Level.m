classdef (Sealed) Level < handle
   %LEVEL A single level in the multi-level cycle.
   % This class holds all data and operations pertinent to a single level
   % in the multi-level cycle: right-hand-side, residual computation and
   % single-level processes such as relaxation.
   
   %======================== MEMBERS =================================
   properties (GetAccess = public, SetAccess = public)
      f % Right-hand-side of both the interior equations & B.C.
   end
   
   properties (GetAccess = public, SetAccess = private)
      domainSize % Size of domain
      h % Mesh-size (same in all directions)
      n % Grid array size vector
   end
   properties (GetAccess = private, SetAccess = private)
      coarseLevel % Next-coarser level
      operator % Compute the discrete operator L(u) at this level
      smoother % Relaxation scheme
      interpolator % Interpolates corrections from fineLevel
      restrictor % Restricts residuals to fineLevel
   end
   
   %======================== CONSTRUCTORS ============================
   methods (Access = private)
      function obj = Level(domainSize, n, operator, smoother, ...
            coarseLevel, interpolator, restrictor)
         % Initialize a level.
         obj.domainSize = domainSize;
         obj.n = n+1;
         hVector = domainSize./n;
         if std(hVector) > eps
            error(['Incompatible domain size [%f] and #intervals '...
               '[%d]: meshsize must be the same in all directions'], ...
               domainSize(1), n(1));
         end
         obj.h = hVector;
         obj.f = zeros(obj.n,1);
         obj.operator = operator(obj);
         obj.smoother = smoother;
         obj.coarseLevel = coarseLevel;
         obj.interpolator = interpolator;
         obj.restrictor = restrictor;
      end
   end
   
   methods (Static)
      function obj = newLevel(domainSize, n, operator, smoother, ...
            coarseLevel, interpolator, restrictor)
         % A factory method of the next-finer level over COARSELEVEL,
         % with an NxN grid of a domain of size DOMAINSIZExDOMAINSIZE,
         % discrete operator OPERATOR, a relaxation scheme SMOOTHER and
         % inter-grid transfers INTERPOLATOR and RESTRICTOR. The
         % right-hand-side is initialized to zero.
         obj = Level(domainSize, n, operator, smoother, ...
            coarseLevel, interpolator, restrictor);
      end
      
      function obj = newCoarsestLevel(domainSize, n, operator, smoother)
         % A factory method of the coarsest level, with an NxN grid of a
         % domain of size DOMAINSIZExDOMAINSIZE, a discrete operator
         % OPERATOR and a relaxation scheme SMOOTHER.
         obj = Level(domainSize, n, operator, smoother, [], [], []);
      end
   end
   
   %======================== METHODS =================================
   methods
      function r = residual(obj, u)
         % Compute the residual F-L(U) for a function U at this level.
         r = obj.f - obj.L(u);
      end
      function v = relax(obj, u)
         % Perform a relaxation sweep. Delegates to the smoother with a
         % call-back to this level.
         v = obj.smoother.relax(obj, u);
      end
      function u = interpolate(obj, uc)
         % Interpolate the correction uc from the next-coarser level.
         u = obj.interpolator.interpolate(obj.coarseLevel, obj, uc);
      end
      function fc = restrict(obj, f)
         % Restrict the residual FC to the next-coarser level.
         fc = obj.restrictor.restrict(obj.coarseLevel, obj, f);
      end
      function x = location(obj, i)
         % Return gridpoint locations at indices (I).
         x = obj.h*(i-1);
      end
      function result = L(obj, u)
         % Apply the discrete operator L to a function U.
         result = obj.operator.L(u);
      end
      function handle = plot(obj, u)
         % Plot the discrete function U on the grid of this level.
         x = obj.location(1:obj.n(1));
         handle = plot(x,u,'--+');
      end
   end
   
end