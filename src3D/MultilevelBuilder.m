classdef (Sealed) MultilevelBuilder < handle
   %MULTILEVELBUILDER Constructs the multi-level data structure.
   % This class builds a list of increasingly-finer levels to be used in
   % the multigrid cycle.
   
   %======================== METHODS =================================
   methods
      function levels = build(obj, options)
         % Build the list of levels from options.
         levels = cell(options.numLevels, 1);
         % Coarsest level
         n = options.nCoarsest;
         levels{1} = Level.newCoarsestLevel(options.domainSize, n, ...
            options.operator, options.smoother);
         % Increasingly-finer levels
         for l=2:options.numLevels
            n = 2*n;
            lev = Level.newLevel(options.domainSize, n, ...
               options.operator, options.smoother, ...
               levels{l-1}, options.interpolator, options.restrictor);
            % Initialize finest right-hand side
            if l == options.numLevels
               % Interior RHS
               MultilevelBuilder.setRhsValues(lev, ...
                  2:n(1)-1, 2:n(2)-1, 2:n(3)-1, options.f);
               % Boundary RHS
               MultilevelBuilder.setRhsValues(lev, ...
                  [1 n(1)], 1:n(2), 1:n(3), options.g);
               MultilevelBuilder.setRhsValues(lev, ...
                  1:n(1), [1 n(2)], 1:n(3), options.g);
               MultilevelBuilder.setRhsValues(lev, ...
                  1:n(1), 1:n(2), [1 n(3)], options.g);
            end
            levels{l} = lev;
         end
      end
   end
   %======================== PRIVATE METHODS =========================
   methods (Static, Access = private)
      function setRhsValues(lev, i, j, k, f)
         % Set the values of indices (i,j,k) of a level's RHS vector to
         % the function f, evaluated at the corresponding gridpoint
         % locations.
         [xInterior,yInterior,zInterior] = lev.location(i,j,k);
         % Convert singleton x,y,z vectors to 3-D matrices
         [X,Y,Z] = ndgrid(xInterior, yInterior, zInterior);
         lev.f(i,j,k) = f(X,Y,Z);
      end
   end

end