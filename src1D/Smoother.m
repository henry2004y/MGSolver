classdef (AllowedSubclasses = {?GaussSeidelSmoother,?JacobiSmoother}) ...
      Smoother < handle
   %SMOOTHER General interface of the smoother abstract class
   %   This is the general interface which allows different implementations
   %   of relaxation schemes.
   
   
   methods(Abstract)
      u = relax(obj,level,u)
   end
   
end

