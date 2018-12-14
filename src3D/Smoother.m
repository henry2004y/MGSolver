classdef (AllowedSubclasses = {?GaussSeidelSmoother,?JacobiSmoother}) ...
      Smoother < handle
   %SMOOTHER General interface of the smoother abstract class
   %   This is the general interface which allows different implementations
   %   of relaxation schemes.
   
   %======================== MEMBERS =================================
   properties(GetAccess = public, SetAccess = private)
      name {char}
   end
   %======================== CONSTRUCTORS ============================
   methods
      function obj = Smoother(name)
         obj.name = name;
      end
   end
   %======================== METHODS =================================
   methods(Abstract)
      u = relax(obj,level,u)
   end
   
   methods
      u = relax_conv(obj,level,u)
   end
   
end

