function f = flops(fl)
%FLOPS Get or set the global flopcount variable.
%   Returns the current flopcount. FLOPS(F) sets flopcount to F.
global flopcount

if nargin == 0
   f = flopcount;
else
   flopcount = fl;
   if nargout == 1; f = fl; end
end

end