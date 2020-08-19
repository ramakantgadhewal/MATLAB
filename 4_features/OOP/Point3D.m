classdef Point3D < Point2D % use '<' to represent inheritance
  properties( Access = private )
    z
  end
  properties(Dependent)
    r
  end
  methods
    function obj = Point3D(x0, y0, z0)
      if nargin == 0
%        obj = obj@Point2D(0.0, 0.0);
        obj.x = 0.0;
        obj.y = 0.0;
        obj.z = 0.0;
      elseif nargin == 3
        obj = obj.Point2D(x0, y0);
        obj.z = z0;
      else
        error('Error: Wrong input arguments.');
      end
    end
    function r = get.r(obj)
      % calculate the values of dependent properties
      r = sqrt(obj.x^2 + obj.y^2 + obj.z^2);
      fprintf('  Attention: get.r in Point3D called.\n\n');
    end
    function print(obj)
      print@Point2D(obj);
      disp(['z =', num2str(str.z)]);
    end
  end
end
