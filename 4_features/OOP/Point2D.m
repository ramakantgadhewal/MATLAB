classdef Point2D < handle
  properties( Access = private )
    % provide default values
    x, y
  end
  properties(Dependent)
    r
  end
  methods
    function obj = Point2D(x0, y0)
      % Constructor
      if nargin == 0
        obj.x = 0.0;
        obj.y = 0.0;
      elseif nargin == 2
        obj.x = x0;
        obj.y = y0;
      else
        error('Error: Wrong input arguments.');
      end
    end
    function r = get.r(obj)
      % calculate the values of dependent properties
      r = sqrt(obj.x^2 + obj.y^2);
      fprintf('  Attention: get.r in Point2D called.\n\n');
    end
    function normalize(obj)
      if obj.r > 0
        obj.x = obj.x/obj.r;
        obj.y = obj.y/obj.r;
      else
        error('Error: %s (%f,%f) cannot be normalized.', ...
            class(obj), obj.x, obj.y);
      end
    end
    function print(obj)
      disp(['x =', num2str(str.x)]);
      disp(['y =', num2str(str.y)]);
    end
  end
end
