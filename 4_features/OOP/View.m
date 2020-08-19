classdef View < handle
  properties
    hFig
    hEdit
  end
  properties(Dependent)
    text
  end
  methods
    function obj = View()
      % the two properties are both objects
      obj.hFig = figure();
      obj.hEdit = uicontrol('style', 'edit', 'parent', obj.hFig)
    end
    function str = get.text(obj)
      str = get(obj.hEdit, 'String');
    end
  end
end
