classdef CsPoint
    % This class should combine together MeasEntity and MoveEntity, because
    % they naturally go together. 
    % 
    % We should still only have to provide xr, yr. Perhaps at this point
    % the easiest thing to do is let ME and MVE be function handles. 
   
    
    properties
        me_factory;
        mve_factory;
    end
    
    methods
        function obj = CsPoint(me_factory, mve_factory)
           obj.me_factory = me_factory;
           obj.mve_factory = mve_factory;
        end
        
        function ME = meas_entity(self, xr, yr, N, ind)
            ME = self.me_factory.meas_entity(xr, yr, N, ind)
        end
        
        function MVE = move_entity(self, xr, yr, N, ind)
            MVE = self.mv_factory.move_entity(xr, yr, N, ind)
        end
    end
    
end

