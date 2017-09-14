classdef tracer < handle
    properties
        id
    end
    
    methods
        function obj = tracer(angle, node, elem, neigh)
            obj.id = TracerWrapper('new', angle, node, elem, neigh);
        end
        
        function delete(obj)
            TracerWrapper('delete', obj.id);
        end
        
        function disp(obj)
            TracerWrapper('disp', obj.id);
        end
    end
    
end

