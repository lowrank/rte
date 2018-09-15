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
        
        function disp(obj, l)
            if (nargin < 2)
                TracerWrapper('disp', obj.id, -1);
            else
                TracerWrapper('disp', obj.id, l);
            end
        end
        
        function v = boundary_transport(obj, node, elem, sigma_a, u)
            v = TracerWrapper('boundary_transport', obj.id, node, elem, sigma_a, u);
        end
    end
    
end

