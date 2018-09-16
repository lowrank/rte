classdef tracer < handle
    properties
        id
    end
    
    methods
        % constructor.
        function obj = tracer(angle, node, elem, neigh)
            obj.id = TracerWrapper('new', angle, node, elem, neigh);
        end
        
        % destructor.
        function delete(obj)
            TracerWrapper('delete', obj.id);
        end
        
        % display the memory usage for the whole object. It might be quite
        % incorrect though. It also displays the ray path for specific
        % node. It is for debug purpose.
        function disp(obj, l, m)
            if (nargin < 3)
                TracerWrapper('disp', obj.id, -1, -1);
            else
                if (l >= 0) && (m >= 0)
                    TracerWrapper('disp', obj.id, l, m);
                else
                    error('Tracer:disp', 'Invalid input, non-negative required.');
                end
            end
        end
        
        % returns the solution from boundary contribution for pure
        % transport.
        function v = boundary_transport(obj, node, elem, sigma_a, u)
            v = TracerWrapper('boundary_transport', obj.id, node, elem, sigma_a, u);
        end
        
        % returns the solution from interior source contribution for pure
        % transport.
        function v = interior_transport(obj, node, elem, sigma_a, u)
            v = TracerWrapper('interior_transport', obj.id, node, elem, sigma_a, u);
        end
    end
    
end

