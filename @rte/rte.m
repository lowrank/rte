classdef rte < handle
    
    properties
        % structs
        rays
        
        % vectors
        nodes
        elems
        neigs
        edges
        segms
        angles
        bc
        
        boundarySource
        sigmaT
        sigmaS
        
        % scalars
        nAngle
        nPoint
    end
    
    methods
        function obj = rte(opt)
            if isfield(opt, 'angle')
                obj.nAngle = opt.angle;
            else 
                obj.nAngle = 128;
            end
            
            if isfield(opt, 'nodes') && isfield(opt, 'minArea')
                meshGenerator(obj, opt.nodes, opt.minArea);
            else 
                % default choice.
                meshGenerator(obj, [0 0; 1 0; 1 1; 0 1]', 0.001); 
            end
            
            obj.angles  = linspace(0, 2*pi, obj.nAngle + 1) ;
            obj.angles  = obj.angles(1:end-1);
            
            obj.rays    = tracer(obj.angles, obj.nodes, obj.elems, obj.neigs);
            
            obj.bc = unique(obj.segms);            
        end
        
        function meshGenerator(obj, nodes, minArea)
            hull = reshape(nodes, 2 * size(nodes, 2), 1);
            idx  = circshift(reshape(repmat(0:size(nodes, 2)-1, 2, 1),...
                    2 * size(nodes, 2), 1), [-1,1]);
                
            mesh = TriangleMesh();
            mesh.set_points_tri(hull);
            mesh.set_facets_tri(idx);
            mesh = mesh.build_tri(); % not ready to go
            mesh = mesh.refine_tri(sprintf('q34.0a%f', minArea));  

            [obj.nodes, obj.segms, obj.elems, obj.edges, obj.neigs] = ...
                mesh.getData_tri();

            obj.nPoint = size(obj.nodes, 2);
            
        end
        
        function setBoundaryCondition(obj, boundaryFunctionHandle)
            obj.boundarySource = zeros(obj.nAngle, obj.nPoint);
            
            for j = 1:length(obj.bc)
                for i = 1:obj.nAngle
                    obj.boundarySource(i, obj.bc(j)) = ...
                        boundaryFunctionHandle( ... 
                        obj.nodes(1, obj.bc(j)), ...
                        obj.nodes(2, obj.bc(j)), ...
                        obj.angles(i) ...
                        );
                end
            end
        end        
    end
    
end

