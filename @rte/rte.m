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
    
    methods (Access = public)
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
        
        function setCoefficents(obj, sigmaTFunctionHandle, sigmaSFunctionHandle)
            obj.sigmaT = sigmaTFunctionHandle(obj.nodes);
            obj.sigmaS = sigmaSFunctionHandle(obj.nodes);
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
        
        function setInteriorCondition(obj, interiorFunctionHandle)
            obj.interiorSource = zeros(obj.nAngle, obj,nPoint);
            
            for j = 1:obj.nPoint
                for i = 1:obj.nAngle
                    obj.interiorSource(i, j) = ...
                        interiorFunctionHandle( ... 
                        obj.nodes(1, j), ...
                        obj.nodes(2, j), ...
                        obj.angles(i));
                end
            end
        end
        
        function ForwardSolve(obj)
            % in the real case, the momentum is not independent quantity.
            % We still have to solve a large system to get the iteration
            % converged.
            tic; u = obj.rays.boundary_transport(...
                obj.nodes, obj.elems, obj.sigmaT, obj.boundarySource); toc;
            boundaryFluence = (sum(u, 1)/ obj.nAngle)';
            
            
            forwardMap = @(X) ( X - obj.sigmaS' .* (sum( obj.rays.interior_transport(...
                obj.nodes, obj.elems, obj.sigmaT,  repmat(X', obj.nAngle, 1)) ,1 )/ obj.nAngle)' );
            
            tic; x = gmres(forwardMap, boundaryFluence, 10, 1e-10, 400); toc;
            
            trisurf(obj.elems(1:3,:)', obj.nodes(1,:), obj.nodes(2,:), ...
                x', 'EdgeColor', 'None');
            
            shading interp;view(2);colorbar;colormap jet;
            
        end
    end
    
    methods (Access = private)
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
    end
    
    methods (Static)
        % computes the Fourier transform directly.
        %
        % caution: 128 directions might only resolve anisotropy at g = 0.8
        % and 256 directions only can resolve g = 0.9. For larger g, it
        % should be noted that such forward peaking case can be well
        % approximated by some other methods and not necessary to pursue
        % this approach any more.
        function f = HenyeyGreenstein(g, nAngle)
            theta = linspace(0, 2 * pi, nAngle + 1);
            theta = theta(1:end- 1);
            f = fft(1/(2*pi) * (1 - g^2) ./ (1 + g^2 - 2 * g * cos(theta)));
        end
    end
        
    
end

