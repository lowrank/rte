classdef rte < handle
    % Radiative transfer class for 2D.
    properties (Access = public)
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
        dtheta
        g
    end
    
    properties (Access = public)
        pFHG
    end
    
    methods (Access = public)
        function obj = rte(opt)
            if isfield(opt, 'angle')
                obj.nAngle = opt.angle;
            else 
                obj.nAngle = 64;
            end
            
            if isfield(opt, 'nodes') && isfield(opt, 'minArea')
                meshGenerator(obj, opt.nodes, opt.minArea);
            else 
                % default choice.
                meshGenerator(obj, [0 0; 1 0; 1 1; 0 1]', 0.001); 
            end
            
            obj.angles  = linspace(0, 2*pi, obj.nAngle + 1) ;
            obj.angles  = obj.angles(1:end-1);
            obj.dtheta  = 2 * pi / obj.nAngle;
            
            obj.rays    = tracer(obj.angles, obj.nodes, obj.elems, obj.neigs);
            
            obj.bc = unique(obj.segms);            
            
            if isfield(opt, 'anisotropy')
                obj.g = opt.anisotropy;
            else
                obj.g = 0.0;
            end
            
            
            obj.pFHG = obj.HenyeyGreenstein(obj.g, obj.nAngle);
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
        
        function setAnisotropy(obj, gprime)
            obj.g = gprime;
            obj.pFHG = obj.HenyeyGreenstein(obj.g, obj.nAngle);
        end
        
        function x = ForwardSolve(obj)
            % There are redundant operations like reshaping in the
            % implementation which can be opt out.
            boundaryContrib = obj.rays.boundary_transport(...
                obj.nodes, obj.elems, obj.sigmaT, obj.boundarySource); 
            
            forwardMap = @(X) (X - obj.mapping(reshape(X, obj.nAngle, obj.nPoint)));
            
            tic; x = gmres(forwardMap, reshape(boundaryContrib, obj.nAngle * obj.nPoint, 1), 10, 1e-10, 400); toc;
            
            x = reshape(x, obj.nAngle, obj.nPoint);
            
        end
        
        function plot(obj, x)
            if (size(x, 2) > 1)
                trisurf(obj.elems(1:3,:)', obj.nodes(1,:), obj.nodes(2,:), ...
                sum(x, 1)/obj.nAngle, 'EdgeColor', 'None');
                shading interp;view(2);colorbar;colormap jet;
            else
                trisurf(obj.elems(1:3,:)', obj.nodes(1,:), obj.nodes(2,:), ...
                x, 'EdgeColor', 'None');
                shading interp;view(2);colorbar;colormap jet;
            end
                
            
        end
    end
    
    methods (Access = public)
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
        
        function Y = mapping(obj, X)
            % input X is a large column vector now.
            rX  = reshape(X, obj.nAngle, obj.nPoint);
            sca = ifft(bsxfun(@times, obj.pFHG', fft(rX))) * obj.dtheta;
            sca = bsxfun(@times, obj.sigmaS', sca');
            Y = reshape(  obj.rays.interior_transport(...
                obj.nodes, obj.elems, obj.sigmaT,  sca'), obj.nAngle * obj.nPoint, 1);
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
        
        % question: how to make this more efficient?
        function f = HenyeyGreenstein(g, nAngle)
            theta = linspace(0, 2 * pi, nAngle + 1);
            theta = theta(1:end- 1);
            f = fft(1/(2*pi) * (1 - g^2) ./ (1 + g^2 - 2 * g * cos(theta)));
        end
    end
        
    
end

