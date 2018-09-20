clear;clc;

opt = struct('anisotropy', 0.5, 'angle', 16, ...
    'nodes', [0 0; 1 0; 1 1; 0 1]', 'minArea', 0.004);

obj = rte(opt);

f = @(x,y,v) (1);
sigmaS = @(x) (1 + 0.4* x(1, :));
sigmaT = @(x) (0.4 *  x(1, :) + 1 + 0.2 * x(2,:));

obj.setBoundaryCondition(f);
obj.setCoefficents(sigmaT, sigmaS);

x1 = obj.ForwardSolve();

% then we compute the adjoint equation.


g = @(x,y ,v) (x);
obj.setBoundaryCondition(g);
obj.setCoefficents(sigmaT, sigmaS);

x2 = obj.ForwardSolve();

% flip directions.
nAngle = obj.nAngle;

x3 = [x2(nAngle/2+1:end, :) ; x2(1:nAngle/2, :)]; % flipped solution (true solution)

% on boundary.
Lbd = size(obj.segms, 2);

s = 0;


for i = 1:Lbd
    cur_edge = obj.segms(:, i);
    if obj.nodes(1, cur_edge(1)) == 1.0 && obj.nodes(1, cur_edge(2)) == 1.0
        % on the edge of right side
        n = [1, 0];
    elseif obj.nodes(1, cur_edge(1)) == 0.0 && obj.nodes(1, cur_edge(2)) == 0.0
        % on the edge of left side
        n = [-1, 0];
        
    elseif obj.nodes(2, cur_edge(1)) == 1.0 && obj.nodes(2, cur_edge(2)) == 1.0
        % on the edge of top side
        n = [0, 1];
    elseif obj.nodes(2, cur_edge(1)) == 0.0 && obj.nodes(2, cur_edge(2)) == 0.0
        % on the edge of bottom side
        n = [0, -1];
    end
    
    l = norm(obj.nodes(:, cur_edge(1))   - obj.nodes(:, cur_edge(2)));
    for j = 1:nAngle
        
        w = [cos((j-1) * obj.dtheta), sin((j-1) * obj.dtheta)];
        
        s = s + 0.5 * x1(j, cur_edge(1)) * x3(j, cur_edge(1)) * (n * w') * obj.dtheta * l;
        s = s + 0.5 * x1(j, cur_edge(2)) * x3(j, cur_edge(2)) * (n * w') * obj.dtheta * l; 
        
    end
    
end

disp(sprintf('error of adjoint %f\n', s));


