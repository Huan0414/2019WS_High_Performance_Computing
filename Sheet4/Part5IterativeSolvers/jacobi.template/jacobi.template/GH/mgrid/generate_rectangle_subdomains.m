function [ g ] = generate_rectangle_subdomains( sx, sy )
%Generates the geometrical description of rectangle [0,sx]x[0,sy] 
% into sx x sy subdomains
% 
%% geometric description of domain,
%  see "Decomposed Geometry Data Structure"
%  https://de.mathworks.com/help/pde/ug/create-geometry-at-the-command-line.html#bulfvw4-1


nsubs  = sx*sy;
nedges = (sx)*(sy+1)+(sx+1)*sy;
g = zeros(nedges,7);
g(:,1) = 2;               % all edges are straight lines

e = 1;                    % edge index

% first the edges in x direction
for y=0:sy
    for x=1:sx
        g(e,2) = x-1;             % v_1x
        g(e,3) = x;               % v_2x
        g(e,4) = y;               % v_1y
        g(e,5) = y;               % v_2y
        g(e,6) = x + y*sx;        % subdomain_left
        g(e,7) = x + (y-1)*sx;    % subdomain_right
        e = e+1;
    end
end

% second the edges in y direction
for y=1:sy
    xsub = 0;
    for x=0:sx
        g(e,2) = x;               % v_1x
        g(e,3) = x;               % v_2x
        g(e,4) = y-1;             % v_1y
        g(e,5) = y;               % v_2y
        g(e,6) = xsub;
        xsub = x + 1 + (y-1)*sx;
        g(e,7) = xsub;
        e = e+1;
    end
    g(e-1,7) = 0;
end

% set all subdomain indices (left/right) with nsubs<0 or nd>nsubs to 0.
idx = (g(:,6)<0) | (nsubs<g(:,6));
g(idx,6) = 0;
idx = (g(:,7)<0) | (nsubs<g(:,7));
g(idx,7) = 0;

g=g.';        % Dimensions!

%%
% pdegplot(g,'EdgeLabels','on','FaceLabels','on')
% [p,e,t] = initmesh(g,'hmax',0.1);
% pdemesh(p,e,t)


