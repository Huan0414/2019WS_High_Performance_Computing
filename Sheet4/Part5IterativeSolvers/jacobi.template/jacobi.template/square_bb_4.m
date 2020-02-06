% Square: 
%   flatpak run org.octave.Octave <filename>
%      or
%   octave --no-window-system --no-gui  -qf <filename>

clear all
clc
% %% L-shape
% g=[2 0 2 0 0 1 0;        % #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
%    2 2 2 0 1 1 0;
%    2 2 1 1 0.5 1 0;
%    2 1 1 0.5 2 1 0;
%    2 1 0 2 2 1 0;
%    2 0 0 2 0 1 0]';

%% square
% g=[2 0 1 0 0 1 0;        % #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
%    2 1 1 0 1 1 0;
%    2 1 0 1 1 1 0;
%    2 0 0 1 0 1 0]';

% %% 2 squares
% g=[2 0 1 0 0 1 0;        % 1 #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
%    2 1 1 0 1 1 2;
%    2 1 0 1 1 1 0;
%    2 0 0 1 0 1 0;
%    2 1 2 0 0 2 0;        % 2 #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
%    2 2 2 0 1 2 0;
%    2 2 1 1 1 2 0
%    ]';

%% 4 squares
g=[2 0 1 0 0 1 0;        % 1 #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
   2 1 1 0 1 1 2;
   2 1 0 1 1 1 3;
   2 0 0 1 0 1 0;
   2 1 2 0 0 2 0;        % 2 #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
   2 2 2 0 1 2 0;
   2 2 1 1 1 2 4;
%    2 1 1 1 0 2 1;
%    2 0 1 1 1 3 1;        % 3 #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
   2 1 1 1 2 3 4;
   2 1 0 2 2 3 0;
   2 0 0 2 1 3 0;
%    2 1 2 1 1 4 2;        % 4 #vertices,v_1x, v_2x, v_1y, v_2y, subdomain_left, subdomain_right
   2 2 2 1 2 4 0;
   2 2 1 2 2 4 0
%    2 1 1 2 1 4 3
   ]';

%% Generate mesh from geometry
%
[p,e,t] = initmesh(g,'hmax',1);  % works correctly
% p(1,15) = 1.51;                  %% angle in trangle > pi/2 ==> now the second refinement produces irregular meshes!
% p(1,15) = 1.7;                   %% angle in trangle > pi/2 ==> now the second refinement produces irregular meshes!

%  ?? 
%  https://de.mathworks.com/help/pde/ug/mesh-data-pet-triples.html
%  generateMesh(...)
%  mesh2Pet(...)
%
% [p,e,t] = initmesh(g); % problems in solution after 2 refinements
% [p,e,t] = initmesh(g,'hmax',0.5); % problems in solution after 2 refinements (peaks with h=0.5, oscillations in (1,1) for h=0.1
% [p,e,t] = initmesh(g,'hmax',0.1/4); % no problems in solution with 0 refinemnet steps

%% Show mesh
pdemesh(p,e,t)
% pdemesh(p,e,t,'NodeLabels','on')

%% Improve mesh
% min(pdetriq(p,t))
% p = jigglemesh(p,e,t,'opt','minimum','iter',inf); 
% min(pdetriq(p,t))
% pdemesh(p,e,t)

%% Refine mesh, see comments in  "Generate mesh from geometry"
%
% nrefine=8;
nrefine=2;                 % 
for k=1:nrefine
    [p,e,t] = refinemesh(g,p,e,t);
%     p = jigglemesh(p,e,t,'opt','minimum','iter',inf);  % improve mesh
    min(pdetriq(p,t))
    fprintf('refinement: %i  nodes: %i     triangles: %i \n', k, size(p,2), size(t,2))
end
% figure; pdemesh(p,e,t,'NodeLabels','on')
%

%% GH
% output from <https://de.mathworks.com/help/pde/ug/initmesh.html initmesh>
%
% coordinates  p: [2][nnode]
% connectivity t: [4][nelem]   with  t(4,:) are the subdomain numbers
% edges        e: [7][nedges]  boundary edges
%                              e([1,2],:) - start/end vertex of edge
%                              e([3,4],:) - start/end values
%                              e(5,:)     - segment number
%                              e([6,7],:) - left/right subdomain

ascii_write_mesh( p, t, e, mfilename);

ascii_write_subdomains( p, t, e, mfilename);


% tmp=t(1:3,:)

