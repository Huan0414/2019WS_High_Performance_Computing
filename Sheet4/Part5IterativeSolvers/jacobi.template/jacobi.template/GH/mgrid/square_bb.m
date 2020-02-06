function square_bb(sx,sy)
% Square: 
%   flatpak run org.octave.Octave <filename>
%      or
%   octave --no-window-system --no-gui  -qf <filename>

if nargin<2 
    sx = 2;
    sy = 2;
end

g = generate_rectangle_subdomains(sx,sy)
% pdegplot(g,'EdgeLabels','on','FaceLabels','on')

%% Generate mesh from geometry
%
% [p,e,t] = initmesh(g,'hmax',1);  % works correctly
[p,e,t] = initmesh(g,'hmax',0.1);  % works ??
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
nrefine=0;                 % 
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

fname = [mfilename,'_',num2str(sx*sy)];
ascii_write_mesh( p, t, e, fname);
ascii_write_subdomains( p, t, e, fname);

% tmp=t(1:3,:)

