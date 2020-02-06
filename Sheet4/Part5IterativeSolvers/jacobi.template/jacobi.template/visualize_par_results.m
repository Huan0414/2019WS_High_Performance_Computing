%% Visualize results
%
%   flatpak run org.octave.Octave <filename>
%      or
%   octave --no-window-system --no-gui  -qf <filename>
%
%      or
%   
%   matlab -nosplash -nodesktop -r 'try visualize_par_results(4); catch; end; quit'
%
function visualize_par_results(nprocs)
%%
if nargin<1 
    nprocs = 4;
end
fprintf('# procs = %d\n',nprocs)

pre  = 'uv_';
post = '.txt';

xc = []; nnodes = [];
ia = []; nelems = [];
v  = [];
node_offset = 0;
elem_offset = 0;
for rank=0:nprocs-1
    fname = [pre,num2str(rank,'%2u'),post];
    [lxc,lia,lv] = ascii_read_meshvector(fname);
%     whos lxc lia lv
    nnodes = [nnodes size(lxc,1)];
    nelems = [nelems size(lia,1)];
    %[xc,ia,v]
    xc = [xc; lxc];
    v  = [v ; lv ];
    ia = [ia; lia+node_offset];
%     node_offset
%     lia = lia + node_offset
%     ia = [ia; lia];
    % index offsets for next subdomain
    node_offset = node_offset + nnodes(end);
    elem_offset = elem_offset + nelems(end);
end

% fname = 'uv.txt';
% [xc,ia,v] = ascii_read_meshvector(fname);

h = trisurf(ia, xc(:,1), xc(:,2), v);
xlabel('x'),ylabel('y'),zlabel('z')

shading interp

waitfor(h)                     % wait for closing the figure
