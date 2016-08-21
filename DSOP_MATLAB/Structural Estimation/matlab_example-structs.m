% Taken from iso2mesh.sourceforge.net/cgi-bin/index.cgi?jsonlab/Doc/Examples

a = struct('node',[1 9 10; 2 1 1.2], 'elem',[9 1;1 2;2 3],...
           'face',[9 01 2; 1 2 3; NaN,Inf,-Inf], 'author','FangQ');

% Now can reference the elements like so:

a.node
% ans =
% 1.0000    9.0000   10.0000
% 2.0000    1.0000    1.2000

a.elem
% ans =
% 9     1
% 1     2
% 2     3



















%
