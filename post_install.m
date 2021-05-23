function post_install(desc)
%disp('$$$ post_install called');
%instDir=[desc.dir '/inst']
%mkdir(instDir);
copyfile('*.m', desc.dir);
copyfile('octave/pdepe.m', desc.dir);
end