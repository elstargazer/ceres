function FillConfigTemplate(config_template_filename,mesh_filename,append)

[path,name,ext] = fileparts(config_template_filename);

in_old = fopen(config_template_filename,'r');
in_new = fopen([path '/' name '_' append ext],'w');

str = fgetl(in_old);
fprintf(in_new,'%s',str);

while (str~=-1)
    str = fgetl(in_old); 
    key = 'mesh_filename';
    if strncmpi(str,key',numel(key))
        
        fprintf(in_new,[key ' = "' mesh_filename '";\n']);
    else        
        fprintf(in_new,'%s\n',str);
    end   
end


fclose(in_old);
fclose(in_new);