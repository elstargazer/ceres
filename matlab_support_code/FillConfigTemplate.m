function FillConfigTemplate(config_template_filename,mesh_filename,append)

[path,name,ext] = fileparts(config_template_filename);

in_old = fopen(config_template_filename,'r');
in_new = fopen([path '/' name '_' append ext],'w');

str = fgetl(in_old);
fprintf(in_new,'%s',str);

while (str~=-1)
    str = fgetl(in_old); 
    key_mesh = 'mesh_filename';
    key_out  = 'output_folder';
    if strncmpi(str,key_mesh,numel(key_mesh))    
        fprintf(in_new,[key_mesh ' = "' mesh_filename(4:end) '";\n']);
    elseif strncmpi(str,key_out,numel(key_out))
            fprintf(in_new,[key_out ' = "output/output_' append '/";\n']);
            mkdir(['../output/output_' append])

    else
        fprintf(in_new,'%s\n',str);
        
        
    end   
end


fclose(in_old);
fclose(in_new);