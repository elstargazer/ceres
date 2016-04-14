function new_complete_path = FillConfigTemplate(config_template_filename,mesh_filename,append,runname)

[path,name,ext] = fileparts(config_template_filename);

mkdir([path '/' runname]);
new_complete_path = [path '/' runname '/' name '_' append ext];

in_old = fopen(config_template_filename,'r');
in_new = fopen(new_complete_path,'w');

str = fgetl(in_old);
fprintf(in_new,'%s',str);

key_mesh = 'mesh_filename';
key_out  = 'output_folder';

while (str~=-1)
    str = fgetl(in_old);
    if strncmpi(str,key_mesh,numel(key_mesh))
        fprintf(in_new,[key_mesh ' = "' mesh_filename(1:end) '";\n']);
    elseif strncmpi(str,key_out,numel(key_out))
        fprintf(in_new,[key_out ' = "output/' runname '/output_' append '/";\n']);
        mkdir(['../output/' runname '/output_' append])
    else
        fprintf(in_new,'%s\n',str);
        
        
    end
end

fclose(in_old);
fclose(in_new);

