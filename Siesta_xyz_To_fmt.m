% for scanning trajectory from ab inito calculations
% Converting to xyz format

clear();
clc();
% %file name
% file_inp='admp_M11.log';
% %no of atoms
% atoms_total=45;

%file name
file_inp='C8_Dia_new.out';
%no of atoms
atoms_total=8;


%output file
file_out=sprintf('%s_fmt.xyz',file_inp);
fod=fopen(file_out,'w');

fid=fopen(file_inp,'r');  
count=0;

while 1
    tline = fgets(fid);
    if ~ischar(tline),   break,   end

%    disp(tline);
    count_check=0;
    if strfind(tline,'outcoor: Atomic coordinates (Ang):')
        fprintf(fod,'%d\n\n',atoms_total);
        for i=1:1:atoms_total
            
%             
%             atom_no=fscanf(fid,'%d',1);
            
            %scanning X
            strx=fscanf(fid,'%s',1);
%             strx=strrep(data,'D+','E+');
%             strx=strrep(strx,'D-','E-');
            x=str2double(strx);
            %scanning Y
            stry=fscanf(fid,'%s',1);
            y=str2double(stry);
            %scanning Z
            strz=fscanf(fid,'%s',1);
            z=str2double(strz);
            
            dummy=fscanf(fid,'%s',1);
            at_type=fscanf(fid,'%s',1);
            
            
            if at_type=='C'
                fprintf(fod,'C%8.3f%8.3f%8.3f\n',x,y,z);
            end
            dummy=fscanf(fid,'%s',1);
            count_check=count_check+1;
                        
        end
%         fprintf(fod,'\n');
        count=count+1

    end
%     display('step %d done',count);
    
end
fclose(fid);


fclose(fod);