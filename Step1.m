%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 1                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Step1(basepath, pathtodata21, pathtodata, years, latlonbounds, yearmetadata, fullextract, slaoradt)
%% Initialize
addpath(strcat(basepath, 'FUNCTIONS'));
%% Set loop for desired number of years
for i=1:length(years)
    output_dir= strcat(basepath);
    if str2double(years(i)) == 2021
        input_str = pathtodata21;
        input_dir = convertStringsToChars(input_str);
    else
        direct = pathtodata;
        if exist(direct)==7 %7 means name is a folder
            input_str = ([pathtodata   '/']);%uigetdir(direct,'DIRECTORY OF AVISO SLA OR MADT FILES');
            input_dir = convertStringsToChars(input_str);
        else
            input_dir = uigetdir('','DIRECTORY OF AVISO SLA OR MADT FILES');
        end
        uvdir = input_dir;
    end
    months = ['01','02','02','04','05','06','07','08','09','10','11','12'];
    NL = latlonbounds(1);
    SL = latlonbounds(2);
    WL = latlonbounds(4);
    EL = latlonbounds(3);
    Begin = yearmetadata(1);
    End = yearmetadata(2);
    extrac = fullextract;
    yrb = fix(Begin/10000); datebeg = datenum(num2str(Begin),'yyyymmdd');
    yre = fix(End/10000); dateend = datenum(num2str(End),'yyyymmdd');
    if extrac ==1
        Extraction_Type = 'Basic Extraction';
    elseif extrac==2
        Extraction_Type = 'Full Extraction';
    end
    Directory_Of_NetCDF_Files = input_dir; Directory_Of_Extracted_Files = output_dir;
    Northern_Limit = NL; Southern_Limit = SL; Eastern_Limit = EL; Western_Limit = WL;
    Beginning_Date = Begin; Final_Date = End;
    Date_Of_Extraction = datestr(now);
    savedest = strcat(basepath, 'EXTRACTION/ConfigFile_fulltime.mat');
    save  (savedest, 'Date_Of_Extraction', 'Extraction_Type', 'Directory_Of_NetCDF_Files', ...
        'Directory_Of_Extracted_Files', 'Northern_Limit', 'Southern_Limit', 'Eastern_Limit', 'Western_Limit', 'Beginning_Date', 'Final_Date');
    
    %% Loop through an individual year at this point
    for yr =1:length(years)
        yrstr=num2str(years(yr));
        list_SSH = [dir([input_str yrstr  '/*.nc'])]; % make sure this directory exists
        list_UV = [dir([ input_str yrstr '/*.nc'])]; % this one too
        cd([input_str '/' yrstr '/'])
        for j = 1:length(list_SSH)
            % Display file treated name and create associated filename
            filename_SSH = list_SSH(j,:);
            filename_UV = list_UV(j,:);
            ADT=double(ncread(filename_SSH.name,'adt'));
            U=double(ncread(filename_UV.name,'ugos'));
            V=double(ncread(filename_UV.name,'vgos'));
            NbLongitudes=double(ncread(filename_SSH.name,'longitude'));
            NbLatitudes=double(ncread(filename_SSH.name,'latitude'));
            if any(NbLongitudes>180)
                NbLongitudes = NbLongitudes-360;
            end
            indlon = find(NbLongitudes>=WL & NbLongitudes<=EL);
            indlat = find(NbLatitudes>=SL & NbLatitudes<=NL);
            X = NbLongitudes(indlon);Y = NbLatitudes(indlat);
           
            f = repmat(coriolis(Y)',length(X),1);
            g = 9.81;
            % Create Matrices of distance between 2 points in x and y direction
            delta_x=NaN(length(X),length(Y));
            delta_y=NaN(length(X),length(Y));
            delta_x(:,1) = ac_distance(Y(1),X(3),Y(1),X(1))*1000;
            delta_x(:,end) =  ac_distance(Y(end),X(3),Y(end),X(1))*1000;
            dy = Y(2)-Y(1);
            delta_y(:,1) =  ac_distance(Y(1)+dy,X(1),Y(1)-dy,X(1))*1000;
            dy = Y(end)-Y(end-1);
            delta_y(:,1) =  ac_distance(Y(end)+dy,X(1),Y(end)-dy,X(1))*1000;
            for m=2:length(Y)-1
                delta_x(:,m) = ac_distance(Y(m), X(3),Y(m), X(1))*1000;
                delta_y(:,m) =  ac_distance(Y(m+1),X(1),Y(m-1),X(1))*1000;
            end
            % Reduce SSH to area of interest
            SSH = ADT(indlon,indlat);
            % compute geostrophic velocity
            U = U(indlon,indlat);
            V = V(indlon,indlat);
            
            %%
            % compute Derivative of veloticy
            
            % calcul de dU/dx
            dU = NaN*U;
            dU(2:end-1,:) = (U(3:end,:)-U(1:end-2,:));
            dX = repmat(delta_x(1,:),[size(U,1) 1]);
            Ux = dU./dX;
            
            %calcul de dV/dx
            dV = NaN*V;
            dV(2:end-1,:) = (V(3:end,:)-V(1:end-2,:));
            Vx = dV./dX;
            
            % calcul de dU/dy
            dU = NaN*U;
            dU(:,2:end-1) = (U(:,3:end)-U(:,1:end-2));
            dy = (Y(3)-Y(1))*ac_distance(0,0,1,0)*1000;
            dY = repmat(dy,[size(U,1) size(U,2)]);
            Uy = dU./dY;
            
            % calcul de dV/dy
            dV = NaN*V;
            dV(:,2:end-1) = (V(:,3:end)-V(:,1:end-2));
            Vy = dV./dY;
            
            % Vorticity
            Vorticity=Vx-Uy;
            
            % Shear, normal componant of strain and Okubo-Weiss parameter
            Ss=Vx+Uy;
            Sn=Ux-Vy;
            OW=Sn.^2+Ss.^2-Vorticity.^2;
            
            %Norm of speed and EKE
            Speed=sqrt(U.^2+V.^2);
            EKE=(U.^2+V.^2)./2;
            %
            % Take dtime values and store it in a datenum format
            date_num=datenum(filename_SSH.name(end-19:end-12),'yyyymmdd');
            % Store variables in a .mat file
            filename_out = [output_dir 'EXTRACTION/INDIVIDUAL_FILES/' 'ADT_uv_' filename_SSH.name(end-19:end-12)];
            save(filename_out,'X','Y','date_num','U','V','Vorticity','Ss','Sn','OW','Speed','EKE','SSH');
            
        end
    end
end
end


