%== get initial wall time
time0=clock();
format long;

out2 = fopen('area_frac.out','w');

% - - Simulation cell parameters: 16

Nx = 64;
Ny = 64;
NxNy= Nx*Ny;
dx = 0.5;
dy = 0.5;

% - - - Time integration parameters:

nstep =     5000;
nprint=      100;
dtime =    0.005;
ttime =      0.0;

% - - - Material Parameters

mobil =  5.0;
grcoef = 0.1;

%

% - - - Generate initial grain_structure
%
iflag  = 2;
isolve = 2;

[etas,ngrain,glist] = init_grain_micro(Nx,Ny,dx,dy, iflag,isolve);

%
% - - - Get Laplacian templet
%

[grad] =laplacian(Nx,Ny,dx,dy);

%
% Evolve:
%

eta = zeros(Nx*Ny,1);

for istep = 1:nstep

    ttime = ttime + dtime;

    for igrain = 1:ngrain

        if(glist(igrain) == 1)

            eta = etas(:, igrain);

            %
            % --Derivative of free energy
            %

    
            dfdeta = free_energ_fd_ca_v2(Nx,Ny, ngrain,etas,eta,igrain);

            %
            % - - - Time integration:
            %
            
            eta = eta - dtime*mobil*(dfdeta - grcoef*grad*eta);
            
            % - - for small deviations:
            
            inrange = (eta >= 0.9999);
            eta(inrange) = 0.9999;
            
            inrange = (eta < 0.00001);
            eta(inrange) = 0.00001;
            % - -
            etas(  :,igrain) =eta;
            % - - Check volume fraction of current grain:
            grain_sum = 0.0;
            grain_sum = sum(eta)/NxNy;
            if(grain_sum <= 0.001)
                glist(igrain) =0;                              
                fprintf('grain: %5d is eliminated \n',igrain);
            end
            % - - - -

        end %if
    end %igrain

    if(mod(istep,nprint)==0)
        fprintf('done step: %5d\n',istep);
        %fname1 =sprintf(’time_%d .out’, istep);
        %out1 = fopen(fname1,’w’);
        %for i=1:Nx
        %for j=1:Ny
        %ii = (i - 1)*Nx+j
        % gg= 0 . 0;
        %for igrain =1:ngrain
        %gg = gg +etas(ii,igrain)^2;
        %end
        %fprintf(out1,’%5d %5d %14 . 6e\n’, i,j,gg);
        %end
        %end
        %fclose(out1);
        
        % - - - write vtk file & calculate area fraction of grains:


        eta2=zeros(Nx,Ny); 
        fprintf(out2,'%14.6e',ttime); 
        for igrain=1:ngrain
            ncount=0;
            for i=1:Nx
                for j=1:Ny
                    ii = (i - 1)*Nx+j; 
                    eta2(i,j) =eta2(i,j)+etas(ii, igrain)^2;
 
                    if(etas(ii,igrain) >= 0.5)
                        ncount=ncount+1;
                    end
 
                end
            end
            ncount=ncount/(NxNy);
            fprintf(out2,'%14.6e',ncount);
        end
        fprintf(out2,'\n');
 
        % - -
 
        write_vtk_grid_values(Nx,Ny,dx, dy,istep,eta2);
 
    end %end if

end %istep
 
% - - - calculate compute time: 159
compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d\n', compute_time);



