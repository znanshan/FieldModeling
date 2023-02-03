%== get initial wall time 
time0=clock(); 
format long; 

out2 = fopen('area_frac.out','w'); 

%-- Simulation cell parameters: 

Nx = 64; 
Ny = 64; 
NxNy = Nx*Ny;
dx = 0.5; 
dy = 0.5; 

%--- Time integration parameters: 

nstep= 100; 
nprint= 10; 
dtime = 0.005; 
ttime= 0.0; 

%--- Material Parameters 

mobil = 5.0; 
grcoef = 0.1; 

%
%--- Generate initial grain_tructure 
%

iflag = 2; 
isolve = 1; 

[etas, ngrain, glist] = init_grain_micro (Nx, Ny, dx, dy, iflag, isolve); 

%
%-- Evolve: 
%

for istep=1:nstep

    ttime=ttime+dtime; 

    for igrain=1:ngrain 

        if (glist (igrain) ==1) 

            for i=1: Nx 
                for j=1: Ny 

                    eta(i,j) =etas(i,j,igrain); 

                end 
            end 

            for i=1:Nx 
                for j=1:Ny 

                    jp=j+1; 
                    jm=j-1; 

                    ip=i+1; 
                    im=i-1; 

                    jp=j+1; 
                    jm=j-1; 

                    ip=i+1; 
                    im=i-1; 

                    if(im == 0) 
                        im = Nx; 
                    end 
                    if (ip == (Nx+1)) 
                        ip=1; 
                    end 

                    if (jm == 0)
                        jm = Ny; 
                    end 

                    if (jp == (Ny+1)) 
                        jp=1; 
                    end 

                    hne=eta(ip, j);
                    hnw=eta(im, j); 
                    hns=eta(i, jm); 
                    hnn=eta(i, jp); 
                    hnc=eta(i, j); 

                    lap_eta(i,j) = (hnw + hne + hns + hnn -4.0*hnc)/(dx*dy); 
% 
%--Derivative of free energy 
% 
                    dfdeta=free_energ_fd_ca_vl (i, j, ngrain, etas, eta, igrain); 

%
%---Time integration: 
%
                    eta(i,j) = eta(i,j) - dtime*mobil* (dfdeta - grcoef*lap_eta(i,j)); 
%%%% 修改
%-- for small deviations: 
                    if (eta (i, j) >= 0.9999)
                        eta (i, j) = 0.9999; 
                    end 

                    if (eta (i, j) < 0.00001)
                        eta (i, j) = 0.00001; 
                    end 
                end %j 
            end %i 

%--

            grain_sum = 0.0; 
            for i=1:Nx 
                for j=1:Ny 
                    etas (i, j, igrain) =eta (i, j); 
                    grain_sum = grain_sum + eta (i, j); 
                end 
            end 

%-- Check volume fraction of current grain: 

            grain_sum = grain_sum/NxNy; 

            if (grain_sum <=0.001)
                glist(igrain)=0; 
                fprintf('grain: %5d is eliminated\n', igrain); 
            end 
%---- 
        end%if 
    end%Bigrain 

    if (mod(istep, nprint) == 0) 

        fprintf('done step: %5d\n', istep); 

%fnamel= sprintf('time_%d.out', istep);  
%out1 = fopen(fnamel, 'w'); 

%for i=1: Nx 
%for j=1: Ny 

%gg=0.0;
%for igrain=1:ngrain 
%gg=gg +etas (i, j, igrain)^2;
%end 

%fprintf (out1,'%5d %5d %14.6e\n',i,j,gg); 
%end
%end

%fclose(out1); 

%--- write vtk file & calculate area fraction of grains: 

        eta2=zeros(Nx, Ny); 
    
        fprintf(out2, '%14.6e', ttime); 
    
        for igrain=1:ngrain 
            ncount=0;        
            for i=1:Nx 
                for j=1:Ny         
     
                    eta2 (i, j) = eta2 (i, j) + etas (i, j, igrain) ^2; 

                    if(etas (i, j, igrain) >=0.5) 
                        ncount=ncount+1; 
                    end 

                end
            end             
            ncount=ncount/ (NxNy); 
            fprintf(out2, '%14.6e',ncount); 
        end 
        fprintf(out2, '\n'); 

%--

        write_vtk_grid_values (Nx, Ny, dx, dy, istep, eta2); 

    end %end if 
end %istep 

%--- calculate compute time: 
compute_time = etime (clock(), time0); 
fprintf('Compute Time: %10d\n', compute_time); 