function  [etas,ngrain,glist] = init_grain_micro(Nx,Ny,dx,dy,iflag,isolve)

    format long;

    if(iflag == 2)
       in = fopen('grain_25.inp','r');
    end

    % - - - - - - - - - - - - -
    % generate two grains
    % - - - - - - - - - - - - -
    
    if(iflag == 1)
        ngrain =2;
        % -etas(, 1)  first grain
        % -etas(, 2)  second grain
        x0 = Nx/2;
        y0 = Ny/2;
        radius = 14.0;
        grain
        
        for i=1:Nx
            for j=1:Ny
                ii= (i - 1)*Nx+j;
                if(isolve == 2)
                    etas(ii,1)=1.0;
                    etas(ii,2)=0.0;
                else
                    etas(i,j,1) =1.0;
                    etas(i,j,2) =0.0;
                end
                xlength =sqrt((i -x0)^2+(j -y0)^2); 
                if(xlength <= radius)
%%%%修改
                    if(isolve == 2)
                        etas(ii,1)=0.0;
                        etas(ii,2)=1.0;
                    else
                    etas(i,j,1)=0.0;
                    etas(i,j,2)=1.0;
                    end
                end %if
            end %j
        end %i
    end %iflag
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % generate polycrsytal microstructure %
    % 	
    if(iflag == 2)
        % --- -- -–
        % read the data generated from voroni_1 .m
        % --- -- -–
        twopi=8.0*atan(1.0);
        epsilon=1.0e-4;
        ndime=2;
        nvpoin = fscanf(in,'%d',1);
        nvnode = fscanf(in,'%d',1);
        nvelem = fscanf(in,'%d',1);
        ngrain = fscanf(in,'%d',1);
        for ipoin=1:nvpoin
            jpoin = fscanf(in,'%d',1);
            dummy = fscanf(in,'%lf %lf',[2,1]); 
            for idime=1:ndime
                vcord(jpoin,idime)=dummy(idime);
            end
        end
        
        for ielem=1:nvelem
            jelem= fscanf(in,'%d',1);
            dummy= fscanf(in,'%d',[nvnode+1,1]);
            for inode=1:nvnode+1
                vlnods(ielem,inode)=dummy(inode);
            end
        end
        %---------------------------
        for ielem=1:nvelem
            jnode=0;
            for inode=1:nvnode
                knode=vlnods(ielem,inode);
                if(knode ~= 0)
                    jnode = jnode +1;
                end
            end
            nnode2(ielem)=jnode;
        end
        %-------
        % form the grid
        %------
        for i=1:Nx
            gx(i) =i*dx;
        end
        for j=1:Ny
            gy(j)=j*dy;
        end
        % - - - - - - - - -
        % initilize order parameters
        % - - - - - - - - -
        for i=1:Nx
            for j=1:Ny
                for igrain=1:ngrain
        
                    if(isolve == 2)
                        ii= (i- 1)*Nx+j;
                        etas(ii,igrain)=0.0;
                    else
                        etas(i,j,igrain) = 0.0;
                    end
                end %igrain
            end %j
        end %i
        
        %--
        %--
        for i=1:Nx
            for j=1:Ny
        
                ii= (i- 1)*Nx+j;
        
                for ielem=1:nvelem
        
                    igrain=vlnods(ielem,nvnode+1);
        
                    theta =0.0;
        
                    mnode = nnode2(ielem);
        
                        for inode=1:mnode
        
                            knode=vlnods(ielem,inode);
                            xv1 =vcord(knode,1);
                            yv1 =vcord(knode,2);
        
                            jnode=vlnods(ielem,inode+1);
                            if(inode == mnode)
                                jnode =vlnods(ielem,1);
                            end
                            xv2=vcord(jnode,1);
                            yv2=vcord(jnode,2);
                            p1x =xv1-gx(i);
                            p1y =yv1-gy(j);
                            
                            p2x =xv2-gx(i);
                            p2y =yv2-gy(j);
                            
                            x1=sqrt(p1x*p1x+p1y*p1y);
                            x2=sqrt(p2x*p2x+p2y*p2y);

                            if(x1*x2 <= epsilon)
                                theta=twopi;
                            else
                                tx1= ((p1x*p2x+p1y*p2y)/(x1*x2)); 
                            end

                            if(abs(tx1) >=1.0)
                                tx1=0.9999999999;
                            end

                            theta=theta+acos(tx1);
                        end %inode

                        if(abs(theta-twopi) <= epsilon)

                        if(isolve == 2)
                            etas(ii,igrain) =1.0;
                        else
                            etas(i,j,igrain)=1.0;
                        end

                    end

                end % ielem

            end %i
        end %j
%------

    end %iflag

%- -- initialize  glist:

    for igrain=1:ngrain
    glist(igrain) =1.0;
    end

end %endfunction



