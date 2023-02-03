function  [dfdeta] = free_energ_fd_ca_v2(Nx,Ny,ngrain, etas,eta,igrain)
    format long; 
    A=1.0;
    B=1.0; 
    NxNy = Nx*Ny;
    sum=zeros(NxNy,1); 
    for jgrain=1:ngrain
        if(jgrain ~= igrain)
            sum = sum +etas(  :,jgrain) .^2;
        end
    end

    dfdeta = A*(2.0*B* eta.*sum+eta.^3 - eta);

end %endfunction
