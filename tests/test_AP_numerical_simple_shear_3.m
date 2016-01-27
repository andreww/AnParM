function [Xs, Ls, ts] = test_AP_numerical_simple_shear_3

    % time step
    dt = 0.02165 ;

    [Xs, Ls, ts] = AP_path_integrate( [0.0 1.0 0.0]', 0.0, dt*21,...
        dt, 0.1, @AP_velocity_simple_shear);
    
    Ls = Ls(:,:,2:22);
    size(Ls)
    Ls
    
    [texture, time] = AP_calculate_texture(2000, 3.5, ...
        [0.3333 0.6667 1.0], Ls, dt, 21);
    
    time
    % Rotate for plotting
    [texture]=AP_rotate_texture_Euler(texture(:,:,22),0, 90.0, 0) ;
    
    % output the final texture    
    MVT_write_VPSC_file('test_3.out', ...
              texture, 'Simple shear output') ;
      
    % plot with MTEX
    MVT_olivine_pole_from_vpsc('test_3.out','scale',[0 7], ...
          'writefile','test_3','png')
end