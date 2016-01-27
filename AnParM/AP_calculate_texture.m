% AP_calculate_texture - Use the AnPAR model to evaluate texture evolution
%                        given material parameters and deformation history.
%
% // Part of AnParM - A MATLAB toolkit for the fast analytical modelling of  //
% //                         crystal deformation                             //
%  
%  Usage:
%
%   [texture, time] = AP_calculate_texture(ngrains,rn,tau,vgrad,dt,nsteps)
%        Input parameters:
%           ngrains : Number of grains in the calculation.
%                rn : Power law coefficient to apply (usually 3.5).
%               tau : CRSSs of the slip systems. This is a 1 x N vector where
%                     the length dictates the number of systems to be updated
%             vgrad : Either (3 x 3) velocity gradient tensor representing the
%                     deformation to be applied at each time step or a 
%                     (3 x 3 x nsteps) list of tensors, one to be applied 
%                     at each timestep. These tensors should be on a
%                     consistent, external, frame of reference. The texture
%                     is returned on this reference frame.
%                dt : the timestep for each strain increment.
%            nsteps : number of timesteps
%                 
%
%        Output parameters:
%           texture : tAn 3 x ngrains x nsteps+1 list of Euler angles in
%                     degrees, representing the orientation of the 
%                     crystals at each timestep.
%              time : time of each texture sample.
%
%
% Reference.
% ~~~~~~~~~~
%
%  Goulding, NG, Ribe, NM, Castelnau, O., Walker A. and Wookey, J. Analytical
%     parameterization of self-consistent polycrystal mechanics: Fast calculation 
%     of upper mantle anisotropy. Geophys. J. Int., in press.
%
% See also: 

% Copyright (c) 2016, Neil Goulding, Neil Ribe, Olivier Castelnau, 
%                     Andrew Walker and James Wookey
%
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright 
%      notice, this list of conditions and the following disclaimer in the 
%      documentation and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names of its 
%      contributors may be used to endorse or promote products derived from 
%      this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function [texture, time] = AP_calculate_texture(ngrains,rn,tau,vgrad,dt,...
                               nsteps)

    % Input checks for vgrad
    vgrad_size = size(vgrad);
    assert(vgrad_size(1) == 3, 'AP:argChk', ... 
        'vgrad must be 3 x 3 or 3 x 3 x nsteps. Rank 1 wrong');
    assert(vgrad_size(2) == 3, 'AP:argChk', ... 
        'vgrad must be 3 x 3 or 3 x 3 x nsteps. Rank 2 wrong.');
    if length(vgrad_size) == 2
        vgrad_it = 0;
    elseif length(vgrad_size) == 3
        assert(vgrad_size(3) == nsteps, 'AP:argChk', ... 
        'vgrad must be 3 x 3 x nsteps if it varies with time')
        vgrad_it = 1;
    else
        error('AP:argChk', ... 
        'vgrad must be 3 x 3 or 3 x 3 x nsteps');
    end


    % Empty arrays for output
    texture = zeros(3, ngrains, nsteps+1);
    time = zeros(nsteps+1, 1);

    % Set up initial FSE and texture
    FSE = eye(3);
    texture_now = MVT_make_random_texture( ngrains ) ;
    size(texture_now)
    
    % Store the initial texture and time
    time(1) = 0.0;
    texture(:,:, 1) = texture_now;
    
    % Loop over timesteps updating the FSE and texture at each point,
    % remembering that we need to update the texture on using vgrad on the
    % frame of the FSE. 
    for istep = 1:nsteps

        if vgrad_it
            vgrad_step = vgrad(:,:,istep);
        else
            vgrad_step = vgrad;
        end
        
        % Update the FSE for this strain incrememnt
        [FSE] = update_finitestrain(FSE,vgrad_step,dt); 
        
        % To update the texture at the current point 
        % for this we need the lengths of the FSE
        % and the rotation onto the reference frame
        % - i.e. the eigenvalues and vectors respctivly.
        [c, R] = FSEprincipAxes(FSE);
      
        % Put the Velocity gradient tensor onto the frame of reference.
        % of the FSE
        vgradR = R'*vgrad_step*R;
      
        % calculate log ratios of the FSE axis lenghts
        r12 = log(c(1)/c(2)) ;
        r23 = log(c(2)/c(3)) ;
        r13 = log(c(1)/c(3)) ; 
      
        % Update the texture here...

        [texture_now] = ...
            AP_Ol_texture_update(texture_now,rn,tau,vgradR,r12,r23,r13,dt);
      
        % Store the texture and time
        time(istep+1) = 0.0;
        texture(:,:, istep+1) = texture_now;
      
    end
    
end


function [c, R] = FSEprincipAxes(F)
    % Return the lenghts of the principal 
    % axes, c, and the orentations, R, of the 
    % finite strain ellipsoid assoceated 
    % with the deformation gradient tensor, F.
    
    % We need the Eigenvectors and Eigenvalues
    % of the left stretch tensor, V, which is 
    % derived from the Left Cauchy-Green deformation
    % tensor. And we need these to be sorted.

    % Left Cauchy-Green deformation tensor:
    B = F*F';
    
    % Now, B = V^2, so we can could take the
    % (matrix) square root sqrtm before finding the
    % eigenvalues and vectors, but it is probably 
    % quicker to find the eigenvectors and eigenvalues
    % of B. These are also the eigenvectors of V and the
    % eignevalues of V are just the (element wise) square 
    % root of the eigenvalues of B.
    
    [EIVEC, EIVAL] = eig(B);
    c_RAW = [sqrt(EIVAL(1,1)) sqrt(EIVAL(2,2)) sqrt(EIVAL(3,3))] ;
    [c, IND] = sort(c_RAW,2, 'descend') ;
    R = EIVEC ; % for dimensioning
    for i=1:3
        R(:,i) = EIVEC(:,IND(i)) ;
    end

end

function [FSE] = update_finitestrain(FSE,vgrad,t) 
   
   % McKenzie '79 eq 8-10
   A = eye(3) - (0.5*t)*vgrad;
   B = eye(3) + (0.5*t)*vgrad;
   FSE = inv(A)*B*FSE;
   
end   