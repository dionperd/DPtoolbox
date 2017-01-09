function R = rotation(N,phi)

%It rotates an nD matrix applying N*(N-1)/2 subsequent rotations in all of
%its 2D subspaces by angles that are given in phi.

R = eye(N);

ii=0;

%For every possible 2D subspace:
for i1 = 1:N-1;
    
    for i2 = i1+1:N;
        
        ii=ii+1; %increase rotation counter
        
        thisR = eye(N); %initialize this rotation matrix as the unit matrix
        
        %calculate cosine and sine of the angle
        c = cosd(phi(ii));
        s = sind(phi(ii));
        
        %create the rotation matrix
        thisR(i1,i1) = c;
        thisR(i2,i2) = c;
        thisR(i1,i2) = s;
        thisR(i2,i1) = -s;
        
        %multiply with R in order to add this rotation
        R = R*thisR;
        
    end
end

