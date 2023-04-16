Lx = 2*pi;
Ly = pi;
Lz = 1;

nx = 256;
ny = 256;
nz = 128;

x_coord = linspace(0, Lx, nx);
y_coord = linspace(0, Ly, ny);
z_coord = linspace(0, Lz, nz);


src_x = 0.3*Lx;
src_y = 0.5*Ly;
src_z = 0.5*Lz;
sigma = 1e-2;

for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            theta(i,j,k) = exp(-( (x_coord(i)-src_x)^2 + (y_coord(j)-src_y)^2 + (z_coord(k)-src_z)^2)/ (2*sigma^2));
        end
    end
end

% contour(z_coord, y_coord, squeeze(theta(72, :, :)))

fileID = fopen('/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/inputs/256_test/theta.00000000','w');
fwrite(fileID,theta,'double');
fclose(fileID);
