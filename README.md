# kde2d-free-energy
Weighted 2D Gaussian kernel density estimator to turn scattered, weighted samples (x y w) into a smooth free-energy map F(x,y) = −kBT ln p(x,y) on a regular grid. kernel.sh sets up the grid, chooses the kernel widths by Silverman's rule and writes the map (fe2d.txt).
