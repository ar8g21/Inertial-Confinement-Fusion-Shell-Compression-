#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include <math.h>
#include <complex.h> 

const double R0 = 0.0450;               // Radius of the liquid droplet (m)
const double R1 = 0.0300;               // Radius of second liquid droplet (m)
const double x0 = 0.10;	                // Position of the liquid droplet on the x-axis (m)
const double rho_l = 1000.0;	        // Liquid density (kg/m^3)
const double rho_g = 1.2;		// Gas density (kg/m^3)
const double mu_l = 1e-3;		// Liquid viscosity (Pa s)
const double mu_g = 1.5e-5;		// Gas viscosity (Pa s)
const double sigma = 0.0072;	        // Surface tension (N/m^2)
const double Lambda = (-0.916);         // Sink Strength (m^2/s)
const double var = (Lambda/(4*M_PI));   // Numerical constant used in finding velocity components
const double e = 0.90;                  // eccentricity

// Maximum error in velocity used for mesh refinement
const double uemax = 0.1;

// Size of the simulation box
const double boxSize = 0.20;

// We are setting default values, but these are typically set using the optional arguments:
int LEVEL = 6;		// Maximum level of refinement
double tMax = 0.004;   // Maximum time of the simulation 

// Defining some helper variables
scalar f0[];
scalar u0[];
double xvel;
double yvel;
scalar ux[];
scalar uy[];

// Boundary conditions temporarily set to default, should be checked.
u.n[left] = dirichlet(0);
u.t[left] = dirichlet(0);
p[left] = neumann(0);
f[left] = f0[];
u.n[right] = dirichlet(0);
p[right] = neumann(0);

// Initial event that preprocesses and is completed before the intmain() function
event init (i = 0){
	if (!restore( file = "restart")) {
		const double sa = 0.0900;                            // semi major axis (m)
		const double sb = (sa * sqrt(1-sq(e))) ;           // semi minor axis (m)
		const double A = (1/(sa*sa));                      // Factor that multiplies the x components of the ellipse
		const double B = (1/(sb*sb));                      // Factor that multiplies the y components of the ellipse 

		// Refine grid near the droplet 
		refine ( ((A * sq(x-x0)) + (B * sq(y)) - 1 < 0 
		&& level< LEVEL)
		&& (sq(x-x0) + sq(y) - sq(1.2*R1) < 0 && level < LEVEL));
		
		// Set volume fraction to 1 inside the radius R0
		fraction (f0, A * sq(x-x0) + B * sq(y) - 1 < 0
		&& sq(x-x0) + sq(y) - sq(R1) > 0);
		
		// The following two lines were used in the atomization example. Not sure if they 
		// are needed.
		f0.refine = f0.prolongation = fraction_refine;
		restriction ({f0});
		
		// Definining the x and y components of velocities respectively, as a function of 
		// source/sink strength and radius, then setting an initial velocity equal to U0 
		// everywhere where we have defined f0=1 and initialising volume fraction and 
		// velocity field, respectively 
		
		foreach() {
		
		  xvel = (1 * (x-x0) * var) / (pow(sq(x-x0)+sq(y), 1.5));
	          yvel = (1 * (y) * var) / (pow(sq(x-x0)+sq(y), 1.5));
		  
		  ux[] = f0[]*xvel;
	          uy[] = f0[]*yvel;
	          
	          }
	          
	        foreach() {
	          
	          f[] = f0[];
		  u.x[] = ux[];
		  u.y[] = uy[];

		}
	}
}

event binary_output (t += 0.1){
	// Volume fraction output, full resolution
	static FILE * f_bin = fopen("f.bin", "a");
	output_matrix (f, f_bin, n=(1 << LEVEL), linear = true);
}

event images (t += 0.000001){
	// output separate png images in a subfolder
	char* fileName; // Make variable for filename
	if(0 > asprintf(&fileName, "png/f_%03d.png", (int)(t*10000))) return 1; // Open file
	output_ppm (f, file=fileName, linear=true, n=1024, min=0, max=1); // Write image to file
	free(fileName); // Free the file
}

event dump_output (t += 0.1){
	dump ();
}

// Apply grid refinement every timestep
event adapt (i++) {
	adapt_wavelet ({f,u}, (double[]){0.001, uemax, uemax}, LEVEL);
}

event report (i += 5){
	printf ("i = %d t = %g\n", i, t);
}

event end (t = tMax){
	printf ("i = %d t = %g\n", i, t);
}

int main ( int argc, char * argv[]) {
	// Currently two optional arguments can be added:
	if (argc > 1)	
		// The first one is the maximum time for the simulation,
		tMax = atof(argv[1]);
	if (argc > 2)
		// and the second one is the maximum level of refinement.
		LEVEL = atoi(argv[2]);

	init_grid(256);
	size (boxSize);
	
	// liquid properties
	rho1 = rho_l;
	mu1 = mu_l;
	// gas properties
	rho2 = rho_g;
	mu2 = mu_g;
	// surface tension
	//f.sigma = sigma;

	// Run the simulation
	run();
}

