//  Created by Ben McMurtry on 04/11/2016.

/*
 This program was written to solve the gravitational three-body problem.
 The program is designed so that:
 • Initial conditions are read from a text file INPUT_FILE, with one line for each body comprising
 values representing the name (a string), mass, three components of position and three
 components of velocity all separated by tabs.
 • The number of bodies is not limited to three and is determined by the number of valid lines
 found in the data file.
 • The user can specify, preferably at run-time so recompilation is not required, whether a
 GNU Scientific Library (or a similar library) routine or your own version of the Verlet
 method is used for the solution.
 • The total energy (potential + kinetic) of the bodies is calculated and used to quantify the
 accuracy of simulation.
 • It uses a Gnuplot to plot orbits in the manner illustrated in the lecture.
 Use your program to:
 • Confirm and plot at least one of the known stable orbit solutions for the three-body problem.
 • Simulate the Sun-Earth-Moon system.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h> /* for PATH_MAX */
#include <math.h>
#include <errno.h>
//#include "myprog.h"

//includes for GSl
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*Definitions for the text file containing the bodies INPUT_FILE*/
#define INPUT_FILE "/Users/benmcmurtry/Desktop/Physics/Year4/ComputationalPhysicsandModelling/PR1/PR1/OrbitData.txt"
#define MAX_FILE_LINE_SIZE 250
#define ITEMS_PER_LINE 8
#define MAX_NAME_SIZE 32

/*Definitions for plotdata function*/
#define PROGNAME "PR1"//delete?
#define VERSION	 "1.0.0"//delete?
#define GNUPLOT_EXE    "/opt/local/bin/gnuplot"
#define GNUPLOT_SCRIPT "/Users/benmcmurtry/Desktop/Physics/Year4/ComputationalPhysicsandModelling/PR1/PR1/myprog_gnuplot.script"
#define GNUPLOT_DATA   "/Users/benmcmurtry/Desktop/Physics/Year4/ComputationalPhysicsandModelling/PR1/PR1/Sun-Moon-Earth.txt"

/*Definitions for Verlet method*/
#define GRAV 6.67408e-11
#define TIMESTEP (60*60*24)
#define N_POINTS 3650

//Make everything static and all that other stuff he said to do. Also check for sun half loops

/*Creates coordinate system*/
typedef enum Coords { x, y, z } coords;

/*Create a struct called body, which has a name, mass, and has a position, velocity and acceleration in all Coords.*/
typedef struct body {
    char name[MAX_NAME_SIZE];
    double mass;
    double r[sizeof(coords)]; /* displacement */
    double v[sizeof(coords)]; /* velocity */
    double a[sizeof(coords)]; /* acceleration */
} Body;

/*Finds the number of valid bodies in INPUT_FILE, so that it can be used to initialise the array of bodies.*/
static int Findbodycount() {
    char line[MAX_FILE_LINE_SIZE];
    FILE *input = fopen( INPUT_FILE, "r" );
    if (!input) {
        fprintf(stderr, "Error: Could not open file '%s'.\n", INPUT_FILE);
        exit(1);
    }
    int Bodycount = 0;
    while (fgets(line, MAX_FILE_LINE_SIZE, input) ) {
        if (line[0] != '#') {
            Bodycount += 1;
        }
    }
    //rewind(input);
    return Bodycount;
}

/*This function initialises an array of (Bodycount) Body structs called bodies, and fills it with the data it reads from INPUT_FILE lines that don't start with '#'. It then prints the body data to stdout*/
Body * makebodies(Bodycount) {
    char line[MAX_FILE_LINE_SIZE];
    char nameBuf[MAX_FILE_LINE_SIZE];
    FILE *input = fopen(INPUT_FILE, "r");
    if (!input) {
        fprintf(stderr, "Error: Could not open file '%s'.\n", INPUT_FILE);
        exit(1); // dont need this here since its in first one.
    }
    printf("Number of bodies is %d, and they have the following properties:\n", Bodycount);
    printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n", "NAME", "MASS", "r.x", "r.y", "r.z", "v.x", "v.y", "v.z");
    
    Body *bodies = malloc(Bodycount * sizeof(Body));
    int bodyN = 0; /* Number of bodies successfully found */
    while (fgets(line, MAX_FILE_LINE_SIZE, input) ) {
        if (line[0] != '#') {
            int nItemsScanned = sscanf(line, "%s %lg %lg %lg %lg %lg %lg %lg",
                                       nameBuf, &bodies[bodyN].mass,
                                       &bodies[bodyN].r[x], &bodies[bodyN].r[y], &bodies[bodyN].r[z],
                                       &bodies[bodyN].v[x], &bodies[bodyN].v[y], &bodies[bodyN].v[z]);
            if (nItemsScanned == ITEMS_PER_LINE) {
                strncpy(bodies[bodyN].name,nameBuf,MAX_NAME_SIZE);
                bodyN++;
            } else {
                fprintf(stderr, "Unknown format: %s\n", line);
            }
        }
    }
    fclose(input);
    for (bodyN = 0; bodyN < Bodycount; bodyN++) {
        printf("%-12s%-12lg%-12lg%-12lg%-12lg%-12lg%-12lg%-12lg\n",
               bodies[bodyN].name, bodies[bodyN].mass,
               bodies[bodyN].r[x], bodies[bodyN].r[y], bodies[bodyN].r[z],
               bodies[bodyN].v[x], bodies[bodyN].v[y], bodies[bodyN].v[z]);
    }
    return bodies;
}

/*Calculates accelerations for all bodies due to all other bodies*/
void calcacc(Body *bodies, int Bodycount) {
    int bodyN, bodyM, coord;
    double separation;
    for(bodyN = 0; bodyN < Bodycount; bodyN++) {
        bodies[bodyN].a[0] = 0.0;
        bodies[bodyN].a[1] = 0.0;
        bodies[bodyN].a[2] = 0.0;
        for (bodyM = 0; bodyM < Bodycount; bodyM++) {
            if (bodyN != bodyM) {
                separation = sqrt(pow((bodies[bodyN].r[x] - bodies[bodyM].r[x]), 2) + pow((bodies[bodyN].r[y] - bodies[bodyM].r[y]), 2) + pow((bodies[bodyN].r[z] - bodies[bodyM].r[z]), 2));
                if (isnan(separation))
                    printf("error");
                //printf("Distance between %s and %s is %lg m\n", bodies[bodyN].name, bodies[bodyM].name, separation);
                for (coord = 0; coord < sizeof(coords) - 1; coord++) {
                    bodies[bodyN].a[coord] += (-GRAV * bodies[bodyM].mass * (bodies[bodyN].r[coord] - bodies[bodyM].r[coord])) / pow(separation, 3);
                }//can i/should i store a separation matrix?
            }
        }
        /*for (coord = 0; coord < sizeof(coords) - 1; coord++) {
            printf("acceleration of body %s in direction %d is %lg\n", bodies[bodyN].name, coord, bodies[bodyN].a[coord]);
        }*/
        //printf("acceleration of body %s in x direction is %lg m/s^2\n", bodies[bodyN].name, bodies[bodyN].a[x]);
        //printf("acceleration of body %s in y direction is %lg m/s^2\n", bodies[bodyN].name, bodies[bodyN].a[y]);
        //printf("acceleration of body %s in z direction is %lg m/s^2\n", bodies[bodyN].name, bodies[bodyN].a[z]);
    }
}

static int Write_VerletStep(const char *filename, Body *bodies, int Bodycount) { //write body positions at times
    FILE *fdat = fopen (filename, "w");
    if (!fdat) {
            //fopen() sets the global variable errno when it fails
            return errno;
    }
    fprintf(fdat, "# Created by %s v%s\n",  PROGNAME, VERSION );
    fprintf(fdat, "# Body 1 pos1, pos2, pos3 etc.n" );
    //calcacc(bodies, Bodycount);
    for (int d = 0; d < N_POINTS; d++) {
        int bodyN;
        for (bodyN = 0; bodyN < Bodycount; bodyN++) { //Update positions and half update velocities
            for (int coord = 0; coord < sizeof(coords) - 1; coord++) {
                bodies[bodyN].r[coord] += (bodies[bodyN].v[coord] * TIMESTEP) + ((bodies[bodyN].a[coord] * TIMESTEP * TIMESTEP) / 2.0);
                bodies[bodyN].v[coord] += ((bodies[bodyN].a[coord] * TIMESTEP) / 2.0);
            }
            /*bodies[bodyN].r[x] += (bodies[bodyN].v[x] * TIMESTEP) + ((bodies[bodyN].a[x] * TIMESTEP * TIMESTEP) / 2.0);
            bodies[bodyN].v[x] += ((bodies[bodyN].a[x] * TIMESTEP) / 2.0);
            bodies[bodyN].r[y] += (bodies[bodyN].v[y] * TIMESTEP) + ((bodies[bodyN].a[y] * TIMESTEP * TIMESTEP) / 2.0);
            bodies[bodyN].v[y] += ((bodies[bodyN].a[y] * TIMESTEP) / 2.0);
            bodies[bodyN].r[z] += (bodies[bodyN].v[z] * TIMESTEP) + ((bodies[bodyN].a[z] * TIMESTEP * TIMESTEP) / 2.0);
            bodies[bodyN].v[z] += ((bodies[bodyN].a[z] * TIMESTEP) / 2.0);*/
            fprintf(fdat, "%-16.8lg%-16.8lg%-16.8lg", bodies[bodyN].r[x], bodies[bodyN].r[y], bodies[bodyN].r[z]);
        }
        calcacc(bodies, Bodycount); //Compute new accelerations and half update velocities
        for (bodyN = 0; bodyN < Bodycount; bodyN++) {
            bodies[bodyN].v[x] += ((bodies[bodyN].a[x] * TIMESTEP) / 2.0);
            bodies[bodyN].v[y] += ((bodies[bodyN].a[y] * TIMESTEP) / 2.0);
            bodies[bodyN].v[z] += ((bodies[bodyN].a[z] * TIMESTEP) / 2.0);
            //fprintf(fdat, "%-16.8lg%-16.8lg%-16.8lg", bodies[bodyN].r[x], bodies[bodyN].r[y], bodies[bodyN].r[z]);
        }
        fprintf(fdat, "\n");
    }
    fclose(fdat);
    return 0;
}

/*int func (double t, const double y[], double f[], void *params) {
    (void)(t); // avoid unused parameter warning
    double mu = *(double *)params;
    f[0] = y[1];
    f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

int jac (double t, const double y[], double *dfdy, double dfdt[], void *params) {
    (void)(t); // avoid unused parameter warning
    double mu = *(double *)params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}*/

int main () {
    
    int Bodycount = Findbodycount(); //Returns integer number of bodies for makebodies
    Body *bodies = makebodies(Bodycount); //Creates an array of length Bodycount of Body structs, and assigns them the values from the INPUT_FILE
    calcacc(bodies, Bodycount); //Calculates the initial accelerations of all the bodies due to all other bodies, and stores them in bodies.
    Write_VerletStep(GNUPLOT_DATA, bodies, Bodycount);//Steps through Verlet method writing body coordinates into GNUPLOT_DATA
    //plotdata();
    free (bodies);  /* XXX use xfree() */
}

/*
 ** This is a simple method to get gnuplot to plot data from within a C program.
 ** It should create 'myprog_gnuplot.dat' and then use the gnuplot script:
 **
 **    ./myprog_gnuplot.script
 **
 ** to plot a figure-of-eight.
 */

/*static int plotdata() {
 char command[PATH_MAX];
 
 printf("%s v%s - Demonstration of plotting data with gnuplot.\n", PROGNAME, VERSION );//this + defines not needed
 
 if (!write_figure_of_eight( GNUPLOT_DATA ) ) {
 
 //Construct the command we will use to invoke gnuplot. Use snprintf() to
 //show how to use fixed length buffers safely in case PATH_MAX isn't enough
 
 snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
 system( command );
 } else {
 fprintf(stderr, "Could not write to file: %s\n", GNUPLOT_DATA );
 return (1);
 }
 return(0);
 }*/

/*static int write_figure_of_eight( const char *filename ) { //write body positions at times
 FILE *fdat = fopen (filename, "w");
 double d;
 if (!fdat) {
 //fopen() sets the global variable errno when it fails
 return errno;
 }
 fprintf(fdat, "# Created by %s v%s\n",  PROGNAME, VERSION );
 fprintf(fdat, "# Phase\tX\tY\n" );
 for ( d = 0.0; d < N_POINTS; d++ ) {
 double phase = ( 2.0 * M_PI * d ) / ( N_POINTS - 1 );
 fprintf(fdat, "%.8e\t%.8e\t%.8e\n", phase, sin(2.0*phase), cos(1.0*phase) );
 }
 fclose(fdat);
 return 0;
 }*/
