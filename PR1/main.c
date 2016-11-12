//  Created by Ben McMurtry on 04/11/2016.

/*This program was written to solve the gravitational three-body problem.
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
 • Simulate the Sun-Earth-Moon system.*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h> /* for PATH_MAX */
#include <math.h>
#include <errno.h>
#include "myprog.h"

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
#define GNUPLOT_DATA   "myprog_gnuplot.dat"

/*Definitions for Verlet method*/
#define GRAV 6.67408e-11
#define TIMESTEP 86400
#define N_POINTS 100

/*Creates coordinate system*/
typedef enum Coords { x, y, z } coords;
/*Create a struct called body, which has a name, mass, and 3 dimensions of coordinates, velocity and acceleration*/
typedef struct body {
    char name[MAX_NAME_SIZE];
    double mass;
    double r[sizeof(coords)]; /* displacement */
    double v[sizeof(coords)]; /* velocity */
    double a[sizeof(coords)]; /* acceleration */
} Body;

/*Finds the number of vald bodies in INPUT_FILE, so that it can be used to initialise the array of bodies*/
static int Findbodycount() {
    char line[MAX_FILE_LINE_SIZE];
    FILE *input = fopen( INPUT_FILE, "r" );
    if (!input) {
        fprintf(stderr, "Error: Could not open file '%s'.\n", INPUT_FILE);
        exit(1);
    }
    int nBodiesScanned = 0;
    while (fgets(line, MAX_FILE_LINE_SIZE, input) ) {
        if (line[0] != '#') {
            nBodiesScanned += 1;
        }
    }
    return nBodiesScanned;
}

/*This functions reads the body data from the file INPUT_FILE, and puts the data into a number of Body structs equal to the number of lines in the file that don't start with '#', then prints the body data to stdout*/
Body * makebodies(nBodiesScanned) {
    char line[MAX_FILE_LINE_SIZE];
    char nameBuf[MAX_FILE_LINE_SIZE];
    FILE *input = fopen( INPUT_FILE, "r" );
    if (!input) {
        fprintf(stderr, "Error: Could not open file '%s'.\n", INPUT_FILE);
        exit(1);
    }
    
    printf("Number of bodies is %d, and they have:\n", nBodiesScanned);
    printf("%-12s%-12s%-12s%-12s%-12s%-12s%-12s%-12s\n", "Name", "Mass", "r.x", "r.y", "r.z", "v.x", "v.y", "v.z");
    
    Body *bodies = malloc(nBodiesScanned * sizeof(Body));
    rewind(input);
    
    int bodyN = 0; /* Number of bodies successfully found */
    while (fgets(line, MAX_FILE_LINE_SIZE, input) ) {
        if (line[0] != '#') {
            int nItemsScanned = sscanf(line,"%s %lg %lg %lg %lg %lg %lg %lg",
                                       nameBuf, &bodies[bodyN].mass,
                                       &bodies[bodyN].r[x], &bodies[bodyN].r[y], &bodies[bodyN].r[z],
                                       &bodies[bodyN].v[x], &bodies[bodyN].v[y], &bodies[bodyN].v[z]);
            if (nItemsScanned == ITEMS_PER_LINE) {
                strncpy(bodies[bodyN].name,nameBuf,MAX_NAME_SIZE);
                bodyN++;
            } else {
                fprintf(stderr, "Unknown format: %s\n",line);
            }
        }
    }
    fclose(input);
    for (bodyN = 0; bodyN < nBodiesScanned; bodyN++) {
        printf("%-12s%-12lg%-12lg%-12lg%-12lg%-12lg%-12lg%-12lg\n",
               bodies[bodyN].name, bodies[bodyN].mass,
               bodies[bodyN].r[x], bodies[bodyN].r[y], bodies[bodyN].r[z],
               bodies[bodyN].v[x], bodies[bodyN].v[y], bodies[bodyN].v[z]);
    }
    return bodies;
}

double separation(Body * bodies, int N, int M) {
    
    return separation;
}


void calcacc(Body *bodies, int nBodiesScanned) {
    int bodyN, bodyM, coord;
    double separation;
    for(bodyN = 0; bodyN < nBodiesScanned; bodyN++) {
        bodies[bodyN].a[0] = 0.0;
        bodies[bodyN].a[1] = 0.0;
        bodies[bodyN].a[2] = 0.0;
        for (bodyM = 0; bodyM < nBodiesScanned; bodyM++) {
            if (bodyN != bodyM) {
                separation = fabs(sqrt(pow((bodies[bodyN].r[x] - bodies[bodyM].r[x]), 2) - pow((bodies[bodyN].r[y] - bodies[bodyM].r[y]), 2) - pow((bodies[bodyN].r[z] - bodies[bodyM].r[z]), 2)));
                printf("Distance between %s and %s is %lg\n", bodies[bodyN].name, bodies[bodyM].name, separation);
                for (coord = 0; coord < sizeof(coords) - 1; coord++) {
                    bodies[bodyN].a[coord] += (-GRAV * bodies[bodyM].mass * (bodies[bodyN].r[coord] - bodies[bodyM].r[coord])) / pow(separation, 3);
                    printf("acceleration of body %s in direction %d is %lg\n", bodies[bodyN].name, coord, bodies[bodyN].a[coord]);
                }
            }
        }
    }
}

    /*for (bodyN = 0; bodyN < nBodiesScanned; bodyN++) {
        for (bodyM = 0; bodyM < nBodiesScanned; bodyM++) {
            if (bodyM != bodyN) {
        bodies[bodyN].a[x] += ((bodies[bodyN].r[x] - bodies[bodyM].r[x]) * -GRAV * bodies[bodyM].mass) / pow(separation[bodyN][bodyM], 3);
        bodies[bodyN].a[y] += ((bodies[bodyN].r[y] - bodies[bodyM].r[y]) * -GRAV * bodies[bodyM].mass) / pow(separation[bodyN][bodyM], 3);
        bodies[bodyN].a[z] += ((bodies[bodyN].r[z] - bodies[bodyM].r[z]) * -GRAV * bodies[bodyM].mass) / pow(separation[bodyN][bodyM], 3);
            }
        }
    }
    for (bodyN = 0; bodyN < nBodiesScanned; bodyN++) {
    printf("Acceleration of %s in x direction is %lg\n", bodies[bodyN].name, bodies[bodyN].a[x]);
    printf("Acceleration of %s in y direction is %lg\n", bodies[bodyN].name, bodies[bodyN].a[y]);
    printf("Acceleration of %s in z direction is %lg\n", bodies[bodyN].name, bodies[bodyN].a[z]);
    }
    
    
    bodies[bodyN].r[x] = bodies[bodyN].r[x] + bodies[bodyN].v[x] * TIMESTEP + ((bodies[bodyN].a[x] * TIMESTEP * TIMESTEP) / 2.0);
    //bodies[bodyN].v[x] = bodies[bodyN].v[x] +
    
    //fclose(input);
//}*/

/*
 ** This is a simple method to get gnuplot to plot data from within a C program.
 ** It should create 'myprog_gnuplot.dat' and then use the gnuplot script:
 **
 **    ./myprog_gnuplot.script
 **
 ** to plot a figure-of-eight.
 */
static int plotdata() {
    char command[PATH_MAX];
    
    printf("%s v%s - Demonstration of plotting data with gnuplot.\n", PROGNAME, VERSION );//this + defines not needed
    
    if (!write_figure_of_eight( GNUPLOT_DATA ) ) {
        /*
         ** Construct the command we will use to invoke gnuplot. Use snprintf() to
         ** show how to use fixed length buffers safely in case PATH_MAX isn't enough
         */
        snprintf(command, sizeof(command), "%s %s", GNUPLOT_EXE, GNUPLOT_SCRIPT );
        system( command );
    } else {
        fprintf(stderr, "Could not write to file: %s\n", GNUPLOT_DATA );
        return (1);
    }
    return(0);
}

static int write_figure_of_eight( const char *filename ) { //write body positions at times
    FILE *fdat = fopen (filename, "w");
    double d;
    if (!fdat) {
        /* fopen() sets the global variable errno when it fails */
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
}

int main () {
    
    int nBodiesScanned = Findbodycount();
    
    Body *bodies;
    bodies = makebodies(nBodiesScanned);
    
    calcacc(bodies, nBodiesScanned);
    //plotdata();
    free (bodies);  /* XXX use xfree() */
}
