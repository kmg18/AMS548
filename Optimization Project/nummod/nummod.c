/*********************************************************
  nummod.c
  -------------------
copyright : (C) 2006 by Ryan Brenke and Philip Yang Shen
email : rbrenke@bu.edu yangshen@bu.edu
 *********************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>

#include _MOL_INCLUDE_

#define PI 3.14159265

//Custom Functions
float calculate_dist_between_atoms(struct atom A, struct atom B);
void calculate_force(struct atomgrp* agA,struct atomgrp* agB, struct prms* prms, double* force_in_x, double* force_in_y, double* force_in_z );
void calculate_energy_atom(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms, double* Energy_atom );
double calculate_energy(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms);
void calculatePI();
void energy_function();
void mmc_two (struct atomgrp* agA, struct atomgrp* agB, struct prms* prms);

void print_short_args (char* app_name);
void print_help_args (char* app_name);
void print_args (char* app_name);
void my_rotate_atomgrp (struct atomgrp* prot, struct atomgrp* rotprot, struct rmatrix* rmatrix, struct tvector* center_of_rotation);
void perturb_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, double translate_rangei, double rotate_range, struct tvector* center_of_rotation);

#define MAXSLEN 200
char* ATOM_PRM_FILE;

int main (int argc, char* argv[])
{
	char* app_name = argv[0];
	if (argc < 2)
	{
		print_short_args (app_name);
		print_help_args (app_name);
		exit (EXIT_FAILURE);
	}

	char* rec_ifile; // input file
	char* lig_ifile; // input file
	//char* ofile = (char*) mymalloc (MAXSLEN * sizeof (char));
	ATOM_PRM_FILE = atom_prm_file (ATOM_PRM);


	size_t slen; // string length

	int c;
	while (1)
	{

		c = getopt (argc, argv, "hp:");
		if (c == -1)
			break;
		switch (c)
		{
			case 'h':
				print_args (app_name);
				return 0;
			case 'p':
				slen = strlen (optarg);
				if (slen > MAXSLEN)
				{
					fprintf (stderr, "atom parameter file name %s is too long\n", optarg);
					exit (EXIT_FAILURE);
				}
				ATOM_PRM_FILE = optarg;
				break;
			default:
				break;
		}
	}

	if (optind+1 < argc)
	{
		rec_ifile = argv[optind];
		optind++;
		lig_ifile = argv[optind];
		optind++;
		while (optind < argc)
		{
			printf ("ignored argument: %s\n", argv[optind]);
			optind++;
		}
	}
	else
	{
		print_short_args (app_name);
		print_help_args (app_name);
		exit (EXIT_FAILURE);
	}


	struct prms* prms = read_prms (ATOM_PRM_FILE, _MOL_VERSION_);
	//Read receptor and ligand structure
	struct atomgrp* agA = read_file_atomgrp (rec_ifile, prms);
	struct atomgrp* agB = read_file_atomgrp (lig_ifile, prms);

	
	printf("\n\n Solution for assignment3 Problem c\n");
	mmc_two(agA, agB, prms);

	//Mukesh - assignment1_problem_a
    struct tvector* cmA=center_of_mass(agA);
    struct tvector* cmB=center_of_mass(agB);
    struct tvector* cmAB1 = (struct tvector*) mymalloc (sizeof (struct tvector));
    cmAB1->X = cmA->X - cmB->X;
    cmAB1->Y = cmA->Y - cmB->Y;
    cmAB1->Z = cmA->Z - cmB->Z;

    for (int i = 1; i <= 10; i ++){

	    struct atomgrp* cmA_moved = copy_atomgrp(agA);
	    cmAB1->X = 0.5 + cmAB1->X;
	    cmAB1->Y = 0.5 + cmAB1->Y;
	    cmAB1->Z = 0.5 + cmAB1->Z;
	    translate_atomgrp (cmA_moved, cmA_moved,(cmAB1)); // translate agA
	    char* current_ofile = (char*) mymalloc (100 * sizeof (char));
	    sprintf (current_ofile, "a1p1%d.ms", i);
	    fprint_file_atomgrp(current_ofile, cmA_moved, prms);

    }

	//Mukesh - assignment1_problem_b
	printf("\n\n Solution for assignment1 Problem b\n");	
	double energy = 0;
	for(int i = 0; i< agA->natoms; i++) {	
		int atom_type_A = agA->atoms[i].atom_typen;
		double chrg_A = prms->chrgs[atom_type_A];

		for(int j = 0; j < agB->natoms; j++) {
			int atom_type_B = agB->atoms[j].atom_typen;
			double chrg_B = prms->chrgs[atom_type_B];
			float distance = calculate_dist_between_atoms(agA->atoms[i], agB->atoms[j]);
			energy = energy + (chrg_A * chrg_B)/distance;
		}
	}
	
	printf("Energy : %.3f\n", energy);


	//Mukesh - assignment1_problem_c
	printf("\n\n Solution for assignment1 Problem c\n");
	double* force_in_x = (double*)mymalloc(sizeof(double) * agA->natoms);
	double* force_in_y = (double*)mymalloc(sizeof(double) * agA->natoms);
	double* force_in_z = (double*)mymalloc(sizeof(double) * agA->natoms);
	double force_x = 0;
	double force_y = 0;
	double force_z = 0;

	calculate_force(agA, agB, prms, force_in_x, force_in_y, force_in_z);
	
	for(int i = 0; i< agA->natoms; i++ )
	{
		force_x = force_x + force_in_x[i];
		force_y = force_y + force_in_y[i];
		force_z = force_z + force_in_z[i];

	}
	printf("Force along X is = %f \n, Force along Y is = %f \n, Force along Z is = %f \n", force_x, force_y, force_z);


	//Mukesh - assignment1_problem_d
	printf("\n\n Solution for assignment1 Problem d\n");
	double* energy_array = (double*)mymalloc(sizeof(double) * agA->natoms);
	double* Energy_array_new = (double*)mymalloc(sizeof(double) * agA->natoms);
	calculate_energy_atom(agA, agB, prms, energy_array);
	struct tvector* cmAB = (struct tvector*) mymalloc (sizeof (struct tvector));
	cmAB->X = 0.01;
	cmAB->Y = 0;
	cmAB->Z = 0;
	struct atomgrp* agA_copy = copy_atomgrp(agA);
	translate_atomgrp (agA_copy, agA_copy, cmAB);
	calculate_energy_atom(agA_copy, agB, prms, Energy_array_new);
	
	double gradient_in_x = 0;
	double gradient_in_y = 0;
	double gradient_in_z = 0;

	for(int i = 0; i< agA->natoms; i++ )
	{
		gradient_in_x = gradient_in_x + ((energy_array[i] - Energy_array_new[i])/ 0.01);
	}
	cmAB->X = 0;
	cmAB->Y = 0.01;
	cmAB->Z = 0;
	agA_copy = copy_atomgrp(agA);
	translate_atomgrp (agA_copy, agA_copy, cmAB);
	calculate_energy_atom(agA_copy, agB, prms, Energy_array_new);
	
	for(int i = 0; i< agA->natoms; i++ )
	{
		gradient_in_y = gradient_in_y + ((energy_array[i] - Energy_array_new[i])/ 0.01);
	}

	cmAB->X = 0;
	cmAB->Y = 0;
	cmAB->Z = 0.01;
	agA_copy = copy_atomgrp(agA);
	translate_atomgrp (agA_copy, agA_copy, cmAB);
	calculate_energy_atom(agA_copy, agB, prms, Energy_array_new);
	
	for(int i = 0; i< agA->natoms; i++ )
	{
		gradient_in_z = gradient_in_z + ((energy_array[i] - Energy_array_new[i])/ 0.01);
	}
	printf("gradient = %f, %f, %f \n",gradient_in_x, gradient_in_y, gradient_in_z);


	//Mukesh - assignment1_problem_e
	printf("\n\n Solution for assignment1 Problem e\n");
	double Energy_shift[200];
	struct tvector* center_mass_a= center_of_mass(agA);
	struct tvector* center_mass_b= center_of_mass(agB);
	struct tvector* cmAB2 = (struct tvector*) mymalloc (sizeof (struct tvector));
	cmAB2->X = 0.01*(center_mass_a->X - center_mass_b->X);
	cmAB2->Y = 0.01*(center_mass_a->Y - center_mass_b->Y);
	cmAB2->Z = 0.01*(center_mass_a->Z - center_mass_b->Z);
	struct atomgrp* agA_copy1 = copy_atomgrp(agA);
	for(int i = 0; i< 200; i++ )
	{
		translate_atomgrp (agA_copy1, agA_copy1, cmAB2);
		Energy_shift[i] = calculate_energy(agA_copy1, agB, prms);
		printf(" %f \n",  Energy_shift[i]);
	}


	printf("\n\n Solution for assignment2 Problem a\n");
	calculatePI();
	printf("\n\n 1. Solution for assignment2 Problem b\n");
	energy_function();

	struct tvector* temp_tv = (struct tvector*) mymalloc (sizeof (struct tvector));
	temp_tv->X = 10;
	temp_tv->Y = 0;
	temp_tv->Z = 0;
	struct atomgrp* agA_moved = copy_atomgrp(agA);
	translate_atomgrp (agA_moved, agA_moved,temp_tv); // translate agA 
	char* current_ofile = (char*) mymalloc (100 * sizeof (char)); 
	sprintf (current_ofile, "test%d.ms", 1);
	fprint_file_atomgrp(current_ofile, agA_moved, prms);
	struct tvector* cm=center_of_mass(agA_moved);
	printf("Center of mass %.3f %.3f %.3f\n",cm->X,cm->Y,cm->Z);
	srand(time(NULL));
	double r=(double)(rand())/((RAND_MAX+1.0));
	printf("Random number %.4f\n",r);
	// Evaluate E
	clock_t start, end;
	double elapsed;
	float E=0;

	int i;
	start = clock();
	for (i = 0; i < 10; i++)
	{
		E = complex_energy (agA, agB, prms);
	}
	end = clock();
	elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf ("time: %.3f\n", elapsed);
	printf("Energy %.3f\n",E);
	return EXIT_SUCCESS;
}


void calculatePI() {
	int pass = 0;
	int i = 0;
	while (i < 15000) {
		double y_random = ((double) rand() / (RAND_MAX));
		double x_random = ((double) rand() / (RAND_MAX));
		
		double dist = sqrt(pow(x_random, 2) + pow(y_random, 2));
		if (dist <= 1) {
			pass += 1;
		}
		i += 1;
	}

	double pi_value = (4.0 * pass)/15000;
	printf("The PI value according to Monte Carlo Integration method is: %lf\n", pi_value);
}

double function_value(double input)
{
	double result = (-0.5 * pow(input,2) -0.5 * input - 0.3) * exp(-fabs(input)) + (0.01 * pow(input,2));
	return result;
}

double new_value()
{
	double result = (rand() % (10 + 10 + 1)) - 10;
	return result;

}

double increment(double a, double b) {
	return a+b;
}
void energy_function() {
	double  new_variable, a = 0;
	double  diffa;
	double E = function_value(a);
	double minimum;
	double least = E;
	double p;
	double E_new;
	for (int iter = 0; iter < 20000; iter += 1)
	{
		p = ((double) rand() / (RAND_MAX));
		diffa = new_value()/10;
		new_variable = increment(a, diffa);
		E_new = function_value(new_variable);
		if (E_new < E || p <= exp(E - E_new))
		{
			E = E_new;
			a = new_variable;
		}
		if(E < least)
		{
			least = E;
			minimum = a;
		}
	}
	printf("x value is : %f E value is: %f\n", minimum, least);
}


void mmc_two (struct atomgrp* agA, struct atomgrp* agB, struct prms* prms) {
	struct tvector* vector_new = (struct tvector* ) malloc (sizeof (struct tvector));

	double element0 = complex_energy(agA, agB, prms);
	double elementX = element0;
	for (int i = 0; i < 10000; i++) {
		struct atomgrp* copy_A = copy_atomgrp(agA);
		struct tvector* center_A = center_of_mass(agA);
		double p = ((double) rand() / (RAND_MAX));
		perturb_atomgrp(agA, copy_A, 2, PI/18, center_A);
		double complex_e = complex_energy(copy_A, agB, prms);
		if (complex_e < element0 || p <= exp(element0 - complex_e)) {
			element0 = complex_e;
			agA = copy_A;
			vector_new = center_A;
		}
		if (element0 < elementX) {
			elementX = element0;
		}
	}

	struct tvector* center_agB = center_of_mass(agB);
	printf("energy value is : %f\n", elementX);
	printf("Protein A center of mass is: (%f %f %f)\n", vector_new->X, vector_new->Y, vector_new->Z);
	printf("Protein B center of mass is: (%f %f %f)\n", center_agB->X, center_agB->Y, center_agB->Z);
}



void print_help_args (char* app_name)
{
	fprintf (stderr, "try '%s -h' for a list of arguments\n", app_name);
}

void print_short_args (char* app_name)
{
	fprintf (stderr, "usage: %s [arguments] RECEPTOR LIGAND\n", app_name);
	fprintf (stderr, "print correlation energy of RECEPTOR and LIGAND\n");
}

void print_args (char* app_name)
{
	print_short_args (app_name);

	printf ("\n");
	printf ("arguments:\n");
	printf ("   %-20s Use <atom.prm> as atom parameter file (default: %s)\n", "-p <atom.prm>", ATOM_PRM_FILE);
}

void my_rotate_atomgrp (struct atomgrp* prot, struct atomgrp* rotprot, struct rmatrix* rmatrix, struct tvector* center_of_rotation)
{
	int atomi;
	for (atomi = 0; atomi < prot->natoms; atomi++)
	{
		float X = prot->atoms[atomi].X;
		float Y = prot->atoms[atomi].Y;
		float Z = prot->atoms[atomi].Z;
		rotprot->atoms[atomi].X = rmatrix->a11*(X - center_of_rotation->X) + rmatrix->a12*(Y - center_of_rotation->Y) + rmatrix->a13*(Z - center_of_rotation->Z) + center_of_rotation->X;
		rotprot->atoms[atomi].Y = rmatrix->a21*(X - center_of_rotation->X) + rmatrix->a22*(Y - center_of_rotation->Y) + rmatrix->a23*(Z - center_of_rotation->Z) + center_of_rotation->Y;
		rotprot->atoms[atomi].Z = rmatrix->a31*(X - center_of_rotation->X) + rmatrix->a32*(Y - center_of_rotation->Y) + rmatrix->a33*(Z - center_of_rotation->Z) + center_of_rotation->Z;
	}
}

void perturb_atomgrp (struct atomgrp* ag, struct atomgrp* moved_ag, double translate_range, double rotate_range, struct tvector* center_of_rotation)
{
	struct tvector* temp_tv = (struct tvector*) mymalloc (sizeof (struct tvector));
	struct rmatrix* temp_rmatrix = (struct rmatrix*) mymalloc (sizeof (struct rmatrix));

	double qw,qx,qy,qz; //quaternions
	double r,theta,phi; //intermediate spherical coordinates for uniform sampling

	if(translate_range<0.0)
	{
		printf ("Input error: translational range should be NONNEGATIVE\n");
		exit (EXIT_FAILURE);
		//printf ("Notice: translational range is forced to be ZERO.\n");
		//translate_range=0;	    
	}
	if(rotate_range<0.0)
	{
		printf ("Input error: rotational range should be NONNEGATIVE\n");
		exit (EXIT_FAILURE);
		//printf ("Notice: rotational range is forced to be ZERO.\n");
		//rotate_range=0;
	}
	else
		if(rotate_range>PI)
		{
			printf ("Input error: maximum rotational range should be PI\n");
			exit (EXIT_FAILURE);
			//printf ("Notice: rotational range is forced to be PI.\n");	    
			//rotate_range=PI;
		}


	//translational perturbation
	/*Uniform sampling in a sphere of radius translate_range*/
	/* intermediate spherical coordinates (r,theta,phi) */

	//random number generator: to modify


	r = translate_range * pow((rand() / ((double)RAND_MAX + 1)),1/3.0);
	phi = acos(1-2*(rand() / ((double)RAND_MAX + 1)));
	theta = 2*PI*(rand() / ((double)RAND_MAX + 1));

	temp_tv->X = r * cos(theta) * sin(phi);
	temp_tv->Y = r * sin(theta) * sin(phi);
	temp_tv->Z = r * cos(phi);


	//rotational perturbation
	//Kuffner paper describes how to generate uniform unit quaternions
	//global uniform sampling: max range PI
	//to modify if need ``local'' orietational perturbation
	//uniform sampling in a sphere of exponential coordinates
	//essential: space of exponential coordinates is similar to the Euclidean space of translations

	r = rotate_range * pow((rand() / ((double)RAND_MAX + 1)),1/3.0);
	phi = acos(1-2*(rand() / ((double)RAND_MAX + 1)));
	theta = 2*PI*(rand() / ((double)RAND_MAX + 1));

	//transform into quaternions
	qw = cos(r/2);
	qx = sin(r/2) * cos(theta) * sin(phi);
	qy = sin(r/2) * sin(theta) * sin(phi);
	qz = sin(r/2) * cos(phi);

	//generate rotation matrix
	temp_rmatrix->a11 = 1 - 2 * (pow(qy,2) + pow(qz,2));
	temp_rmatrix->a12 = 2 * (qx*qy - qz*qw);
	temp_rmatrix->a13 = 2 * (qx*qz + qy*qw);

	temp_rmatrix->a21 = 2 * (qx*qy + qz*qw);
	temp_rmatrix->a22 = 1 - 2 * (pow(qx,2) + pow(qz,2));
	temp_rmatrix->a23 = 2 * (qy*qz - qx*qw);

	temp_rmatrix->a31 = 2 * (qx*qz - qy*qw);
	temp_rmatrix->a32 = 2 * (qy*qz + qx*qw);
	temp_rmatrix->a33 = 1 - 2 * (pow(qx,2) + pow(qy,2));

	my_rotate_atomgrp (ag, moved_ag, temp_rmatrix, center_of_rotation);
	translate_atomgrp (moved_ag, moved_ag, temp_tv);
}


float calculate_dist_between_atoms(struct atom A, struct atom B)
{
	float x_diff = A.X - B.X;
	float y_diff = A.Y - B.Y;
	float z_diff = A.Z - B.Z;
	float distance = pow(((x_diff * x_diff) + (y_diff * y_diff)+ (z_diff * z_diff)), 0.5);
	return distance;
}

void calculate_force(struct atomgrp* agA,struct atomgrp* agB, struct prms* prms, double* force_in_x, double* force_in_y, double* force_in_z ) {
	
	for (int i = 0; i < agA->natoms; i++) {	
		
		double forcex = 0;
		double forcey = 0;
		double forcez = 0;
		int atom_A = agA->atoms[i].atom_typen;
		double charge_A = prms->chrgs[atom_A];
		
		for (int j = 0; j < agB->natoms; j++) {
			
			int atom_B = agB->atoms[j].atom_typen;
			double charge_B = prms->chrgs[atom_B];
			float atom_distance = calculate_dist_between_atoms(agA->atoms[i], agB->atoms[j]);
			float distance_y = agA->atoms[i].Y - agB->atoms[j].Y;
			float distance_z = agA->atoms[i].Z - agB->atoms[j].Z;
			float distance_x = agA->atoms[i].X - agB->atoms[j].X;
			
			forcey = forcey + (charge_A*charge_B * distance_y) / ( pow(atom_distance, 3)); 
			forcez = forcez + (charge_A*charge_B * distance_z) / ( pow(atom_distance, 3)); 
			forcex = forcex + (charge_A*charge_B * distance_x) / ( pow(atom_distance, 3)); 
			
		}
		force_in_y[i] = forcey;
		force_in_z[i] = forcez;
		force_in_x[i] = forcex;	
		
	}
}

void calculate_energy_atom(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms, double* Energy_atom )
{

	for(int i = 0; i< agA->natoms; i++ )
	{	
		double energy = 0;
		int atom_type_A = agA->atoms[i].atom_typen;
		double chrg_A = prms->chrgs[atom_type_A];
		for(int j = 0; j < agB->natoms; j++) 
		{
			int atom_type_B = agB->atoms[j].atom_typen;
			double chrg_B = prms->chrgs[atom_type_B];
			float distance = calculate_dist_between_atoms(agA->atoms[i], agB->atoms[j]);
			energy = energy + (chrg_A * chrg_B)/distance;
		}
		Energy_atom[i] = energy;	
	}
}

double calculate_energy(struct atomgrp* agA, struct atomgrp* agB, struct prms* prms)
{
	double energy = 0;
	for(int i = 0; i< agA->natoms; i++ )
	{	
		int atom_type_A = agA->atoms[i].atom_typen;
		double chrg_A = prms->chrgs[atom_type_A];
		for(int j = 0; j < agB->natoms; j++) 
		{
			int atom_type_B = agB->atoms[j].atom_typen;
			double chrg_B = prms->chrgs[atom_type_B];
			float distance = calculate_dist_between_atoms(agA->atoms[i], agB->atoms[j]);
			energy = energy + (chrg_A * chrg_B)/distance;
		}
	}
	return energy;
	
	
}


