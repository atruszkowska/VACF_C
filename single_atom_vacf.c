#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "math.h"

#define MAXNAME 250
#define MAXL 250
#define MAXW 650
#define LONGS 5000

// Functions
int parse_string(char* s, char* a[]);
void char_to_float(char* in[], double out[], int nw);

/* Program for computing components of the velocity autocorrelation 
 * function from a LAMMPS dump file for a single atom
 * Syntax: single_atom_vacf <dumpfile> <file_x> <file_y> <file_z> #atom
 * Notes: - <dumpfile> is the standard LAMMPS .d file; 
 *          <file_x>, <file_y>, and <file_z> are output file names for 
 *          x, y, and z components of the single atom VACF;
 *          #atom is the ID number of the atom in question 
 *        - Number of atoms assumed fixed.
          - Structure of the file - order and names of 
            variables in the dumpfile - is hardcoded.
            To change, need to modify the code.
          - Outputs 3 files - for x, y, and z diretion with VACF 
            inner sum components for the atom i.e. each of the files
            has a form:
            v(t0)v(t0+t0)|v(t1)v(t1+t0)|...|v(tN)v(tN+t0)
                .
                .
            v(t0)v(t0+tM)|v(t1)v(t1+tM)|...|v(tN)v(tN+tM)
            where v is atom's velocity, ts are the time steps,
            N is the total number of time steps and M is N-dt
            i.e. - one time step before N
          - Array sizes above (#define ...) may have to be 
            made larger based on the input file      
   Last modified: December 19 2015
*/     
          
int main(int argc, char *argv[])
{
	// Miscelaneous variables
	char *file_in, *file_out_1, *file_out_2, *file_out_3;
	char *sp;
	char s[MAXL];
	int frame=0, nat, k, itm, ktm, jtm, k2;
	int iat,nw,ind,tmp;
	char *wa[MAXW];
	double w[MAXW];


	// Pointers to files
	FILE *fpi, *fpo1, *fpo2, *fpo3;
	
	// Usage check
	if(argc!=6){
		printf("Error\n");
		printf("Syntax single_atom_vacf <dumpfile> <file_x> <file_y> <file_z> #atom\n");
		exit(1);
	}
	// Get the names and open the dump file
	file_in=argv[1];
	file_out_1=argv[2];
	file_out_2=argv[3];
	file_out_3=argv[4];
	ktm=atoi(argv[5]);
	// 
	fpi=fopen(file_in,"r");

	// Get the number of frames and the number
	// of atoms
	while((sp = fgets(s,MAXL,fpi))!= NULL){
		// Collect the timestep
		if(!strcmp(sp,"ITEM: TIMESTEP\n"))
		    frame++;
		    if(!strcmp(sp,"ITEM: NUMBER OF ATOMS\n"))
		        nat = atoi(fgets(s,MAXL,fpi));
	}
	// Print the information
	printf("Number of frames: %d and number of atoms %d\n",frame, nat);
	printf("Chosen atom: %d \n", ktm);
	rewind(fpi);

	// Allocate space for arrays	 
	// ----------------------------------------------------------
	// VACF_x components
	double** VACF_x= (double**) malloc(frame*sizeof(double*));
	for(k=0;k<frame;k++)
		VACF_x[k]= (double*) malloc(frame*sizeof(double));
	// Initialize
	for(k=0;k<frame;k++)
		for(k2=0;k2<frame;k2++)
		    VACF_x[k][k2]=0.0;
	// VACF_y components
	double** VACF_y= (double**) malloc(frame*sizeof(double*));
	for(k=0;k<frame;k++)
		VACF_y[k]= (double*) malloc(frame*sizeof(double));
	// Initialize
	for(k=0;k<frame;k++)
		for(k2=0;k2<frame;k2++)
		    VACF_y[k][k2]=0.0;
	// VACF_z components
	double** VACF_z= (double**) malloc(frame*sizeof(double*));
	for(k=0;k<frame;k++)
		VACF_z[k]= (double*) malloc(frame*sizeof(double));
	// Initialize
	for(k=0;k<frame;k++)
		for(k2=0;k2<frame;k2++)
		    VACF_z[k][k2]=0.0;
	// Temporary velocity storage
	double** vel_temp= (double**) malloc(frame*sizeof(double*));
	for(k=0;k<frame;k++)
		vel_temp[k]= (double*) malloc(3*sizeof(double));

	// Main computation
	// ----------------------------------------------------------
	// VACF components for a single atom
    // Loop through the whol file, store velocity components 
    // from each frame, than proceed with computations
    // for that one atom
	tmp=0;
	while((sp = fgets(s,MAXL,fpi))!= NULL){
		// Change the header here if needed
		if(!strcmp(sp,"ITEM: ATOMS id x y z vx vy vz c_myPE c_myKE c_myStress[1] c_myStress[2] c_myStress[3] c_myStress[4] c_myStress[5] c_myStress[6] \n")){
		    iat=0;
		    // Get the data as a string, parse it into words and
		    // then convert each word to a float 
		    while(iat<nat){
		        sp = fgets(s,MAXL,fpi);
		        nw=parse_string(s,wa);
		        char_to_float(wa,w,nw);
		        ind=w[0]-1;
		        // If it is the current atom ktm 
		        // save the velocity data, break 
		        // the inner while loop and move on 
		        // to the next frame
		        if(ind+1==ktm){
		            vel_temp[tmp][0]=w[4];
		            vel_temp[tmp][1]=w[5];                            
		            vel_temp[tmp][2]=w[6];
		            tmp++;
		            break;
		        }
		        iat++;
		    }
		}
	}
	for(itm=0;itm<frame;itm++){
	   for(jtm=0;jtm<=frame-itm-1;jtm++){
		    VACF_x[itm][jtm]=vel_temp[jtm][0]*vel_temp[jtm+itm][0];
		    VACF_y[itm][jtm]=vel_temp[jtm][1]*vel_temp[jtm+itm][1];
		    VACF_z[itm][jtm]=vel_temp[jtm][2]*vel_temp[jtm+itm][2];
		}
	}

	// Write output
	// ----------------------------------------------------------
	fpo1=fopen(file_out_1,"w");
	for(jtm=0;jtm<frame;jtm++){
		for(itm=0;itm<frame;itm++) 
		    fprintf(fpo1, "%f ", VACF_x[jtm][itm]);
		fprintf(fpo1,"\n");
	}
	fclose(fpo1);

	fpo2=fopen(file_out_2,"w");
	for(jtm=0;jtm<frame;jtm++){
		for(itm=0;itm<frame;itm++) 
		    fprintf(fpo2, "%f ", VACF_y[jtm][itm]);
		fprintf(fpo2,"\n");
	}
	fclose(fpo2);

	fpo3=fopen(file_out_3,"w");
	for(jtm=0;jtm<frame;jtm++){
		for(itm=0;itm<frame;itm++) 
		    fprintf(fpo3, "%f ", VACF_z[jtm][itm]);
		fprintf(fpo3,"\n");
	}
	fclose(fpo3);

	// Free arrays
	// ----------------------------------------------------------
	for(k=0;k<frame;k++)
		free(VACF_x[k]);
	free(VACF_x);
	for(k=0;k<frame;k++)
		free(VACF_y[k]);
	free(VACF_y);
	for(k=0;k<frame;k++)
		free(VACF_z[k]);
	free(VACF_z);
	for(k=0;k<frame;k++)
		free(vel_temp[k]);
	free(vel_temp);

	// Close the input file
	// ----------------------------------------------------------
	fclose(fpi);
}

/* Function to parse a line of input into an aray of words */
/* s - string to be parsed
 * a - string with parsed elements */
int parse_string(char* s,char* a[])
{
    int nw,j; 
    a[0] = strtok(s," \t\n\r\v\f"); 
    nw = 1;				 
    while((a[nw]= strtok(NULL," \t\n\r\v\f"))!=NULL)
        nw++;
   return nw;
}

/* Function to convert array of words to array
 * of doubles
 * in[] - string with pointers to words
 * out[] - string with doubles */ 
void char_to_float(char* in[], double out[], int nw)
{
    int k;
    for(k=0;k<nw;k++)
        out[k]=atof(in[k]);
}


