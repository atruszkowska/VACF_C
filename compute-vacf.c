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

/* Program for computing velocity autocorrelation 
 * function from a LAMMPS dump file
 * Syntax: compute-vacf <dumpfile> <outputfile>
 * Notes: - <dumpfile> is the standard LAMMPS .d file
 *          and <outputfile> is the name of the file where
 *          the output is supposed to be written
 *        - Number of atoms assumed fixed.
          - Structure of the file - order and names of 
            variables in the dumpfile - is hardcoded.
            To change, need to modify the code.
          - Outputs a file with Time step | vacf_x | vacf_y | vacf_z
          - Array sizes above (#define ...) may have to be 
            made larger based on the input file
   Last modified: December 17 2015      
*/
            
int main(int argc, char *argv[])
{
    // Miscelaneous variables
    char *file_in, *file_out;
    char *sp;
    char s[MAXL];
    int frame=0, nat, k, itm, ktm, jtm, k2;
    int iat,nw,ind,tmp;
    char *wa[MAXW];
    double w[MAXW];
    
    // Pointers to files
    FILE *fpi, *fpo;
   
   	// Usage check
    if(argc!=3){
        printf("Error\n");
        printf("Syntax compute_vacf <dumpfile> <outputfile>\n");
        exit(1);
    }
    
    // Get the file names and open the dump file
    file_in=argv[1];
    file_out=argv[2];
    // 
    fpi=fopen(file_in,"r");
    
    // Get the number of frames and the number
    // of atoms
    while((sp = fgets(s,MAXL,fpi))!= NULL){
        // Count the frames
        if(!strcmp(sp,"ITEM: TIMESTEP\n"))
            frame++;
        if(!strcmp(sp,"ITEM: NUMBER OF ATOMS\n"))
            nat = atoi(fgets(s,MAXL,fpi));  
    }
    printf("Number of frames: %d and number of atoms %d\n",frame, nat);
    rewind(fpi);
    
    // Time step array
    double* time=(double*) malloc(frame*sizeof(double));
    k=0;
    while((sp = fgets(s,MAXL,fpi))!= NULL){
        // Collect the timestep
        if(!strcmp(sp,"ITEM: TIMESTEP\n")){
            sp = fgets(s,MAXL,fpi);
            nw=parse_string(s,wa);
            char_to_float(wa,w,nw);
            time[k]=w[0];
            k++;
        }
    }
    rewind(fpi);
    
    
    // Allocate and initialize arrays
    // ----------------------------------------------------------
    // Total VACF
    double** VACF= (double**) malloc(frame*sizeof(double*));
    for(k=0;k<frame;k++)
        VACF[k]= (double*) malloc(3*sizeof(double));
    // Initialize VACF
    for(k=0;k<frame;k++)
        for(k2=0;k2<3;k2++)
            VACF[k][k2]=0.0;
    // Temporary velocity storage for every atom
    double** vel_temp= (double**) malloc(frame*sizeof(double*));
    for(k=0;k<frame;k++)
        vel_temp[k]= (double*) malloc(3*sizeof(double));
    // Temporary VACF for an atom
    double* Vtemp=(double*) malloc(3*sizeof(double));
   
	// Main computation
   	// ----------------------------------------------------------
    // Loop for each atom, store velocity components 
    // from each frame, than proceed with computations
    // for that one atom
    for(ktm=1;ktm<=nat;ktm++){
        tmp=0;
        // Track progress
        if(ktm%100==0)
            printf("Atom #=%d \n",ktm);            
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
        // Per atom computation
        for(itm=0;itm<frame;itm++){
            // Reinitialize the sum   
            Vtemp[0]=0.0,Vtemp[1]=0.0,Vtemp[2]=0.0;
            for(jtm=0;jtm<=frame-itm-1;jtm++){
                Vtemp[0]=Vtemp[0]+vel_temp[jtm][0]*vel_temp[jtm+itm][0];
                Vtemp[1]=Vtemp[1]+vel_temp[jtm][1]*vel_temp[jtm+itm][1];
                Vtemp[2]=Vtemp[2]+vel_temp[jtm][2]*vel_temp[jtm+itm][2];
            }
            // Add to the VACF matrix
            VACF[itm][0]=VACF[itm][0]+Vtemp[0]/(frame-itm);
            VACF[itm][1]=VACF[itm][1]+Vtemp[1]/(frame-itm);
            VACF[itm][2]=VACF[itm][2]+Vtemp[2]/(frame-itm);
        }
       rewind(fpi);
    }
    
    // Correct for the number of atoms
    for(k=0;k<frame;k++)
        for(k2=0;k2<3;k2++)
            VACF[k][k2]=VACF[k][k2]/nat;
    
	// Write the result 
    // ----------------------------------------------------------
    // Open the file for writing
    fpo=fopen(file_out,"w");
    // Write VACF
    for(jtm=0;jtm<frame;jtm++){
        fprintf(fpo, "%f %f %f %f\n", time[jtm], VACF[jtm][0],VACF[jtm][1],VACF[jtm][2]);
    }
    
    // Free arrays
    // ----------------------------------------------------------
    for(k=0;k<frame;k++)
        free(VACF[k]);
    free(VACF);
    for(k=0;k<frame;k++)
        free(vel_temp[k]);
    free(vel_temp);
    free(Vtemp);
    free(time);
    
    // Close the files
	// ----------------------------------------------------------
    fclose(fpi);
    fclose(fpo);
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


