/**
   @Program : l-chem-dist.c
   @Note    : This will calculate shortest path(chemical distance) between two particles that were neighbour in previous time step
   @Compilation : gcc D2.c -L/Data/.gsl-2.2.1/lib/ -lgsl -lgslcblas -lm

   Here various things are calculated from a dump trajectory of gel.
   
   
   execution: ./a.out N dump-series_file-name finalframe 

**/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define l 2
#define mtot (2*l+1)


#define maxNN 64

#define e0str "e10"

#define str "e0-10_S0-1_B-0.01"




void load_frame(int N,FILE **fp, int iframe, double *L,double *tild, double ***xyz, int **type){

  int i,id, itype,tstep; double t,Li,Lf,x_cor,y_cor,z_cor,vx,vy,vz;

  //FILE **fp=fopen(FILEname,"r");;
  
  rewind(*fp);
  int frame_read=0;
  //Skip frist (iframe-1) frame data
  for(i=0;i<iframe*(N+9);i++)fscanf(*fp,"%*[^\n]\n");
  while (frame_read<1){
    fscanf(*fp,"%*[^\n]\n");
    fscanf(*fp,"%d\n",&tstep);
    //printf("time step=%d\n",tstep);
    for(i=0;i<3;i++)fscanf(*fp,"%*[^\n]\n");
    fscanf(*fp,"%*f %*f %lf\n",&t);*tild=t;
    fscanf(*fp,"%lf %lf %*f\n",&Li,&Lf);(*L)=Lf-Li;

    for(i=0;i<2;i++)fscanf(*fp,"%*[^\n]\n");
    for(i=0;i<N;i++){
      fscanf(*fp,"%d %d %lf %lf %lf %lf %lf %lf\n",&id,&itype,&x_cor,&y_cor,&z_cor,&vx,&vy,&vz);
      id-=1;
      (*xyz)[id][0]=x_cor;
      (*xyz)[id][1]=y_cor;
      (*xyz)[id][2]=z_cor;
      (*type)[id]=itype;
    }
    frame_read++;
  }
  
  //printf("loaded\n");


}

int delta(int i, int j){
  int value;
  if(i==j) value=1;
  else value=0;
  return value;
}



int main(int argc, char *argv[])     
{
  FILE *fp,*fp1;
  static char FILE[1024],FILE1[1024],ifile[1024];
  static double L,iL;
  static int N;  
  static long i,j,k,m;
  static double gmax,t;
  
  static double charge,x_cor,y_cor,z_cor,vx,vy,vz,tild,tild0,tild1,tildp;
  static double xr,yr,zr;
  static int initframe,finalframe,del_frame,id,itype,nx,ny,nz;
  static char char_dummy[100];
  static float dummy_float;
  
  static double Li,Lf,dr,delx;

      
 
  static int ipart,ineighbour,ilist;
    
  static float Rc=1.568;
  
  static int mincycle,maxcycle,maxset,minset;
  static int icycle, tstep;


  static int type1; // type1=1 for A and type2=2 for B particle.
  //====================================================================
    
  
  if(argc==4){
    N=atof(argv[1]);
    //argv[2]: filename to read
    finalframe=atof(argv[3]);
    
  }
  else {
    printf("Give N, filename and iframe as arguments\n");
    return 0;
  }
  
    
  char command[128];
  sprintf(command, "mkdir -p %s",str);
  system(command);
  


  
    //--------------------------------------------------------------------
    //=================defining various array ============================
    double **xyz_1,**xyz_P;
      
    xyz_1=malloc(N*sizeof(double *));
    xyz_P=malloc(N*sizeof(double *));
    for(j=0;j<N;j++){

      xyz_1[j]=malloc(3*sizeof(double));
      xyz_P[j]=malloc(3*sizeof(double));
    }
    
    for(j=0;j<N;j++)
      for(k=0;k<3;k++){

	xyz_1[j][k]=0.0;
	xyz_P[j][k]=0.0;
      }
    
    
    int *type;
    type=malloc(N*sizeof(int));

    
    int **neighbourLIST1;
    neighbourLIST1=malloc(N*sizeof(int *));
    for(j=0;j<N;j++)neighbourLIST1[j]=malloc(maxNN*sizeof(int));
    int *nlist1;
    nlist1=malloc(N*sizeof(int));


    int **neighbourLISTP;
    neighbourLISTP=malloc(N*sizeof(int *));
    for(j=0;j<N;j++)neighbourLISTP[j]=malloc(maxNN*sizeof(int));
    int *nlistP;
    nlistP=malloc(N*sizeof(int));
    


   
    //-------------------------------------------------------------------
    


    int *LIST,*ID;
    LIST =malloc(N*sizeof(int));
    ID   =malloc(N*sizeof(int));
        
    int iroot,inew,n,ncl,k1,k2,ntot;
    int ifinal,dist,iflag;
    
    int *BondBreak;
    BondBreak=malloc(N*sizeof(int));
    //===================================================================

    
    sprintf(FILE,"%s",argv[2]);
    fp=fopen(FILE,"r");
    printf("file open %s\n",FILE);
    
    
    load_frame(N,&fp,finalframe,&L,&tild1,&xyz_1,&type);

    load_frame(N,&fp,finalframe-1,&L,&tildp,&xyz_P,&type);
    
    fclose(fp);
    
    iL=1.0/L;
    
    exit;
      //==========Transform co-ordinate to wrapped one========== 
      for(i=0;i<N;i++){
	for(j=0;j<3;j++){
	
	  
	  while(xyz_1[i][j]<0)xyz_1[i][j]+=L;
	  while(xyz_1[i][j]>L)xyz_1[i][j]-=L;
	  
	  while(xyz_P[i][j]<0)xyz_P[i][j]+=L;
	  while(xyz_P[i][j]>L)xyz_P[i][j]-=L;
	  
	}
      }
      //========================================================= 




        //========== Get the neighbour list of frame 1 ======================/
      // We getting this only for shortest path caluclation, no dr required.
      for(i=0;i<N;i++)nlist1[i]=0;
      
      for(ipart=0;ipart<N-1;ipart++){
	
	for(ineighbour=ipart+1;ineighbour<N;ineighbour++){

  	  xr=xyz_1[ipart][0]-xyz_1[ineighbour][0];
  	  yr=xyz_1[ipart][1]-xyz_1[ineighbour][1];
  	  zr=xyz_1[ipart][2]-xyz_1[ineighbour][2];

	  xr=xr-tild1*round(yr*iL);
  	  xr=xr-L*round(xr*iL);
  	  yr=yr-L*round(yr*iL);
  	  zr=zr-L*round(zr*iL);
  	  dr=sqrt(xr*xr+yr*yr+zr*zr);
	  
	  //if(dr<Rc[type[ipart]][type[ineighbour]]){
	  Rc=(0.88+(type[ipart]-1)*0.02+(type[ineighbour]-1)*0.02)*1.4;
	  if(dr<Rc){
	    

	    
	    neighbourLIST1[ipart][nlist1[ipart]++]=ineighbour;
	    neighbourLIST1[ineighbour][nlist1[ineighbour]++]=ipart;
	    if(nlist1[ineighbour]>maxNN || nlist1[ipart]> maxNN){
		printf("nlist exeeds maxNN, increase maxNN\n");
		return(0);
	      }
	  }
	}
      }
      //========================================================= 

     
        //========== Get the neighbour list of previous frame marked as xyz_P
      // We getting this only for shortest path caluclation.
      for(i=0;i<N;i++)nlistP[i]=0;
      
      for(ipart=0;ipart<N-1;ipart++){
	
	for(ineighbour=ipart+1;ineighbour<N;ineighbour++){

  	  xr=xyz_P[ipart][0]-xyz_P[ineighbour][0];
	  yr=xyz_P[ipart][1]-xyz_P[ineighbour][1];
	  zr=xyz_P[ipart][2]-xyz_P[ineighbour][2];

	  xr=xr-tildp*round(yr*iL);
  	  xr=xr-L*round(xr*iL);
  	  yr=yr-L*round(yr*iL);
  	  zr=zr-L*round(zr*iL);
  	  dr=sqrt(xr*xr+yr*yr+zr*zr);
	  
	  
	  Rc=(0.88+(type[ipart]-1)*0.02+(type[ineighbour]-1)*0.02)*1.4;
	  if(dr<Rc){
	    neighbourLISTP[ipart][nlistP[ipart]++]=ineighbour;
	    neighbourLISTP[ineighbour][nlistP[ineighbour]++]=ipart;
	    if(nlistP[ineighbour]>maxNN || nlistP[ipart]> maxNN){
		printf("nlist exeeds maxNN, increase maxNN\n");
		return(0);
	      }
	  }
	}
      }
      //========================================================= 
 


      
      
      //=========================================================
      //                Shortest path calculation
      //=========================================================


      
      
      // Let's count cluster for the frame!!!
      
      for(ipart=0;ipart<N;ipart++)BondBreak[ipart]=0;
      
      for(ipart=0;ipart<N;ipart++){
	for(ineighbour=0;ineighbour<nlistP[ipart];ineighbour++){
	  iroot= ipart;
	  ifinal=neighbourLISTP[ipart][ineighbour];
	  
	  
	  
	  ntot=0;ncl=0;
	  for(i=0;i<N;i++)ID[i]=-1;
	  
	  if(ID[iroot]==-1){
	    ncl++;
	    ID[iroot]=ncl;
	    n=1;
	    LIST[n]=iroot;
	    k1=n;
	    k2=k1;
	    dist=0;
	    iflag=1;
	    while(k1<=k2 && iflag==1){
	      dist++;
	      for(k=k1;k<=k2;k++){
		iroot=LIST[k];
		
		for(j=0;j<nlist1[iroot];j++){
		  inew=neighbourLIST1[iroot][j];
		  
		  if(ID[inew]==-1){

		    ID[inew]=ncl;
		    n++;
		    LIST[n]=inew;
		    
		    if(inew==ifinal)iflag=0;   
		  }
		}//loop j
	      }//loop k
	      k1=k2+1;
	      k2=n;
	    }//while
	    
	  }//if
	  
	  if(dist>BondBreak[ipart])BondBreak[ipart]=dist;
	  if(dist>BondBreak[ifinal])BondBreak[ifinal]=dist;
	  //if(dist==32)printf("ipart=%d ifinal=%d\n",ipart,ifinal);
	  
	}//ineighbour
      }//ipart
      
      
      


      

      //==========================================================
      
      
    
      int dump=0;
      if(dump==1){
	sprintf(FILE,"%s/configuration_l_frame%d.dump",str,finalframe);
	fp=fopen(FILE,"w");
	fprintf(fp,"%d\n",N);
	fprintf(fp,"Lattice=\" %e 0.0 0.0 %e %e 0.0 0.0 0.0 %e \" Properties=pos:R:3:D2:R:1:q2:R:1:Nb:I:1:l:I:1\n",L,0.0,L,L);
	for(i=0;i<N;i++){
	  //fprintf(fp,"%f %f %f %f %f %d %d\n",xyz_1[i][0],xyz_1[i][1],xyz_1[i][2],D2Narray[i],qlocNarray[i],nlist1[i],BondBreak[i]);
	}
	fclose(fp);
      }
      else{
	sprintf(FILE,"%s/l_frame%d.dump",str,finalframe);
	fp=fopen(FILE,"w");
	for(i=0;i<N;i++){
	  fprintf(fp,"%d\n",BondBreak[i]);
	}
	fclose(fp);



      }

      
      //*/
      
      

      //}//icycle

//}//iset

    
    
    

      //free(xyz_0);
 free(xyz_1);
 free(xyz_P);
    
    return 0;
} //end main








//     0.88 1.232
//     0.90 1.260 
//     0.92 1.288
//     0.94 1.316
//     0.96 1.344
//     0.98 1.372
//     1.00 1.400  
//     0.92 1.288
//     0.94 1.316
//     0.96 1.344
//     0.98 1.372
//     1.00 1.400  
//     1.02 1.428
//     0.96 1.344
//     0.98 1.372
//     1.00 1.400  
//     1.02 1.428
//     1.04 1.456
//     1.00 1.400    
//     1.02 1.428
//     1.04 1.456
//     1.06 1.484
//     1.04 1.456
//     1.06 1.484
//     1.08 1.512
//     1.08 1.512
//     1.10 1.540
//     1.12 1.568


