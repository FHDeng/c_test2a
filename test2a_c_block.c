
#include <stdio.h>
#include <math.h>

/*This is the rewriting of Rocco's matlab code to C*/
#define xnum 151
#define ynum 61
#define nb 5 //number of block

double Max (double *a)
{
    double maxinum=*a;
    int m,n,r;
    for(m=0;m<ynum;m++)
        for(n=0;n<xnum;n++)
    {
        r=m*xnum;
        if(*(a+r+n)>maxinum)
        {
            maxinum=*(a+r+n);
        }
    }
    return maxinum;
}

double Max1 (double a, double b)
{
    if(a>b)
         return a;
    else
         return b;
}

double Min (double *a)
{
    double minimum=*a;
    int m,n,r;
    for(m=0;m<ynum;m++)
        for(n=0;n<xnum;n++)
    {
        r=m*xnum;
        if(*(a+r+n)<minimum)
        {
            minimum=*(a+r+n);
        }
    }
    return minimum;
}

/*double Min (double a[ynum][xnum])
{
    double minimum=a[0][0];
    int m,n;
    for(m=0;m<ynum;m++)
        for(n=0;n<xnum;n++)
    {
        //r=m*xnum;
        if(a[m][n]<minimum)
        {
            minimum=a[m][n];
        }
    }
    return minimum;
}*/

double Min1 (double a, double b)
{
    if(a<b)
         return a;
    else
         return b;
}

void solver(double dx_, double dy_, double dt_, double *h0_, double *h1_,
           double *Kx_, double *Ky_, double *Ss_, double *Q_, double xlS_,
           double xrS_,double v0_,double *x_)
{
    int i,j,r;
    double ka,kb,kc,kd,vxa,vxb,vyc,vyd;
    double kappa_1;

    // internal values following pages 139-140 in Gerya's book
    for(i=1;i<ynum-1;i++)
    {
        for(j=1;j<xnum-1;j++)
        {
            r=i*xnum;
            ka=(*(Kx_+r+j-1)+*(Kx_+r+j))/2;
            kb=(*(Kx_+r+j+1)+*(Kx_+r+j))/2;
            kc=(*(Ky_+r-xnum+j)+*(Ky_+r+j))/2;
            kd=(*(Ky_+r+xnum+j)+*(Ky_+r+j))/2;
            vxa=-ka*(*(h0_+r+j)-*(h0_+r+j-1))/dx_;
            vxb=-kb*(*(h0_+r+j+1)-*(h0_+r+j))/dx_;
            vyc=-kc*(*(h0_+r+j)-*(h0_+r-xnum+j))/dy_;
            vyd=-kd*(*(h0_+r+xnum+j)-*(h0_+r+j))/dy_;
            kappa_1=dt_/(*(Ss_+r+j));
            *(h1_+r+j)=*(h0_+r+j)+*(Q_+r+j)/(*(Ss_+r+j))
            -kappa_1*((vxb-vxa)/dx_+(vyd-vyc)/dy_);
        }

    }

    //fix surface hydraulic values to 0
    for (j=0;j<xnum;j++)
        *(h1_+j)=0;

    //left and right sides boundary condition: no flow
    for (i=1;i<ynum-1;i++)
    {
        r=i*xnum;
        *(h1_+r)=*(h1_+r+1); //left side
        *(h1_+r+xnum-1)=*(h1_+r+xnum-2);//right side
    }

    //Bottom boundary condition: no flow
    for (i=0;i<xnum;i++)
        *(h1_+(ynum-1)*xnum+i)=*(h1_+(ynum-2)*xnum+i);

    //fix vertical flow in the source region at the bottom
    r=xnum*(ynum-1);
    for (j=0;j<xnum;j++)
    {
        if((*(x_+j)>=xlS_)&&(*(x_+j)<=xrS_))
        {
            kc=(*(Ky_+r-xnum+j)+*(Ky_+r+j))/2;
            *(h1_+r+j)=-v0_*dy_/kc+*(h1_+r-xnum+j);
        }
    }
}

void flow(double dx_,double dy_,double *h, double *Kx_, double *Ky_,
          double *v1x_,double *v1y_, double *v1_)
{
    int i,j,r;
    double ka,kb,kc,kd,vxa,vxb,vyc,vyd;
    //internal values following pages 139-140 in Gerya's book
    for(i=1;i<ynum-1;i++)
    {
        for(j=1;j<xnum-1;j++)
        {
            r=i*xnum;
            ka=(*(Kx_+r+j-1)+*(Kx_+r+j))/2;
            kb=(*(Kx_+r+j+1)+*(Kx_+r+j))/2;
            kc=(*(Ky_+r-xnum+j)+*(Ky_+r+j))/2;
            kd=(*(Ky_+r+xnum+j)+*(Ky_+r+j))/2;
            vxa=-ka*(*(h+r+j)-*(h+r+j-1))/dx_;
            vxb=-kb*(*(h+r+j+1)-*(h+r+j))/dx_;
            vyc=-kc*(*(h+r+j)-*(h+r-xnum+j))/dy_;
            vyd=-kd*(*(h+r+xnum+j)-*(h+r+j))/dy_;
            *(v1x_+r+j)=(vxb+vxa)/2;
            *(v1y_+r+j)=(vyd+vyc)/2;
            //*(v1_+r+j)=sqrt((*(v1x_+r+j))*(*(v1x_+r+j))+(*(v1y_+r+j))*(*(v1y_+r+j)));
        }

    }
    //surface values flow
    for(j=1;j<xnum-1;j++)
    {
        ka=(*(Kx_+j-1)+*(Kx_+j))*0.5;
        kb=(*(Kx_+j+1)+*(Kx_+j))*0.5;
        kd=(*(Ky_+j)+*(Ky_+xnum+j))*0.5;
        vxa=-ka*(0-0)/dx_;
        vxb=-kb*(0-0)/dx_;
        vyd=-kd*(*(h+xnum+j)-0)/dy_;
        *(v1x_+j)=0.5*(vxb+vxa);
        *(v1y_+j)=vyd;
        //*(v1_+j)=sqrt((*(v1x_+j))*(*(v1x_+j))+(*(v1y_+j))*(*(v1y_+j)));
    }
    //surface//the left up corner// I think we can calculate the four corners in the end
    kb=(*(Kx_+1)+*(Kx_))*0.5;
    kd=(*(Ky_+xnum)+*(Ky_))*0.5;
    vxb=-kb*(*(h+1)-*(h))/dx_;
    vyd=-kd*(*(h+xnum)-*(h))/dy_;
    *(v1x_)=vxb;
    *(v1y_)=vyd;
    //*(v1_)=sqrt((*v1x_)*(*v1x_)+(*v1y_)*(*v1y_));
    //surface//the right up corner
    ka=(*(Kx_+xnum-2)+*(Kx_+xnum-1))*0.5;
    kd=(*(Ky_+xnum+xnum-1)+*(Ky_+xnum-1))*0.5;
    vxa=-ka*(*(h+xnum-1)-*(h+xnum-2))/dx_;
    vyd=-kd*(*(h+xnum+xnum-1)-*(h+xnum-1))/dy_;
    *(v1x_+xnum-1)=vxa;
    *(v1y_+xnum-1)=vyd;

    //bottom values flow
    r=(ynum-1)*xnum;
    for(j=1;j<xnum-1;j++)
    {

        ka=(*(Kx_+r+j-1)+*(Kx_+r+j))*0.5;
        kb=(*(Kx_+r+j+1)+*(Kx_+r+j))*0.5;
        kc=(*(Ky_+r-xnum+j)+*(Ky_+r+j))*0.5;
        vxa=-ka*(*(h+r+j)-*(h+r+j-1))/dx_;
        vxb=-kb*(*(h+r+j+1)-*(h+r+j))/dx_;
        vyc=-kc*(*(h+r+j)-*(h+r-xnum+j))/dy_;
        *(v1x_+r+j)=(vxa+vxb)*0.5;
        *(v1y_+r+j)=vyc;
    }
    //bottom//left corner
    kb=(*(Kx_+r+1)+*(Kx_+r))*0.5;
    kc=(*(Ky_+r)+*(Ky_+r-xnum))*0.5;
    vxb=-kb*(*(h+r+1)-*(h+r))/dx_;
    vyc=-kc*(*(h+r)-*(h+r-xnum))/dy_;
    *(v1x_+r)=vxb;
    *(v1y_+r)=vyc;
    //bottom//right corner
    ka=(*(Kx_+r+xnum-1)+*(Kx_+r+xnum-2))*0.5;
    kc=(*(Ky_+r+xnum-1)+*(Ky_+r-xnum+xnum-1))*0.5;
    vxa=-ka*(*(h+r+xnum-1)-*(h+r+xnum-2))/dx_;
    vyc=-kc*(*(h+r+xnum-1)-*(h+r-xnum+xnum-1))/dy_;
    *(v1x_+r+xnum-1)=vxa;
    *(v1y_+r+xnum-1)=vyc;

    //left
    for(i=1;i<ynum-1;i++)
    {
        r=i*xnum;
        kb=(*(Kx_+r+1)+*(Kx_+r))*0.5;
        kc=(*(Ky_+r)+*(Ky_+r-xnum))*0.5;
        kd=(*(Ky_+r+xnum)+*(Ky_+r))*0.5;
        vxb=-kb*(*(h+r+1)-*(h+r))/dx_;
        vyc=-kc*(*(h+r)-*(h+r-xnum))/dy_;
        vyd=-kd*(*(h+r+xnum)-*(h+r))/dy_;
        *(v1x_+r)=vxb;
        *(v1y_+r)=(vyc+vyd)*0.5;
    }
    //left//up corner
    kb=(*(Kx_+1)+*(Kx_))*0.5;
    kd=(*(Ky_+xnum)+*(Ky_))*0.5;
    vxb=-kb*(*(h+1)-*(h))/dx_;
    vyd=-kd*(*(h+xnum)-*(h))/dy_;
    *(v1x_)=vxb;
    *(v1y_)=vyd;
    //left//bottom corner
    r=xnum*(ynum-1);
    kb=(*(Kx_+r+1)+*(Kx_+r))*0.5;
    kc=(*(Ky_+r)+*(Ky_+r-xnum))*0.5;
    vxb=-kb*(*(h+r+1)-*(h+r))/dx_;
    vyc=-kc*(*(h+r)-*(h+r-xnum))/dy_;
    *(v1x_+r)=vxb;
    *(v1y_+r)=vyc;

    //right
    for(i=1;i<ynum-1;i++)
    {
        r=i*xnum;
        ka=(*(Kx_+r+xnum-1)+*(Kx_+r+xnum-2))*0.5;
        kc=(*(Ky_+r+xnum-1)+*(Ky_+r-xnum+xnum-1))*0.5;
        kd=(*(Ky_+r+xnum-1)+*(Ky_+r+xnum+xnum-1))*0.5;
        vxa=-ka*(*(h+r+xnum-1)-*(h+r+xnum-2))/dx_;
        vyc=-kc*(*(h+r+xnum-1)-*(h+r-xnum+xnum-1))/dy_;
        vyd=-kd*(*(h+r+xnum+xnum-1)-*(h+r+xnum-1))/dy_;
        *(v1x_+r+xnum-1)=vxa;
        *(v1y_+r+xnum-1)=(vyc+vyd)*0.5;
    }
    //right//top corner
    ka=(*(Kx_+xnum-1)+*(Kx_+xnum-2))*0.5;
    kd=(*(Ky_+xnum+xnum-1)+*(Ky_+xnum-1))*0.5;
    vxa=-ka*(*(h+xnum-1)-*(h+xnum-2))/dx_;
    vyd=-kd*(*(h+xnum+xnum-1)-*(h+xnum-1))/dy_;
    *(v1x_+xnum-1)=vxa;
    *(v1y_+xnum-1)=vyd;
    //right//bottom corner
    r=xnum*(ynum-1);
    ka=(*(Kx_+r+xnum-1)+*(Kx_+r+xnum-2))*0.5;
    kc=(*(Ky_+r+xnum-1)+*(Ky_+r-xnum+xnum-1))*0.5;
    vxa=-ka*(*(h+r+xnum-1)-*(h+r+xnum-2))/dx_;
    vyc=-kc*(*(h+r+xnum-1)-*(h+r-xnum+xnum-1))/dy_;
    *(v1x_+r+xnum-1)=vxa;
    *(v1y_+r+xnum-1)=vyc;

    //calculate v1
    for(i=0;i<ynum;i++)
        for(j=0;j<xnum;j++)
    {
        r=i*xnum;
        *(v1_+r+j)=sqrt(*(v1x_+r+j)*(*(v1x_+r+j))+*(v1y_+r+j)*(*(v1y_+r+j)));
    }

}

void main()
{
    int i=0,j=0,r=0;
    //Model size
    double xsize=5e5;
    double ysize=5e4;
    //Number of time steps
    int maxtnum=200000;
    //grid size
    double dx=xsize/(xnum-1);
    double dy=ysize/(ynum-1);
    //Grid step
    double xstp,ystp;
    double x[xnum],y[ynum];

    xstp=xsize/(xnum-1); //horizontal
    ystp=ysize/(ynum-1); //vertical
    for(i=0;i<xnum;i++)
    {
        *(x+i)=xstp*i;
    }
    for(i=0;i<ynum;i++)
    {
        *(y+i)=ystp*i;
    }

    //Initial hydraulic uniform head value
    double H0=0;
    double h0[ynum][xnum];
    double h1[ynum][xnum];
    double *h0p,*h1p;
    h0p=&h0[0][0];
    h1p=&h1[0][0];

    //background permeability for x and y direction.
    double kxbg=1;
    double kybg=1.5;

    //Source at the bottom
    double Ls=3500;
    double xlS=0.5*xsize-0.5*Ls;
    double xrS=0.5*xsize+0.5*Ls;
    double v0=100/Ls;
    //printf("\n%f,%f",xlS,xrS);

    //high K blocks
    //int nb=5; //number of blocks
    double Kbx[nb],Kby[nb]; //K for different block
    double xc[nb],yc[nb];//positions of the center points of blocks
    double L[nb],W[nb];//length and width of block
    double xlb[nb],xrb[nb],ytb[nb],ybb[nb];//edges of the blocks
    //Block a
    xc[0]=xsize*0.5-3e4;
    yc[0]=25e3;
    L[0]=28e3;W[0]=10e3;
    Kbx[0]=kxbg*5;
    Kby[0]=kybg*5;

    //block b
    xc[1]=xsize*0.5+5e3;
    yc[1]=20e3;
    L[1]=15e3;W[1]=15e3;
    Kbx[1]=kxbg*4;
    Kby[1]=kybg*4;

    //block c
    xc[2]=xsize*0.5+50e3;
    yc[2]=18e3;
    L[2]=12e3;W[2]=12e3;
    Kbx[2]=kxbg*15;
    Kby[2]=kybg*5;

    //block d
    xc[3]=xsize*0.5+50e3+80e3;
    yc[3]=15e3;
    L[3]=12e3;W[3]=12e3;
    Kbx[3]=kxbg*15;
    Kby[3]=kybg*15;

    //block e
    xc[4]=xsize*0.5+110e3;
    yc[4]=14e3;
    L[4]=11e3;W[4]=11e3;
    Kbx[4]=kxbg*25;
    Kby[4]=kybg*25;

    for(i=0;i<nb;i++)
    {
        xlb[i]=xc[i]-L[i]*0.5;
        xrb[i]=xc[i]+L[i]*0.5;
        ytb[i]=yc[i]-W[i]*0.5;
        ybb[i]=yc[i]+W[i]*0.5;
        //printf("%f\t%f\t%f\t%f\t\%f\t\%f\n",xlb[i],xrb[i],ytb[i],ybb[i],Kbx[i],Kby[i]);
    }

    //Specific storage uniform......
    double ss=1;
    double Ss[ynum][xnum],Kx[ynum][xnum],Ky[ynum][xnum],Q[ynum][xnum];
    double v1[ynum][xnum],v1x[ynum][xnum],v1y[ynum][xnum];

    for (i=0;i<ynum;i++)
    {
        for(j=0;j<xnum;j++)
        {
            h0[i][j]=0;
            h1[i][j]=0;
            Ss[i][j]=ss;
            Kx[i][j]=kxbg;
            Ky[i][j]=kybg;
            Q[i][j]=0;
            v1[0][0]=v1x[0][0]=v1y[0][0]=0;
        }
    }

    //set the K for different block
    for(r=0;r<nb;r++)
    {
        for(i=0;i<ynum;i++)
        {
            for(j=0;j<xnum;j++)
            {
                if((y[i]<=ybb[r])&&(y[i]>=ytb[r])&&(x[j]>=xlb[r])&&(x[j]<=xrb[r]))
                {
                    Kx[i][j]=Kbx[r];
                    Ky[i][j]=Kby[r];
                }

            }

        }
         //printf("%d\t",r);

    }

    //diffusivity for computing time step
    double kappax=Max(&Kx[0][0])/Min(&Ss[0][0]);
    double kappay=Max(&Ky[0][0])/Min(&Ss[0][0]);
    double kappa=Max1(kappax,kappay);

    //time step limit, [book]introduction to numerical geodynamic modelling p134, equation (10.4)
    double dtexp=100*100.0/(3*kappa);
    //time step
    double dt=1.0*dtexp;
    double time=0;
    int it;

    for(it=1;it<=maxtnum;it++)
    {
         solver(dx,dy,dt,h0p,h1p,&Kx[0][0],&Ky[0][0],&Ss[0][0],&Q[0][0],xlS,xrS,v0,&x[0]);
         //h0p=h1p;
         flow(dx,dy,&h1[0][0],&Kx[0][0],&Ky[0][0],&v1x[0][0],&v1y[0][0],&v1[0][0]);
         for(i=0;i<ynum;i++)
            for(j=0;j<xnum;j++)
            h0[i][j]=h1[i][j];
         time=time+dt;
    }

    //calculate the min of head and v1 for reference
    /*double min_h,min_v1;
    min_h=Min(&h1[0][0]);
    min_v1=Min(&v1[0][0]);
    printf("%g\n%g\n",min_h,min_v1);*/


    FILE *fp1,*fp2;
    FILE *fp3,*fp4;
    FILE *fp5,*fp6,*fp7,*fp8;
    fp1=fopen("head_h1.txt","w");
    fp2=fopen("v1.txt","w");
    fp3=fopen("v1x.txt","w");
    fp4=fopen("v1y.txt","w");
    fp5=fopen("Kx.txt","w");
    fp6=fopen("Ky.txt","w");
    fp7=fopen("xy.txt","w"); //storage x and y grid step
    fp8=fopen("paramaters.txt","w");

    for(i=0;i<ynum;i++)
    {
        for(j=0;j<xnum;j++)
        {
            //sometimes, if the results are to small (like 1e-100), we should be careful when write the result to file. maybe we can try %g
            fprintf(fp1,"%.100lf\t",h1[i][j]); //%15.10lf
            fprintf(fp2,"%.100lf\t",v1[i][j]);
            fprintf(fp3,"%.100lf\t",v1x[i][j]);
            fprintf(fp4,"%.100lf\t",v1y[i][j]);
            fprintf(fp5,"%f\t",Kx[i][j]);
            fprintf(fp6,"%f\t",Ky[i][j]);
        }
        fprintf(fp1,"\n");
        fprintf(fp2,"\n");
        fprintf(fp3,"\n");
        fprintf(fp4,"\n");
        fprintf(fp5,"\n");
        fprintf(fp6,"\n");

    }

    for(i=0;i<ynum;i++)
        fprintf(fp7,"%f\t",y[i]);
    fprintf(fp7,"\n");

    for(j=0;j<xnum;j++)
        fprintf(fp7,"%f\t",x[j]);

    fprintf(fp8,"%lf\t%lf\t%d\t%d\t%d\t%lf",xsize,ysize,xnum,ynum,it-1,time);

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    fclose(fp6);
    fclose(fp7);
    fclose(fp8);

    //plot the results in matlab


}
