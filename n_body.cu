#include <iostream>
#include <cuda_runtime.h>
#include<cmath>


const double NEWTON_G = 6.67384e-11;
const double SOFTENING =  1e-9f;

__constant__ double NEWTON_GG = 6.67384e-11;
__constant__ double SOFTENINGG =  1e-9f;

void writeSoA(double** f, int B, int size, const char *filename){
	FILE* file;
    file=fopen(filename, "w");
	
	fprintf(file, "%d\n", B);
	for (int i = 0; i < 5; i++){
		for (int j = 0; j < size/5; j++){
			fprintf(file, "%lf ", ((*f)[j + i*(size/5)]));
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void ReadSoA(double** f, int* B, const char *filename) {
    FILE *fp;
    fp = fopen(filename, "r");
    fscanf(fp, "%d", B);
    int size = (*B) * 5; // 5 atributos: masa, pos_x, pos_y, vel_x, vel_y
    
    double* F = new double[size];
    int i;
    for (i = 0; i < size; i++) {
    	fscanf(fp, "%lf ", &(F[i]));
    }
    *f = F;
    fclose(fp);
}


void printArray(int size, double *arr) {
    printf("[");
	for (int i = 0; i < size; i++) {
		printf("%lf ", arr[i]);
	}
    printf("]\n");
}


void N_body_CPU(int size, double delta_t, double *f, double *fout, int T){
    //printf("N_body_CPU\n");
    
    for (int body_i = 0; body_i < size/5; body_i++){ // para cada cuerpo
        //printf("body %d ", body_i);
        //if (body_i == 95 && T == 0) printf("(%lf, %lf, %lf, %lf, %lf)\n", f[body_i], f[body_i + (size/5)], f[body_i + (size/5)*2], f[body_i + (size/5)*3], f[body_i + (size/5)*4]);
        double mass1  = f[body_i];
        double x1 = f[body_i + (size/5)];
        double y1 = f[body_i + (size/5)*2];
        double vx1 = f[body_i + (size/5)*3];
        double vy1 = f[body_i + (size/5)*4];

    
        double Fx = 0;
        double Fy = 0;
        for (int j = 0; j < size/5; j++){ // comparar con otros cuerpos
            if (j == body_i) continue; // creo que puedo obviarlo pues el radio seria cero (aunque nos da division por cero)
            double mass2  = f[j];
            double x2 = f[j + (size/5)];
            double y2 = f[j + (size/5)*2];


            double distance =  sqrt( pow(x2-x1, 2) + pow(y2-y1, 2) + pow(SOFTENING, 2) ); 
            //if (body_i == 0) printf("distance: %lf\n", distance);
            double angle = atan((y2-y1)/(x2-x1));

            Fx +=  NEWTON_G*mass1*mass2/(pow(distance, 2)) * cos(angle); 
            Fy +=  NEWTON_G*mass1*mass2/(pow(distance, 2)) * sin(angle); 
            
        }
        

        double new_vx1 = vx1 + Fx*delta_t/mass1;
        double new_vy1 = vy1 + Fy*delta_t/mass1;

    
        // a futuro, usar otro arreglo para la masa, pues no cambia
        fout[body_i] = mass1; 
        fout[body_i + (size/5)]   = x1 + new_vx1*delta_t; //new x
        fout[body_i + (size/5)*2] = y1 + new_vy1*delta_t ; //new  y
        fout[body_i + (size/5)*3] = new_vx1; //new vx
        fout[body_i + (size/5)*4] = new_vy1; //new vy
    }
}

__global__ void N_body_GPU(int size, double delta_t, double *f, double *fout){
    int body_i= threadIdx.x + blockDim.x*blockIdx.x;
    if (body_i < size/5){

        double mass1  = f[body_i];
        double x1 = f[body_i + (size/5)];
        double y1 = f[body_i + (size/5)*2];
        double vx1 = f[body_i + (size/5)*3];
        double vy1 = f[body_i + (size/5)*4];
        double mass2,x2,y2,distance,angle,new_vx1,new_vy1;
        double Fx = 0;
        double Fy = 0;
        for (int j = 0; j < size/5; j++){ // comparar con otros cuerpos
                if (j != body_i){
                mass2  = f[j];
                x2 = f[j + (size/5)];
                y2 = f[j + (size/5)*2];

                distance =  sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + SOFTENINGG); 
                //printf("distance: %lf\n", distance);
                angle = atan((y2-y1)/(x2-x1));

                Fx +=  NEWTON_GG*mass2/(pow(distance, 2)) * cos(angle); 
                Fy +=  NEWTON_GG*mass2/(pow(distance, 2)) * sin(angle); 
            }
        }
        //printf("F: %lf\n", F);

        new_vx1 = vx1 + Fx*delta_t;
        new_vy1 = vy1 + Fy*delta_t;

        //printf("F*delta_t/mass1: %lf\n", F*delta_t/mass1);
        //printf("new_vx1: %lf\n", new_vx1);
        //printf("new_vy1: %lf\n", new_vy1);
    
        // a futuro, usar otro arreglo para la masa, pues no cambia
        fout[body_i] = mass1; 
        fout[body_i + (size/5)]   = x1 + new_vx1*delta_t; //new x
        fout[body_i + (size/5)*2] = y1 + new_vy1*delta_t ; //new  y
        fout[body_i + (size/5)*3] = new_vx1; //new vx
        fout[body_i + (size/5)*4] = new_vy1; //new vy
    }
}


__global__ void N_body_GPU_F(int size, double delta_t, double *f, double *fout,int T){
	int body_i= threadIdx.x + blockDim.x*blockIdx.x;
	if (body_i<size/5){
		extern __shared__ double datos[];
		// 5 atributos: masa
		datos[body_i			]= f[body_i];
		datos[body_i+ (size/5)	]= f[body_i + (size/5)];
       	datos[body_i+ (size/5)*2]= f[body_i + (size/5)*2];
		double autx,auty,rx,ry,vx,vy;
        vx=f[body_i+ (size/5)*3];
        vy=f[body_i+ (size/5)*4];
        double angle;
        double Ax,Ay;
		fout[body_i]=datos[body_i];
        for (int t = 0; t < T; ++t){
            __syncthreads();
            Ax=0.0,Ay=0.0;
            autx=datos[body_i+ (size/5)  ];
            auty=datos[body_i+ (size/5)*2];
            auty=datos[body_i+ (size/5)*2];
    		for (int i = 0; i < size/5; ++i){
                if (i!=body_i){
        			rx=autx-datos[i+ (size/5)  ];
        			ry=auty-datos[i+ (size/5)*2];
        			angle=atan(ry/rx);
                    rx=datos[i]/sqrt(rx*rx+ry*ry+SOFTENINGG);
        			Ax += rx*cos(angle);
        			Ay += rx*sin(angle);
                }
    		}
    		Ax*=NEWTON_GG*delta_t;
    		Ay*=NEWTON_GG*delta_t;
    		datos[body_i+ (size/5)  ]=autx+vx*delta_t+Ax*delta_t;
    		datos[body_i+ (size/5)*2]=auty+vy*delta_t+Ay*delta_t;
    		vx=Ax+vx;
    		vy=Ay+vy;
        }
        fout[body_i+ (size/5)  ]=autx+vx*delta_t+Ax*delta_t;
        fout[body_i+ (size/5)*2]=auty+vy*delta_t+Ay*delta_t;
        fout[body_i+ (size/5)*3]=Ax+vx;
        fout[body_i+ (size/5)*4]=Ay+vy;
	}
}

int main() {	


    cudaEvent_t ct1, ct2, ct3, ct4;

    clock_t t1, t2;
    double ms;

    char filename[] = "input.txt";
    char filename_out[] = "-CPU-Resultado.txt";
    char filename_aux[30];

    char final[] = "final";
    char directory[] = "data/";
    char directory_aux[30];

	float dt,dt2;

	//int iterator=2;
	int B;
	
    double *f, *fout, *fhost, *fhostout, *faux,*ff;

    int grid_size, block_size = 256;

	ReadSoA(&fhost, &B, filename);
    int size = B*5;

    cudaMalloc((void**)&f, size* sizeof(double));
    cudaMalloc((void**)&ff, size* sizeof(double));
    cudaMalloc((void**)&fout, size* sizeof(double));

    cudaMemcpy(f, fhost, size* sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(ff, fhost, size* sizeof(double), cudaMemcpyHostToDevice);
    int debug = 0;
    
    if (debug){
        printf("B: %d\n",  B );
        printf("size: %d\n", size);
        printArray(size, fhost);
    } 

    double delta_t = 0.01;

    fhostout = new double[size];

    if (debug){
        printArray(size, fhost);
    }

    long T = 10000;
    
    char integer_string[32];
    char integer_string2[32];

    /*** CPU ***/
    int cpu = 1;
    if (cpu){ 

        t1 = clock();
        for (long t = 0; t < T; t++){
            N_body_CPU(size, delta_t, fhost, fhostout, t);
            faux=fhost;
            fhost = fhostout;
            fhostout=faux;
            
            if (t % 1000 == 0 || t == T-1){
                sprintf(integer_string, "%d", t);
                sprintf(integer_string2, "-%d", T);
                strcpy(filename_aux, filename_out);
                strcpy(directory_aux, directory);

                writeSoA(&fhostout, B, size, 
                          strcat(directory_aux, strcat(integer_string, strcat(integer_string2, filename_aux) ) ));
            }
            
            //printArray(size, fhostout);
            //std::cout << "-----------------------" << std::endl;
        }
        t2 = clock();

        if (debug){
            printArray(size, fhost);
        }

        ms = 1000.0 * (double)(t2 - t1) / CLOCKS_PER_SEC;
    	std::cout << "Tiempo CPU  : " << ms << "[ms]" << std::endl;

        //writeSoA(&fhostout, B, size, strcat(final, filename_out) );
    }
    int long_simulation = 0;
    
    if (long_simulation){
        T = 20000*20;
    }


    /***  GPU ***/
    int gpu1 = 1;
    if (gpu1){

        char filename_outGPU[] = "-GPU-Resultado.txt";
        grid_size = (int)ceil((float) B / block_size);

        cudaEventCreate(&ct1);
        cudaEventCreate(&ct2);
        cudaEventRecord(ct1);
        for (long t = 0; t < T; t++){
            N_body_GPU<<<grid_size, block_size>>>(size, delta_t, f, fout);
            faux = fout;
            fout = f;
            f = faux;
            
            if (t % 1000 == 0 || t == T-1){
                sprintf(integer_string, "%d", t);
                sprintf(integer_string2, "-%d", T);
                strcpy(filename_aux, filename_outGPU);
                strcpy(directory_aux, directory);
                cudaMemcpy(fhostout, f, size* sizeof(double), cudaMemcpyDeviceToHost);

                writeSoA(&fhostout, B, size, 
                            strcat(directory_aux, strcat(integer_string, strcat(integer_string2, filename_aux))) );
            }
        }
        cudaEventRecord(ct2);
        cudaEventSynchronize(ct2);
        cudaEventElapsedTime(&dt, ct1, ct2);

        std::cout << "Tiempo GPU  : " << dt << "[ms]" << std::endl;

        cudaMemcpy(fhostout, f, size* sizeof(double), cudaMemcpyDeviceToHost);
        //strcpy(filename_out, "GPU-Resultado.txt");
        //writeSoA(&fhostout, B, size, filename_outGPU);
    }

    /***  GPU Fast ***/
    char filename_outFPU[] = "data/FPU-Resultado.txt";

    grid_size = (int)ceil((float) B / block_size);

    cudaEventCreate(&ct3);
    cudaEventCreate(&ct4);
    cudaEventRecord(ct3);

    N_body_GPU_F<<<grid_size, block_size,B*3* sizeof(double)>>>(size, delta_t, ff, fout,T);

    cudaEventRecord(ct4);
    cudaEventSynchronize(ct4);
    cudaEventElapsedTime(&dt2, ct3, ct4);

    std::cout << "Tiempo GPU-F: " << dt2 << "[ms]" << std::endl;

    cudaMemcpy(fhostout, fout, size* sizeof(double), cudaMemcpyDeviceToHost);

    writeSoA(&fhostout, B, size, filename_outFPU);
    

    delete[] fhostout;

}
