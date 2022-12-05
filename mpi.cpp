using namespace std;
#include <stdio.h>
#include<vector>
#include <mpi.h>
#include <iostream>
#include <array>
#include<cmath>
#include<math.h> 
#include<omp.h>

int N = pow(2,24);
int stages = log2(N);
unsigned int reverseBits(unsigned int n)
{
	unsigned int rev = 0;

	for (int i=0; i<stages; i++){
        rev <<= 1;
        if ((n & 1)==1){
            rev++;
        }
        n>>=1;

    }

	// required number
	return rev;
}

int main(int argc, char** argv){
    int process_Rank, size_Of_Cluster, message_Item, offset, chunksize, leftover;
    int k,j,i;
    double initial_time_1, initial_time_2, initial_time_3, initial_time_4;
    //int N = pow(2,25);
    double pi = 3.1415926535897932384626;
    vector<vector<double>> sample (N, vector<double>(2,0));
    for (int i=0; i<N; i++){
        sample[i][0]=(sin((2*pi*i/N)));
    }
    //vector<vector<double>> sample = {{0,0},{1,0},{0,0},{-1,0},{0,0},{1,0},{0,0},{-1,0}};
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    //int N = sample.size();
    //int wr = size_Of_Cluster-1;//number of workers
    chunksize = N/size_Of_Cluster; 
    int stages = log2(N);
    double twiddle_real[N/2];
    double twiddle_img[N/2];
    double array_real[N];
    double array_img[N];
    initial_time_1 = MPI_Wtime();
    if (process_Rank == 0){
        //double pi = 3.1415926535897932384626;
        
        double del_theta = 2*pi/N;
        double c_del = cos(del_theta);
        double s_del = sin(del_theta);
        double sample_real[N];
        double sample_complex[N];
        twiddle_real[0] = 1;
        twiddle_img[0] = 0;
        for (int m =1; m<N/2; m++){ //twiddle factors
            twiddle_real[m] = c_del*twiddle_real[m-1] - s_del*twiddle_img[m-1];
            twiddle_img[m] = c_del*twiddle_img[m-1] + s_del*twiddle_real[m-1];
            //cout<< tf[m][0]<<endl;
        }
        for (int i = 0; i<N; i++){
        sample_real[i]= sample[i][0];
        sample_complex[i]= sample[i][1];
        
        }
        for (int i=0; i<N; i++){
            int index = reverseBits(i);
            sample[i][0]= sample_real[index];
            sample[i][1]= sample_complex[index];
        }
        //diveide the sample into real and imaginary array
        for (int i=0;i<N;i++){
            array_real[i]= sample[i][0];
            array_img[i]= sample[i][1];
        }
        ;  
    }
    initial_time_2 = MPI_Wtime();
    MPI_Bcast(&twiddle_real, N/2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&twiddle_img, N/2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    initial_time_3 = MPI_Wtime();
    double para_time;
    for (k=0; k<log2(N/size_Of_Cluster); k++){//only go through the stages where the number of blocks is more than 
        int num_blks_stage_proc = chunksize/pow(2,k+1);
        int num_bf_block = pow(2,k);
        int bf_span = pow(2,k);
        int blk_step = pow(2, k+1);
        int bf_step=  1;
        int twiddle_index_step = N/pow(2,k+1);
        //create a buffer to hold the chunk of real and imginary data
        double sub_array_real[chunksize];
        double sub_array_img[chunksize];
        //scatter the data to all process
        MPI_Scatter(&array_real, chunksize, MPI_DOUBLE, &sub_array_real,chunksize,MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(&array_img, chunksize, MPI_DOUBLE, &sub_array_img,chunksize,MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //perform the serial fft for each block in each processor. 
        for (int m = 0; m<num_blks_stage_proc; m++){ 
            for (int n = 0; n<num_bf_block; n++){
                //cout<< "value n = "<< n<< endl;
                int blk_pntr = m*blk_step;
                int twiddle_index = n*twiddle_index_step;
                double tf_real = twiddle_real[twiddle_index]; // calling the function everytime is not good 
                double tf_img = -twiddle_img[twiddle_index]; // add the twiddle factor function 
                int input_1 = blk_pntr + n*bf_step;
                int input_2 = input_1 + bf_span;
                double bf_input_1_real =  sub_array_real[input_1];
                double bf_input_1_img =  sub_array_img[input_1];
                double bf_input_2_real = sub_array_real[input_2]*tf_real+ sub_array_img[input_2]*(-tf_img);
                double bf_input_2_img = sub_array_img[input_2]*tf_real - sub_array_real[input_2]*(-tf_img);
                double bf_out_1_real = bf_input_1_real + bf_input_2_real;
                double bf_out_1_img = bf_input_1_img + bf_input_2_img;
                double bf_out_2_real = bf_input_1_real - bf_input_2_real;
                double bf_out_2_img = bf_input_1_img - bf_input_2_img;
                sub_array_real[input_1] =  bf_out_1_real;
                sub_array_img[input_1] =  bf_out_1_img;
                sub_array_real[input_2] = bf_out_2_real;
                sub_array_img[input_2] = bf_out_2_img;
                //cout<<sub_array_real[input_2]<<endl;
            }

        }
        initial_time_4 = MPI_Wtime();
        MPI_Gather(&sub_array_real, chunksize, MPI_DOUBLE, &array_real, chunksize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(&sub_array_img, chunksize, MPI_DOUBLE, &array_img, chunksize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        para_time += MPI_Wtime()-initial_time_4;
    }



    if (process_Rank == 0){
        double initial_time_5 = MPI_Wtime();
        for (int k=log2(N/size_Of_Cluster); k<stages; k++){
            int num_blks_stage_proc= N/(pow(2,k+1));
            int num_bf_block = pow(2,k);
            int bf_span = pow(2,k);
            int blk_step = pow(2, k+1);
            int bf_step=  1;
            int twiddle_index_step = N/pow(2,k+1);
            for (int m = 0; m<num_blks_stage_proc; m++){ 
                for (int n = 0; n<num_bf_block; n++){
                    int blk_pntr = m*blk_step;
                    int twiddle_index = n*twiddle_index_step;
                    double tf_real = twiddle_real[twiddle_index]; // calling the function everytime is not good 
                    double tf_img = -twiddle_img[twiddle_index]; // add the twiddle factor function 
                    int input_1 = blk_pntr + n*bf_step;
                    int input_2 = input_1 + bf_span;
                    double bf_input_1_real =  array_real[input_1];
                    double bf_input_1_img =  array_img[input_1];
                    double bf_input_2_real = array_real[input_2]*tf_real+ array_img[input_2]*(-tf_img);
                    double bf_input_2_img = array_img[input_2]*tf_real - array_real[input_2]*(-tf_img);
                    double bf_out_1_real = bf_input_1_real + bf_input_2_real;
                    double bf_out_1_img = bf_input_1_img + bf_input_2_img;
                    double bf_out_2_real = bf_input_1_real - bf_input_2_real;
                    double bf_out_2_img = bf_input_1_img - bf_input_2_img;
                    array_real[input_1] =  bf_out_1_real;
                    array_img[input_1] =  bf_out_1_img;
                    array_real[input_2] = bf_out_2_real;
                    array_img[input_2] = bf_out_2_img;
                }
            }
        }
        //double initial_time_6 = MPI_Wtime();
        //double computation_time = (initial_time_2-initial_time_1)+(initial_time_4-initial_time_3)+(initial_time_6-initial_time_5);
        //double communication_time = (initial_time_3-initial_time_2)+para_time;
        //double computation_time =  (initial_time_6 - initial_time_1)-communication_time;
        //double total_time = (initial_time_6-initial_time_1);
        //double para_time = final_para_time - initial_time;
        //double seq_time = final_time - final_para_time;
        //cout<<"computation time: "<< computation_time<<endl;
        //cout<<"communication time: "<< communication_time<<endl;
        //cout<<"total time: "<< total_time<<endl;

    }
    double initial_time_6 =  MPI_Wtime();
    double total_time = initial_time_6-initial_time_1;
    if (process_Rank==0){
        cout<<total_time<<endl;
    }
    


    MPI_Finalize();
    return 0;
}