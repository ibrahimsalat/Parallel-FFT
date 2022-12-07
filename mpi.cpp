using namespace std;
#include <stdio.h>
#include<vector>
#include <mpi.h>
#include <iostream>
#include <array>
#include<cmath>
#include<math.h> 
#include<omp.h>

int N = pow(2,8);
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
    double initial_time_1, initial_time_2, initial_time_3, initial_time_4, initial_time_5, initial_time_6, initial_time_7, initial_time_8, test_time_1;
    //int N = pow(2,25);
    float pi = 3.1415926535897932384626;
    vector<vector<float>> sample (N, vector<float>(2,0));
    for (int i=0; i<N; i++){
        sample[i][0]=(sin((2*pi*i/N)));
    }
    //vector<vector<float>> sample = {{0,0},{1,0},{0,0},{-1,0},{0,0},{1,0},{0,0},{-1,0}};
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    //int N = sample.size();
    //int wr = size_Of_Cluster-1;//number of workers

    chunksize = N/size_Of_Cluster; 
    int stages = log2(N);
    float twiddle_real[N/2];
    float twiddle_img[N/2];
    float array_real[N];
    float array_img[N];
    float sub_array_real[chunksize];
    float sub_array_img[chunksize];
    int sub_bitrev_index[chunksize];
    //int bitrev_index[N];
    int *bitrev_index;
    bitrev_index = new int[N];
    initial_time_1 = MPI_Wtime();

    for (int i= process_Rank*chunksize; i<(process_Rank*chunksize)+chunksize; i++){
        sub_bitrev_index[i]= reverseBits(i);
    }
    initial_time_7 = MPI_Wtime();

    MPI_Gather(&sub_bitrev_index, chunksize, MPI_FLOAT, &bitrev_index, chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    initial_time_8 = MPI_Wtime();
    if (process_Rank == 0){
        //float pi = 3.1415926535897932384626;
        
        float del_theta = 2*pi/N;
        float c_del = cos(del_theta);
        float s_del = sin(del_theta);
        
        twiddle_real[0] = 1;
        twiddle_img[0] = 0;
        for (int m =1; m<N/2; m++){ //twiddle factors
            twiddle_real[m] = c_del*twiddle_real[m-1] - s_del*twiddle_img[m-1];
            twiddle_img[m] = c_del*twiddle_img[m-1] + s_del*twiddle_real[m-1];
        }
        
        test_time_1 = MPI_Wtime();
        for (int i=0; i<N; i++){
            
            int index = bitrev_index[i];
            array_real[i]= sample[index][0];
            array_img[i]= sample[index][1];
        }
        test_time_1= MPI_Wtime()- test_time_1;
        delete[] bitrev_index;
    }
    
    initial_time_2 = MPI_Wtime();
    MPI_Bcast(&twiddle_real, N/2, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&twiddle_img, N/2, MPI_FLOAT, 0, MPI_COMM_WORLD);


    //scatter the data to all process
    MPI_Scatter(&array_real, chunksize, MPI_FLOAT, &sub_array_real,chunksize,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(&array_img, chunksize, MPI_FLOAT, &sub_array_img,chunksize,MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    initial_time_3 = MPI_Wtime();
    for (k=0; k<log2(N/size_Of_Cluster); k++){//only go through the stages where the number of blocks is more than 
        int num_blks_stage_proc = chunksize/pow(2,k+1);
        int num_bf_block = pow(2,k);
        int bf_span = pow(2,k);
        int blk_step = pow(2, k+1);
        int bf_step=  1;
        int twiddle_index_step = N/pow(2,k+1);
        //create a buffer to hold the chunk of real and imginary data
        
        //perform the serial fft for each block in each processor. 
        for (int m = 0; m<num_blks_stage_proc; m++){ 
            for (int n = 0; n<num_bf_block; n++){
                //cout<< "value n = "<< n<< endl;
                int blk_pntr = m*blk_step;
                int twiddle_index = n*twiddle_index_step;
                float tf_real = twiddle_real[twiddle_index]; // calling the function everytime is not good 
                float tf_img = -twiddle_img[twiddle_index]; // add the twiddle factor function 
                int input_1 = blk_pntr + n*bf_step;
                int input_2 = input_1 + bf_span;
                float bf_input_1_real =  sub_array_real[input_1];
                float bf_input_1_img =  sub_array_img[input_1];
                float bf_input_2_real = sub_array_real[input_2]*tf_real+ sub_array_img[input_2]*(-tf_img);
                float bf_input_2_img = sub_array_img[input_2]*tf_real - sub_array_real[input_2]*(-tf_img);
                float bf_out_1_real = bf_input_1_real + bf_input_2_real;
                float bf_out_1_img = bf_input_1_img + bf_input_2_img;
                float bf_out_2_real = bf_input_1_real - bf_input_2_real;
                float bf_out_2_img = bf_input_1_img - bf_input_2_img;
                sub_array_real[input_1] =  bf_out_1_real;
                sub_array_img[input_1] =  bf_out_1_img;
                sub_array_real[input_2] = bf_out_2_real;
                sub_array_img[input_2] = bf_out_2_img;
                //cout<<sub_array_real[input_2]<<endl;
            }

        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    initial_time_4 = MPI_Wtime();
    MPI_Gather(&sub_array_real, chunksize, MPI_FLOAT, &array_real, chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sub_array_img, chunksize, MPI_FLOAT, &array_img, chunksize, MPI_FLOAT, 0, MPI_COMM_WORLD);
    initial_time_5 = MPI_Wtime();

    //delete[] bitrev_index;

    if (process_Rank == 0){
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
                    float tf_real = twiddle_real[twiddle_index]; // calling the function everytime is not good 
                    float tf_img = -twiddle_img[twiddle_index]; // add the twiddle factor function 
                    int input_1 = blk_pntr + n*bf_step;
                    int input_2 = input_1 + bf_span;
                    float bf_input_1_real =  array_real[input_1];
                    float bf_input_1_img =  array_img[input_1];
                    float bf_input_2_real = array_real[input_2]*tf_real+ array_img[input_2]*(-tf_img);
                    float bf_input_2_img = array_img[input_2]*tf_real - array_real[input_2]*(-tf_img);
                    float bf_out_1_real = bf_input_1_real + bf_input_2_real;
                    float bf_out_1_img = bf_input_1_img + bf_input_2_img;
                    float bf_out_2_real = bf_input_1_real - bf_input_2_real;
                    float bf_out_2_img = bf_input_1_img - bf_input_2_img;
                    array_real[input_1] =  bf_out_1_real;
                    array_img[input_1] =  bf_out_1_img;
                    array_real[input_2] = bf_out_2_real;
                    array_img[input_2] = bf_out_2_img;
                }
            }
        }
        //float initial_time_6 = MPI_Wtime();
        //float computation_time = (initial_time_2-initial_time_1)+(initial_time_4-initial_time_3)+(initial_time_6-initial_time_5);
        //float communication_time = (initial_time_3-initial_time_2)+para_time;
        //float computation_time =  (initial_time_6 - initial_time_1)-communication_time;
        //float total_time = (initial_time_6-initial_time_1);
        //float para_time = final_para_time - initial_time;
        //float seq_time = final_time - final_para_time;
        //cout<<"computation time: "<< computation_time<<endl;
        //cout<<"communication time: "<< communication_time<<endl;
        //cout<<"total time: "<< total_time<<endl;
        initial_time_6 =  MPI_Wtime();
    }
    
    double total_time = initial_time_6-initial_time_1;
    double communication_time = (initial_time_3-initial_time_2)+(initial_time_5-initial_time_4)+(initial_time_8-initial_time_7);
    double parallel_time = (initial_time_4-initial_time_3);
    double serial_time = (initial_time_6-initial_time_5);
    double set_up_time = (initial_time_7-initial_time_1)+(initial_time_2-initial_time_8);
    double compute_time =  (initial_time_7-initial_time_1)+(initial_time_2-initial_time_8)+(initial_time_4-initial_time_3)+(initial_time_6-initial_time_5);
    if (process_Rank==0){
        cout<<"total time:          "<< total_time<<endl;
        cout<<"communication time:  "<< communication_time<<endl;
        cout<<"compute time:        "<<compute_time<<endl;
        cout<<"para compute time:   "<<parallel_time<<endl;
        cout<<"serial compute time: "<<serial_time<<endl;
        cout<<"setup compute time:  "<<set_up_time<<endl;
        cout<<"test time:           "<<test_time_1<<endl;

    }
    


    MPI_Finalize();
    return 0;
}