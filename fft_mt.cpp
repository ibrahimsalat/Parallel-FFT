#include<iostream>
#include <array>
#include <vector>
#include<cmath>
#include<complex>
#include<math.h> 
#include<omp.h> 
using namespace std;
/*The bit Reversed function that calculates the bit reversed number of a given number n*/
unsigned int reverseBits(unsigned int n, int stages)
{   
	unsigned int rev = 0;
    for (int i=0; i<stages; i++){
        rev <<= 1;
        if ((n & 1)==1){
            rev++;
        }
        n>>=1;

    }
    return rev;
}
/*The main Cooley-Turkey FFT that takes an input sample and outputs the FFT of the sample*/
vector<vector<double>> fft(int stages, vector<vector<double>> sample, int num_of_threads){
    /*Initialization of constansts and arrays*/
    int N = pow(2,stages);
    double pi = 3.1415926535897932384626;
    vector<double> omega_cos (N/2, 0);
    vector<double> omega_sin (N/2, 0);
    vector<vector<double>> tf (N/2, vector<double>(2,0));
    
    /*calculates the twiddle factors and stores them in arrays of real and imaginary*/
    double del_theta = 2*pi/N;
    double c_del = cos(del_theta);
    double s_del = sin(del_theta);
    omp_set_num_threads(num_of_threads);

    tf[0][0] = 1;
    tf[0][1] = 0;
    #pragma omp parallel for
        for (int m =1; m<N/2; m++){ //twiddle factors
            tf[m][0] = c_del*tf[m-1][0] - s_del*tf[m-1][1];
            tf[m][1] = c_del*tf[m-1][1] + s_del*tf[m-1][0];

    }
    
    for (int k=0; k<stages; k++){//loops for stages 
        /*calculation of FFT parameters*/
        int num_blks_stages= N/(pow(2,k+1));
        int num_bf_block = pow(2,k);        
        int bf_span = pow(2,k);
        int blk_step = pow(2, k+1);
        int bf_step=  1;
        int twiddle_index_step = N/pow(2,k+1);
        #pragma omp parallel for collapse(2)
            for (int m = 0; m<num_blks_stages; m++){ //loops over blocks in a stage
                for (int n = 0; n<num_bf_block; n++){//loops over butterfly per block
                    int blk_pntr = m*blk_step;
                    int twiddle_index = n*twiddle_index_step;
                    double tf_real = tf[twiddle_index][0]; 
                    double tf_img = -tf[twiddle_index][1]; // add the twiddle factor
                    double input_1 = blk_pntr + n*bf_step;
                    double input_2 = input_1 + bf_span;
                    double bf_input_1_real =  sample[input_1][0];
                    double bf_input_1_img =  sample[input_1][1];
                    double bf_input_2_real = sample[input_2][0]*tf_real+ sample[input_2][1]*(-tf_img);
                    double bf_input_2_img = sample[input_2][1]*tf_real - sample[input_2][0]*(-tf_img);
                    double bf_out_1_real = bf_input_1_real + bf_input_2_real;
                    double bf_out_1_img = bf_input_1_img + bf_input_2_img;
                    double bf_out_2_real = bf_input_1_real - bf_input_2_real;
                    double bf_out_2_img = bf_input_1_img - bf_input_2_img;
                    sample[input_1][0] =  bf_out_1_real;
                    sample[input_1][1] =  bf_out_1_img;
                    sample[input_2][0] = bf_out_2_real;
                    sample[input_2][1] = bf_out_2_img;
                }
            }
    }
    return sample;
}

//The main method
int main(int argc, char *argv[]){
    int stages = atoi(argv[1]);
    int num_of_threads = atoi(argv[2]);
    double pi = 3.1415926535897932384626;
    double initial, final;
    int N = pow(2,stages);
    vector<vector<double>> sample (N, vector<double>(2,0));
    vector<double> sample_real (N, 0);
    vector<double> sample_complex (N, 0);
    for (int i=0; i<N; i++){
        sample[i][0]=(sin((2*pi*i/N)));

    }
    omp_set_num_threads(num_of_threads);
    initial = omp_get_wtime();
    for (int i = 0; i<N; i++){
        sample_real[i]= sample[i][0];
        sample_complex[i]= sample[i][1];
        
    }
    /*prepares the sample in a bit reversed order*/
    #pragma omp parallel for
        for (int i=0; i<N; i++){
            int index = reverseBits(i, stages);
            sample[i][0]= sample_real[index];
            sample[i][1]= sample_complex[index];
        }

    
    fft(stages, sample, num_of_threads);
    final = omp_get_wtime();
    double time = final -initial;
    cout<< "time taken with "<< num_of_threads<< " threads:"<<endl;
    cout<<time<<endl;
    return 0;
}