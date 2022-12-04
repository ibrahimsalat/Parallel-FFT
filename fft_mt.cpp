#include<iostream>
#include <array>
#include <vector>
#include<cmath>
#include<complex>
#include<math.h> 
#include<omp.h> 
using namespace std;

int N = 10000000;
vector<vector<double>> sample (N, vector<double>(2,0));

//double N = sample.size();
double pi = 3.1415926535897932384626;
int stages = log2(N);
vector<double> omega_cos (N/2, 0);
vector<double> omega_sin (N/2, 0);
vector<double> c (N/2, 0);
vector<double> s (N/2, 0);
vector<double> k (stages, 0);
vector<double> sample_real (N, 0);
vector<double> sample_complex (N, 0);
vector<vector<double>> tf (N/2, vector<double>(2,0));
vector<vector<double>> w_index (N, vector<double>(stages,0));



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

vector<vector<double>> fft(){
    int n;
    double initial, final;
    double del_theta = 2*pi/N;
    double c_del = cos(del_theta);
    double s_del = sin(del_theta);
    int num_of_threads = 24;
    omp_set_num_threads(num_of_threads);
    initial = omp_get_wtime();

    tf[0][0] = 1;
    tf[0][1] = 0;
    #pragma omp parallel for
        for (int m =1; m<N/2; m++){ //twiddle factors
            tf[m][0] = c_del*tf[m-1][0] - s_del*tf[m-1][1];
            tf[m][1] = c_del*tf[m-1][1] + s_del*tf[m-1][0];
            //cout<< tf[m][0]<<endl;
    }
    
    for (int k=0; k<stages; k++){
        int num_blks_stages= N/(pow(2,k+1));
        int num_bf_block = pow(2,k);        
        int bf_span = pow(2,k);
        int blk_step = pow(2, k+1);
        int bf_step=  1;
        int twiddle_index_step = N/pow(2,k+1);
        #pragma omp parallel for collapse(2)
            for (int m = 0; m<num_blks_stages; m++){ 
                for (int n = 0; n<num_bf_block; n++){
                    //cout<< "value n = "<< n<< endl;
                    int blk_pntr = m*blk_step;
                    int twiddle_index = n*twiddle_index_step;
                    double tf_real = tf[twiddle_index][0]; // calling the function everytime is not good 
                    double tf_img = -tf[twiddle_index][1]; // add the twiddle factor function 
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
    final = omp_get_wtime();
    double time = final - initial;
    cout<< "time taken with "<< num_of_threads<< " threads:"<<endl;
    cout<< time<< endl;
    return sample;
}

//The main method
int main(){
    //double pi = 3.1415926535897932384626;
    for (int i=0; i<N; i++){
        sample[i][0]=(sin((2*pi*i/N)));

    }

    fft();
    //implement scheduling, miltithread the bit reverse. multithread twiddle factors
    //make graphs  
    return 0;
}