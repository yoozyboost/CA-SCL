#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stack>
#include <set>
#include "PolarCode.h"

using namespace std;

#define CODE_RATE 0.5
#define N 1024
#define k 496
#define r 16
// #define list_size 4
#define n (int)log2(N)

#define PI M_PI

#define SNR_MIN -3				
#define SNR_MAX 3
#define EbN0_MIN -3
#define EbN0_MAX -1.5				
#define SNR_INTERVAL 1
#define EbN0_INTERVAL 0.25
#define TRIAL_MAX 200000
#define TRIAL_MIN 1000			
#define TRIAL_INTERVAL 100		
#define ERROR_MAX 10000

int main(int argc, char *argv[]){

    long int trial;
	double snr;//signal noise ratio	
    double ebn0;				    
    double BER; //bit error rate
    double BLER; //block error rate

    double EbN0dB;
    double EbN0;
    double EsN0;

    double Eb;
    double Es;
    double N0;

    double sigma_w; //variance of AWGN

    int error_count;
    int block_error_count;

    int construction_mode;
    double design_snr = 0.0;
    int list_size;
    int crc_mode;


    srand((unsigned) time(NULL));

    construction_mode = atoi(argv[1]);
    design_snr = atof(argv[2]);
    list_size = atof(argv[3]);
    crc_mode = atof(argv[4]);

    int flag;
    int error_count_to_print;    

    char *sFile1 = new char[200];
    char *sFile2 = new char[200];

    if(crc_mode == 0){
    
        if(construction_mode == 0){
            sprintf(sFile1, "BPSK_polar_BER_N=%d_L=%d_BTC_designsnr=%.2f.txt",N,list_size,design_snr);
            sprintf(sFile2, "BPSK_polar_BLER_N=%d_L=%d_BTC_designsnr=%.2f.txt",N,list_size,design_snr);
        }
        else if(construction_mode == 1){
            sprintf(sFile1, "BPSK_polar_BER_N=%d_IGA_designsnr=%.2f.txt",N,design_snr);
            sprintf(sFile2, "BPSK_polar_BLER_N=%d_IGA_designsnr=%.2f.txt",N,design_snr);
        }
    // else if(construction_mode == 2){

    // sprintf(sFile1, "BPSK_polar_BER_N=%d_RCA_designsnr=%.2f.txt",N,design_snr);
    // sprintf(sFile2, "BPSK_polar_BLER_N=%d_RCA_designsnr=%.2f.txt",N,design_snr);

    // }
        else{
            printf("mode is invalid");
            exit(0);
        }
    }
    else{
        if(construction_mode == 0){
            sprintf(sFile1, "BPSK_polar_BER_N=%d_L=%d_CRC-%d_BTC_designsnr=%.2f.txt",N,list_size,crc_mode,design_snr);
            sprintf(sFile2, "BPSK_polar_BLER_N=%d_L=%d_CRC-%d_BTC_designsnr=%.2f.txt",N,list_size,crc_mode,design_snr);
        
        }
        else if(construction_mode == 1){
            sprintf(sFile1, "BPSK_polar_BER_N=%d_IGA_designsnr=%.2f.txt",N,design_snr);
            sprintf(sFile2, "BPSK_polar_BLER_N=%d_IGA_designsnr=%.2f.txt",N,design_snr);
        }
    // else if(construction_mode == 2){

    // sprintf(sFile1, "BPSK_polar_BER_N=%d_RCA_designsnr=%.2f.txt",N,design_snr);
    // sprintf(sFile2, "BPSK_polar_BLER_N=%d_RCA_designsnr=%.2f.txt",N,design_snr);

    // }
        else{
            printf("mode is invalid");
            exit(0);
        }

    }

    FILE *fp1,*fp2;
    fp1 = fopen(sFile1, "w");
    fp2 = fopen(sFile2, "w");

    PolarCode polarcode(N,k,construction_mode,design_snr,crc_mode);


    for( ebn0 = EbN0_MIN; ebn0 <= EbN0_MAX ; ebn0 += EbN0_INTERVAL ){

        sigma_w = sqrt(1) / sqrt( ( 2.0  *  (CODE_RATE) *
                            pow(10.0, ( (double)ebn0 / 10.0 ))));

        // Es = 1.0;
        // Eb = Es / (CODE_RATE); 
        // EbN0 = pow(10.0,0.1 * ebn0);
        // N0 = Eb / EbN0;
        // sigma_w = sqrt(N0);

		error_count = 0;
        block_error_count = 0;

        for( trial = 0 ; trial < TRIAL_MAX; trial++ ){


            std::vector<uint8_t> info_bits = polarcode.dataGenerator(k);
            
            std::vector<uint8_t> code_bits = polarcode.encode(info_bits);

            std::vector<int8_t> data2 = polarcode.BPSKmodulation(code_bits);
            std::vector<double> received_signal = polarcode.AWGNchannel(data2,sigma_w);
            // std::vector<uint8_t> received_signal = polarcode.BPSKdemodulation(data3);
 

            std::vector<double> p0(N), p1(N);
            std::vector<double> llr(N);

            for (uint16_t i = 0; i < N; ++i) {
                // p1.at(i) = exp(-(received_signal.at(i) + 1.0)*(received_signal.at(i) + 1.0 )/(2.0*sigma_w*sigma_w))/sqrt(2.0 * M_PI * sigma_w * sigma_w);
                // p0.at(i) = exp(-(received_signal.at(i) - 1.0)*(received_signal.at(i) - 1.0)/(2.0*sigma_w*sigma_w))/sqrt(2.0 * M_PI * sigma_w * sigma_w);
                llr.at(i) = 4 * received_signal.at(i) / (sigma_w * sigma_w);
            }

            // std::vector<uint8_t> decoded_info_bits = polarcode.decode_scl_p1(p1, p0, list_size);   
            std::vector<uint8_t> decoded_info_bits = polarcode.decode_scl_llr(llr, list_size);   

            error_count += polarcode.errorCount(info_bits,decoded_info_bits);
            block_error_count += polarcode.blockErrorCount(info_bits,decoded_info_bits);

            if( ( trial > TRIAL_MIN ) && ( error_count > ERROR_MAX ) ){
				trial++; break;
			}

			// show the progress
			if( trial % TRIAL_INTERVAL == 0 ){
			fprintf( stderr, " %f %8ld    (%d)\r", ebn0, trial, error_count);
			}

        }//end of one trial


        long all_bits = k * trial;
        long all_block = trial;
        BER = (double)error_count/all_bits;
        BLER = (double)block_error_count/all_block;
        printf("EbN0=%g sigma=%g error_count = %d all = %ld BER=%g\n",ebn0,sigma_w,error_count,all_bits,BER);
        printf("EbN0=%g sigma=%g block error count = %d all block = %ld BLER=%g\n",ebn0,sigma_w,block_error_count,trial,BLER);
        fprintf(fp1,"%3f %g\n",ebn0,BER);
        fprintf(fp2,"%3f %g\n",ebn0,BLER);

    }

    fclose(fp1);
    fclose(fp2);

    printf("finished!");

    return 0;
}