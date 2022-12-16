#ifndef POLARC_POLARCODE_H
#define POLARC_POLARCODE_H


#include <cstdint>
#include <vector>
#include <math.h>
#include <stack>          // std::stack

class PolarCode{

    public:

        PolarCode(uint16_t N,uint16_t K,uint16_t construction_mode,double design_snr,int crc_mode):
            _N(N), _K(K),_design_snr(design_snr),_block_length(N),_info_length(K),_crc_mode(crc_mode),_construction_mode(construction_mode)
        {
            _n = (uint16_t)log2(N);
            _bit_rev_order.resize(_N);
            _frozen_bits.resize(_N);
            _channel_order_descending.resize(_N);

            crc_type();
            create_bit_rev_order();
            initialize_frozen_bits();
            
            // getFrozenBitBTC(_design_snr);

        }

        //basic function for simulation
        std::vector<uint8_t> dataGenerator(uint16_t data_length);
        std::vector<int8_t> BPSKmodulation(std::vector<uint8_t> data_to_modulate);
        std::vector<double> AWGNchannel(std::vector<int8_t> send_data,double dispersion);
        std::vector<uint8_t> BPSKdemodulation(std::vector<double> data_to_demodulate);
        uint16_t errorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2);
        uint16_t blockErrorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2);
        
        //another function
        void quickSort(std::vector<double>& number,std::vector<uint16_t>& index,uint16_t start,uint16_t end);

        //polar code
        std::vector<uint8_t> encode(std::vector<uint8_t> info_bits);
        std::vector<uint8_t> decode_scl_p1(std::vector<double> p1, std::vector<double> p0, uint16_t list_size);
        std::vector<uint8_t> decode_scl_llr(std::vector<double> llr, uint16_t list_size);

    private:

        uint16_t _N;
        uint16_t _K;
        uint16_t _n;
        uint16_t _block_length;
        uint16_t _info_length;
        uint16_t _crc_size;
        double _design_snr;

        uint16_t _construction_mode;
        uint16_t _crc_mode;

        std::vector<uint8_t> _crc_array;

        std::vector<uint16_t> _bit_rev_order;
        std::vector<uint16_t> _channel_order_descending;
        std::vector<uint8_t> _frozen_bits;

        std::vector<uint8_t> _sending_bits;
        std::vector<uint8_t> _received_bits;

        void crc_type();
        void initialize_frozen_bits();
        void create_bit_rev_order();

        std::vector<uint8_t> decode_scl();
        bool _llr_based_computation;

        std::vector<std::vector<double *>> _arrayPointer_LLR;
        std::vector<double> _pathMetric_LLR;

        uint16_t _list_size;

        std::stack<uint16_t> _inactivePathIndices;
        std::vector<uint16_t > _activePath;
        std::vector<std::vector<double *>> _arrayPointer_P;
        std::vector<std::vector<uint8_t *>> _arrayPointer_C;
        std::vector<uint8_t *> _arrayPointer_Info;
        std::vector<std::vector<uint16_t>> _pathIndexToArrayIndex;
        std::vector<std::stack<uint16_t>> _inactiveArrayIndices;
        std::vector<std::vector<uint16_t>> _arrayReferenceCount;

        void initializeDataStructures();
        uint16_t assignInitialPath();
        uint16_t clonePath(uint16_t);
        void killPath(uint16_t l);

        double * getArrayPointer_P(uint16_t lambda, uint16_t  l);
        double * getArrayPointer_LLR(uint16_t lambda, uint16_t  l);
        uint8_t * getArrayPointer_C(uint16_t lambda, uint16_t  l);

        void recursivelyCalcP(uint16_t lambda, uint16_t phi);
        void recursivelyCalcLLR(uint16_t lambda, uint16_t phi);
        void recursivelyUpdateC(uint16_t lambda, uint16_t phi);

        void continuePaths_FrozenBit(uint16_t phi);
        void continuePaths_UnfrozenBit(uint16_t phi);

        uint16_t findMostProbablePath(bool check_crc);
        
        bool crc_check(uint8_t * info_bits_padded); 
        
        double calcGamma(double u);
        double calcXi(double gamma);
        double calcXiInverse(double z);
        double bisection(double y,double a,double b,double eps);
        double ff(double x,double y);

};


#endif //POLARC_POLARCODE_H