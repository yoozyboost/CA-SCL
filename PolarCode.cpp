#include "PolarCode.h"
#include <iostream>
#include <cmath>       /* log */
#include <sstream>      // std::stringstream
#include <fstream>
#include <iomanip>      // std::setprecision
#include <random>
#include <algorithm>
#include <chrono>

std::vector<uint8_t> PolarCode::dataGenerator(uint16_t data_length){

    std::vector<uint8_t> generated_data(data_length,0);

    for(uint16_t i=0;i<data_length;i++){
        generated_data.at(i) = (rand()%2);
    }

    return generated_data;
}


std::vector<double> PolarCode::AWGNchannel(std::vector<int8_t> send_data, double dispersion){

    uint16_t data_size = send_data.size();
    std::vector<double> sent_data(data_size);

    double u1,u2;
    double v1,v2;

    double sqrt2inv = 1.0/sqrt(2.0);

    for(uint16_t i = 0; i < data_size; i++ ){
		// generating random values
		do{
			u1 = (double) rand()/RAND_MAX;
		} while( u1 == 0 );
		u2 = (double) rand()/RAND_MAX;

		v1 = sqrt( -2.0 * log( u1 ) );
		v2 = 2.0 * M_PI * u2;

		sent_data.at(i) = send_data.at(i) +  v1 * dispersion * cos(v2) * sqrt2inv;
    }

    return sent_data;

}

std::vector<int8_t> PolarCode::BPSKmodulation(std::vector<uint8_t> data_to_modulate){
    uint16_t data_size = data_to_modulate.size();
    std::vector<int8_t> modulated_data(data_size);

    for(int16_t i=0;i<data_size;i++){
        modulated_data.at(i) = 1 - 2 * data_to_modulate.at(i);
    }

    return modulated_data;

}

std::vector<uint8_t> PolarCode::BPSKdemodulation(std::vector<double> data_to_demodulate){
    uint16_t data_size = data_to_demodulate.size();
    std::vector<uint8_t> received_data(data_size);

    for(int16_t i=0;i<data_size;i++){
        if(data_to_demodulate.at(i) > 0.0){
            received_data.at(i) = 0;
        }
        else{
            received_data.at(i) = 1;
        }
    }

    return received_data;

}

uint16_t PolarCode::errorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2){
    uint16_t data_size = data1.size();
    uint16_t error_count = 0;
    for(uint16_t i=0;i<data_size;i++){
        if(data1.at(i) != data2.at(i)){
            error_count++;
        }
    }
    return error_count;
}

uint16_t PolarCode::blockErrorCount(std::vector<uint8_t> data1,std::vector<uint8_t> data2){
    uint16_t data_size = data1.size();
    uint16_t error_count = 0;
    for(uint16_t i=0;i<data_size;i++){
        if(data1.at(i) != data2.at(i)){
            error_count = 1;
            continue;
        }
    }
    return error_count;
}

void PolarCode::crc_type(){

    if(_crc_mode == 16){
        _crc_size = 16;
        _crc_array.resize(17);
        std::vector<uint8_t> crc = {1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1};
        _crc_array = crc;
    }
    else if(_crc_mode == 8){
        _crc_size = 8;
        _crc_array.resize(9);
        std::vector<uint8_t> crc = {1,1,1,0,1,0,1,0,1};
        _crc_array = crc;
    }
    else if(_crc_mode == 6){
        _crc_size = 6;
        _crc_array.resize(7);
        std::vector<uint8_t> crc = {1,1,0,0,0,0,1};
        _crc_array = crc;
    }
    else if(_crc_mode == 0){
        _crc_size = 0;
    }
    else{
        printf("invalid CRC mode"); 
        exit(0);
    }
}


void PolarCode::create_bit_rev_order() {
    for (uint16_t i = 0; i < _N; ++i) {
        uint16_t to_be_reversed = i;
        _bit_rev_order.at(i) = (uint16_t) ((to_be_reversed & 1) << (_n - 1));
        for (uint8_t j = (uint8_t) (_n - 1); j; --j) {
            to_be_reversed >>= 1;
            _bit_rev_order.at(i) += (to_be_reversed & 1) << (j - 1);
        }
    }
}


std::vector<uint8_t> PolarCode::encode(std::vector<uint8_t> info_bits) {

    std::vector<uint8_t> info_bits_padded(_block_length, 0);
    std::vector<uint8_t> info_crc_bits_padded(_info_length + _crc_mode);
    std::vector<uint8_t> coded_bits(_block_length);
    std::vector<uint8_t> crc_pool(_info_length+_crc_size,0);
    
    for (uint16_t i = 0; i < _info_length; ++i) {
            info_crc_bits_padded.at(i) = info_bits.at(i);
            crc_pool.at(i) = info_bits.at(i);
    }

    if(_crc_size != 0){
        for(int i=0;i<_info_length;i++){
            if(crc_pool.at(i) != 0){
                for(int j=0;j < _crc_size+1;j++){
                    crc_pool.at(i+j) = crc_pool.at(i+j)^ _crc_array.at(j);
                }
            }
        }
        for (uint16_t i = _info_length; i < _info_length + _crc_size; ++i) {
            info_crc_bits_padded.at(i) = crc_pool.at(i);
        }
    }

    uint16_t info_crc_index = 0;

    for(int i = 0; i < _block_length; i++){
        if(_frozen_bits.at(i) == 0){
            info_bits_padded.at(i) = info_crc_bits_padded.at(info_crc_index);
            info_crc_index++;
        }
    }

    for (uint8_t iteration = 0; iteration < _n; ++iteration) {
        uint16_t  increment = (uint16_t) (1 << iteration);
        for (uint16_t j = 0; j < increment; j +=  1) {
            for (uint16_t i = 0; i < _N; i += 2 * increment) {
                info_bits_padded.at(i + j) = (uint8_t)((info_bits_padded.at(i + j) + info_bits_padded.at(i + j + increment)) % 2);
            }
        }
    }

    for (uint16_t i = 0; i < _block_length; ++i) {
        coded_bits.at(i) = info_bits_padded.at(_bit_rev_order.at(i));
    }

    return coded_bits;

}



std::vector<uint8_t> PolarCode::decode_scl_p1(std::vector<double> p1, std::vector<double> p0, uint16_t list_size) {

    _list_size = list_size;
    _llr_based_computation = false;

    initializeDataStructures();

    uint16_t  l = assignInitialPath();

    double * p_0 = getArrayPointer_P(0, l);

    for (uint16_t beta = 0; beta < _block_length; ++beta ) {
        p_0[2*beta] = (double) p0.at(beta);
        p_0[2*beta + 1] = (double) p1.at(beta);
    }

    return decode_scl();

}

std::vector<uint8_t> PolarCode::decode_scl_llr(std::vector<double> llr, uint16_t list_size) {

    _list_size = list_size;

    _llr_based_computation = true;

    initializeDataStructures();

    uint16_t  l = assignInitialPath();

    double * llr_0 = getArrayPointer_LLR(0, l);

    for (uint16_t beta = 0; beta < _block_length; ++beta ) {
        llr_0[beta] = llr.at(beta);
    }

    return decode_scl();

}

std::vector<uint8_t> PolarCode::decode_scl() {


    for (uint16_t phi = 0; phi < _block_length; ++phi ){

        if (_llr_based_computation )
            recursivelyCalcLLR(_n, phi);
        else
            recursivelyCalcP(_n, phi);


        if (_frozen_bits.at(phi) == 1)
            continuePaths_FrozenBit(phi);
        else
            continuePaths_UnfrozenBit(phi);

        if ((phi%2) == 1)
            recursivelyUpdateC(_n, phi);

    }
    uint16_t l = findMostProbablePath((bool) _crc_size);
    
    uint8_t * c_0 = _arrayPointer_Info.at(l);
    std::vector<uint8_t> decoded_info_bits(_info_length);

    std::vector<uint8_t> info_crc_bits(_info_length + _crc_mode);

    uint16_t info_crc_index = 0;

    for (uint16_t beta = 0; beta < _block_length; ++beta ){
        if(_frozen_bits.at(beta) == 0){
            decoded_info_bits.at(info_crc_index) = c_0[beta];
            info_crc_index++;
        }
        if(info_crc_index == _info_length){
            break;
        }
    }

    for (uint16_t s = 0; s < _list_size; ++s) {
        delete[] _arrayPointer_Info.at(s);
        for (uint16_t lambda = 0; lambda < _n + 1; ++lambda) {

            if (_llr_based_computation )
                delete[] _arrayPointer_LLR.at(lambda).at(s);
            else
                delete[] _arrayPointer_P.at(lambda).at(s);
            delete[] _arrayPointer_C.at(lambda).at(s);
        }
    }

    return decoded_info_bits;

}


void PolarCode::initializeDataStructures() {

    while (_inactivePathIndices.size()) {
        _inactivePathIndices.pop();
    };
    _activePath.resize(_list_size);

    if (_llr_based_computation) {
        _pathMetric_LLR.resize(_list_size);
        _arrayPointer_LLR.resize(_n + 1);
        for (int i = 0; i < _n + 1; ++i)
            _arrayPointer_LLR.at(i).resize(_list_size);
    }
    else {
        _arrayPointer_P.resize(_n + 1);
        for (int i = 0; i < _n + 1; ++i)
            _arrayPointer_P.at(i).resize(_list_size);
    }

    _arrayPointer_C.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayPointer_C.at(i).resize(_list_size);

    _arrayPointer_Info.resize(_list_size);

    _pathIndexToArrayIndex.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _pathIndexToArrayIndex.at(i).resize(_list_size);

    _inactiveArrayIndices.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i) {
        while (_inactiveArrayIndices.at(i).size()) {
            _inactiveArrayIndices.at(i).pop();
        };
    }

    _arrayReferenceCount.resize(_n + 1);
    for (int i = 0; i < _n + 1; ++i)
        _arrayReferenceCount.at(i).resize(_list_size);

    for (uint16_t s = 0; s < _list_size; ++s) {
        _arrayPointer_Info.at(s) = new uint8_t[_block_length]();
        for (uint16_t lambda = 0; lambda < _n + 1; ++lambda) {
            if (_llr_based_computation) {
                _arrayPointer_LLR.at(lambda).at(s) = new double[(1 << (_n - lambda))]();
            }
            else {
                _arrayPointer_P.at(lambda).at(s) = new double[2 * (1 << (_n - lambda))]();
            }
            _arrayPointer_C.at(lambda).at(s) = new uint8_t[2 * (1 << (_n - lambda))]();
            _arrayReferenceCount.at(lambda).at(s) = 0;
            _inactiveArrayIndices.at(lambda).push(s);
        }
    }

    for (uint16_t l = 0; l < _list_size; ++l) {
        _activePath.at(l) = 0;
        _inactivePathIndices.push(l);
        if (_llr_based_computation) {
            _pathMetric_LLR.at(l) = 0;
        }
    }
}

uint16_t PolarCode::assignInitialPath() {

    uint16_t  l = _inactivePathIndices.top();
    _inactivePathIndices.pop();
    _activePath.at(l) = 1;
    // Associate arrays with path index
    for (uint16_t lambda = 0; lambda < _n + 1; ++lambda) {
        uint16_t  s = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();
        _pathIndexToArrayIndex.at(lambda).at(l) = s;
        _arrayReferenceCount.at(lambda).at(s) = 1;
    }
    return l;
}

uint16_t PolarCode::clonePath(uint16_t l) {
    uint16_t l_p = _inactivePathIndices.top();
    _inactivePathIndices.pop();
    _activePath.at(l_p) = 1;

    if (_llr_based_computation)
        _pathMetric_LLR.at(l_p) = _pathMetric_LLR.at(l);

    for (uint16_t lambda = 0; lambda < _n + 1; ++lambda ) {
        uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
        _pathIndexToArrayIndex.at(lambda).at(l_p) = s;
        _arrayReferenceCount.at(lambda).at(s)++;
    }
    return l_p;
}

void PolarCode::killPath(uint16_t l) {
    _activePath.at(l) = 0;
    _inactivePathIndices.push(l);
    if (_llr_based_computation )
        _pathMetric_LLR.at(l) = 0;

    for (uint16_t lambda = 0; lambda < _n + 1; ++lambda ) {
        uint16_t s = _pathIndexToArrayIndex.at(lambda).at(l);
        _arrayReferenceCount.at(lambda).at(s)--;
        if (_arrayReferenceCount.at(lambda).at(s) == 0 ) {
            _inactiveArrayIndices.at(lambda).push(s);
        }
    }
}

double * PolarCode::getArrayPointer_P(uint16_t lambda, uint16_t  l) {
    uint16_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint16_t s_p;
    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {
        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_P.at(lambda).at(s_p));
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;
    }
    return _arrayPointer_P.at(lambda).at(s_p);
}

double * PolarCode::getArrayPointer_LLR(uint16_t lambda, uint16_t  l) {
    uint16_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint16_t s_p;
    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {
        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));
        std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) +  (1 << (_n - lambda)),  _arrayPointer_LLR.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;
    }
    return _arrayPointer_LLR.at(lambda).at(s_p);
}


uint8_t * PolarCode::getArrayPointer_C(uint16_t lambda, uint16_t  l) {
    uint16_t  s = _pathIndexToArrayIndex.at(lambda).at(l);
    uint16_t s_p;
    if (_arrayReferenceCount.at(lambda).at(s) == 1) {
        s_p = s;
    }
    else {

        s_p = _inactiveArrayIndices.at(lambda).top();
        _inactiveArrayIndices.at(lambda).pop();

        //copy
        if (_llr_based_computation )
            std::copy(_arrayPointer_LLR.at(lambda).at(s), _arrayPointer_LLR.at(lambda).at(s) +  (1 << (_n - lambda)),  _arrayPointer_LLR.at(lambda).at(s_p));
        else
            std::copy(_arrayPointer_P.at(lambda).at(s), _arrayPointer_P.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_P.at(lambda).at(s_p));

        std::copy(_arrayPointer_C.at(lambda).at(s), _arrayPointer_C.at(lambda).at(s) +  (1 << (_n - lambda + 1)),  _arrayPointer_C.at(lambda).at(s_p));

        _arrayReferenceCount.at(lambda).at(s)--;
        _arrayReferenceCount.at(lambda).at(s_p) = 1;
        _pathIndexToArrayIndex.at(lambda).at(l) = s_p;

    }
    return _arrayPointer_C.at(lambda).at(s_p);
}

void PolarCode::recursivelyCalcP(uint16_t lambda, uint16_t phi) {
    if ( lambda == 0 )
        return;
    uint16_t psi = phi >> 1;
    if ( (phi % 2) == 0)
        recursivelyCalcP(lambda -1, psi);

    double sigma = 0.0f;
    for (uint16_t l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        double * p_lambda = getArrayPointer_P(lambda, l);
        double * p_lambda_1 = getArrayPointer_P(lambda - 1, l);

        uint8_t * c_lambda = getArrayPointer_C(lambda, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            if ( (phi %2) == 0 ){
                p_lambda[2 * beta] = 0.5f * ( p_lambda_1[2*(2*beta)]*p_lambda_1[2*(2*beta+1)]
                                              + p_lambda_1[2*(2*beta) + 1]*p_lambda_1[2*(2*beta+1) + 1]);
                p_lambda[2 * beta + 1] = 0.5f * ( p_lambda_1[2*(2*beta) +1]*p_lambda_1[2*(2*beta+1)]
                                                  + p_lambda_1[2*(2*beta)]*p_lambda_1[2*(2*beta+1) + 1]);
            }
            else {
                uint8_t  u_p = c_lambda[2*beta];
                p_lambda[2 * beta] = 0.5f * p_lambda_1[2*(2*beta) + (u_p % 2)] *   p_lambda_1[2*(2*beta + 1)];
                p_lambda[2 * beta + 1] = 0.5f * p_lambda_1[2*(2*beta) + ((u_p+1) % 2)] *   p_lambda_1[2*(2*beta + 1) + 1];
            }
            sigma = std::max(sigma,  p_lambda[2 * beta]);
            sigma = std::max(sigma,  p_lambda[2 * beta + 1]);


        }
    }

    for (uint16_t l = 0; l < _list_size; ++l) {
        if (sigma == 0) // Typically happens because of undeflow
            break;
        if (_activePath.at(l) == 0)
            continue;
        double *p_lambda = getArrayPointer_P(lambda, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            p_lambda[2 * beta] = p_lambda[2 * beta] / sigma;
            p_lambda[2 * beta + 1] = p_lambda[2 * beta + 1] / sigma;
        }
    }
}

void PolarCode::recursivelyCalcLLR(uint16_t lambda, uint16_t phi) {
    if ( lambda == 0 )
        return;
    uint16_t psi = phi >> 1;
    if ( (phi % 2) == 0)
        recursivelyCalcLLR(lambda -1, psi);

    for (uint16_t l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        double * llr_lambda = getArrayPointer_LLR(lambda, l);
        double * llr_lambda_1 = getArrayPointer_LLR(lambda - 1, l);

        uint8_t * c_lambda = getArrayPointer_C(lambda, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            if ( (phi %2) == 0 ){
                if (40 > std::max(std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]))){
                    llr_lambda[beta] = std::log ( (exp(llr_lambda_1[2 * beta] + llr_lambda_1[2 * beta + 1]) + 1) /
                                                  (exp(llr_lambda_1[2*beta]) + exp(llr_lambda_1[2*beta+1])));
                }
                else {
                    llr_lambda[beta] = (double)  ((llr_lambda_1[2 * beta] < 0) ? -1 : (llr_lambda_1[2 * beta] > 0)) *
                                       ((llr_lambda_1[2 * beta + 1] < 0) ? -1 : (llr_lambda_1[2 * beta + 1] > 0)) *
                                       std::min( std::abs(llr_lambda_1[2 * beta]), std::abs(llr_lambda_1[2 * beta + 1]));
                }
            }
            else {
                uint8_t  u_p = c_lambda[2*beta];
                llr_lambda[beta] = (1 - 2 * u_p) * llr_lambda_1[2*beta] + llr_lambda_1[2*beta + 1];
            }

        }
    }
}

void PolarCode::recursivelyUpdateC(uint16_t lambda, uint16_t phi) {

    uint16_t psi = phi >> 1;
    for (uint16_t l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        uint8_t *c_lambda = getArrayPointer_C(lambda, l);
        uint8_t *c_lambda_1 = getArrayPointer_C(lambda - 1, l);
        for (uint16_t beta = 0; beta < (1 << (_n - lambda)); ++beta) {
            c_lambda_1[2 * (2 * beta) + (psi % 2)] = (uint8_t) ((c_lambda[2 * beta] + c_lambda[2 * beta + 1]) % 2);
            c_lambda_1[2 * (2 * beta + 1) + (psi % 2)] = c_lambda[2 * beta + 1];
        }
    }
    if ( (psi % 2) == 1)
        recursivelyUpdateC((uint16_t) (lambda - 1), psi);

}

void PolarCode::continuePaths_FrozenBit(uint16_t phi) {
    for (uint16_t l = 0; l < _list_size; ++ l) {
        if (_activePath.at(l) == 0)
            continue;
        uint8_t  * c_m = getArrayPointer_C(_n, l);
        c_m[(phi % 2)] = 0; // frozen value assumed to be zero
        if (_llr_based_computation) {
            double *llr_p = getArrayPointer_LLR(_n, l);
            _pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
        }
        _arrayPointer_Info.at(l)[phi] = 0;
    }
}

void PolarCode::continuePaths_UnfrozenBit(uint16_t phi) {

    std::vector<double>  probForks((unsigned long) (2 * _list_size));
    std::vector<double> probabilities;
    std::vector<uint8_t>  contForks((unsigned long) (2 * _list_size));


    uint16_t  i = 0;
    for (unsigned l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0) {
            probForks.at(2 * l) = NAN;
            probForks.at(2 * l + 1) = NAN;
        }
        else {
            if (_llr_based_computation ) {
                double *llr_p = getArrayPointer_LLR(_n, l);
                probForks.at(2 * l) =  - (_pathMetric_LLR.at(l) + log(1 + exp(-llr_p[0])));
                probForks.at(2 * l + 1) = -  (_pathMetric_LLR.at(l) + log(1 + exp(llr_p[0])));
            }
            else {
                double *p_m = getArrayPointer_P(_n, l);
                probForks.at(2 * l) = p_m[0];
                probForks.at(2 * l + 1) = p_m[1];
            }

            probabilities.push_back(probForks.at(2 * l));
            probabilities.push_back(probForks.at(2 * l +1));

            i++;
        }
    }

    uint16_t  rho = _list_size;
    if ( (2*i) < _list_size)
        rho = (uint16_t) 2 * i;

    for (uint8_t l = 0; l < 2 * _list_size; ++l) {
        contForks.at(l) = 0;
    }
    std::sort(probabilities.begin(), probabilities.end(), std::greater<double>());

    double threshold = probabilities.at((unsigned long) (rho - 1));
    uint16_t num_paths_continued = 0;

    for (uint8_t l = 0; l < 2 * _list_size; ++l) {
        if (probForks.at(l) > threshold) {
            contForks.at(l) = 1;
            num_paths_continued++;
        }
        if (num_paths_continued == rho) {
            break;
        }
    }

    if  ( num_paths_continued < rho ) {
        for (uint8_t l = 0; l < 2 * _list_size; ++l) {
            if (probForks.at(l) == threshold) {
                contForks.at(l) = 1;
                num_paths_continued++;
            }
            if (num_paths_continued == rho) {
                break;
            }
        }
    }

    for (unsigned l = 0; l < _list_size; ++l) {
        if (_activePath.at(l) == 0)
            continue;
        if ( contForks.at(2 * l)== 0 && contForks.at(2 * l + 1) == 0 )
            killPath(l);
    }

    for (unsigned l = 0; l < _list_size; ++l) {
        if ( contForks.at(2 * l) == 0 && contForks.at(2 * l + 1) == 0 )
            continue;
        uint8_t * c_m = getArrayPointer_C(_n, l);

        if ( contForks.at(2 * l) == 1 && contForks.at(2 * l + 1) == 1 ) {

            c_m[(phi%2)] = 0;
            uint16_t l_p = clonePath(l);
            c_m = getArrayPointer_C(_n, l_p);
            c_m[(phi%2)] = 1;

            std::copy(_arrayPointer_Info.at(l), _arrayPointer_Info.at(l) +  phi,  _arrayPointer_Info.at(l_p));
            _arrayPointer_Info.at(l)[phi] = 0;
            _arrayPointer_Info.at(l_p)[phi] = 1;

            if (_llr_based_computation ) {
                double *llr_p = getArrayPointer_LLR(_n, l);
                _pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
                llr_p = getArrayPointer_LLR(_n, l_p);
                _pathMetric_LLR.at(l_p) += log(1 + exp(llr_p[0]));
            }

        }
        else {
            if ( contForks.at(2 * l) == 1) {
                c_m[(phi%2)] = 0;
                _arrayPointer_Info.at(l)[phi] = 0;

                if (_llr_based_computation ) {
                    double *llr_p = getArrayPointer_LLR(_n, l);
                    _pathMetric_LLR.at(l) += log(1 + exp(-llr_p[0]));
                }
            }
            else {
                c_m[(phi%2)] = 1;
                _arrayPointer_Info.at(l)[phi] = 1;
                if (_llr_based_computation ) {
                    double *llr_p = getArrayPointer_LLR(_n, l);
                    _pathMetric_LLR.at(l) += log(1 + exp(llr_p[0]));
                }
            }
        }
    }

}

uint16_t PolarCode::findMostProbablePath(bool check_crc) {

    uint16_t  l_p = 0;
    double p_p1 = 0;
    double p_llr = std::numeric_limits<double>::max();
    bool path_with_crc_pass = false;
    for (uint16_t l = 0; l < _list_size; ++l) {

        if (_activePath.at(l) == 0)
            continue;

        if ( (check_crc) && (! crc_check(_arrayPointer_Info.at(l))))
            continue;

        path_with_crc_pass = true;

        if (_llr_based_computation) {
            if (_pathMetric_LLR.at(l) < p_llr ) {
                p_llr = _pathMetric_LLR.at(l);
                l_p  = l;
            }
        }
        else {
            uint8_t * c_m = getArrayPointer_C(_n, l);
            double * p_m = getArrayPointer_P(_n, l);
            if ( p_p1 < p_m[c_m[1]]) {
                l_p = l;
                p_p1 = p_m[c_m[1]];
            }
        }
    }

    if ( path_with_crc_pass)
        return l_p;
    else
        return findMostProbablePath(false);
}

double PolarCode::calcXi(double gamma){
    double xi;

    double threshold0 = 0.2;
    double threshold1 = 0.7;
    double threshold2 = 10.0;

    double a = -0.4527;
    double b = 0.0218;
    double c = 0.86;
    double k0 = 8.554;

    double a0 = -0.002706;
    double a1 = -0.476711;
    double a2 = 0.0512;

    if(gamma <= threshold0){
        xi = -0.5 * gamma + 0.125 * pow(gamma,2.0) - 0.125 * pow(gamma,3.0);
    }
    else if(gamma <= threshold1){
        xi = a0 + a1 * gamma + a2 * pow(gamma,2.0);
    }
    else if(gamma <= threshold2){
        xi = a * pow(gamma,c) + b;
    }
    else{
        xi = -0.25 * gamma + 0.5 * log(M_PI) - 0.5 * log(gamma) + log(1 - 0.25*(M_PI*M_PI)/gamma + k0 / (gamma * gamma)) ;
    }

    return xi;
}

double PolarCode::calcXiInverse(double z){

    double xi_inverse;
    double eps;


    double threshold0 = 0.2;
    double threshold1 = 0.7;
    double threshold2 = 10.0;

    double a = -0.4527;
    double b = 0.0218;
    double c = 0.86;

    double a0 = -0.002706;
    double a1 = -0.476711;
    double a2 = 0.0512;

    double Z0 = -0.5 * threshold0+ 0.125 * pow(threshold0,2.0) - 0.125 * pow(threshold0,3.0);
    double Z1 = a0 + a1 * threshold1 + a2 * pow(threshold1,2.0);
    double Z2 = a * pow(threshold2,c) + b;


    if(z >= Z0){
        xi_inverse = -2.0 * z * pow(z,2.0)  + pow(z,3.0);
    }
    else if(z >= Z1){
        xi_inverse = (-a1 - sqrt(a1 * a1 - 4.0 * a2 * (a0 - z)))/(2.0 * a2);
    }
    else if(z >=Z2){
        xi_inverse = pow((z-b)/(a) , 1.0/c);
    }
    else{

        eps = 1.0e-12;
        if( z < -2000.0 )
            eps = 1.0e-10;

        xi_inverse = bisection( z, 0.0, -4.0 * z, eps );
        if( xi_inverse < -1.0 ) {
            printf("retry\n");
            eps = 1.0e-9;
            xi_inverse = bisection( z, 0.0, -4.0 * z, eps );      
        }
        if( xi_inverse < -1.0 ) {
            printf("retry again\n");
            eps = 1.0e-8;
            xi_inverse = bisection( z, 0.0, -4.0 * z, eps );      
        }
        if( xi_inverse < -1.0 ) {
            printf("Give up\n");
            exit(0);
        }
    }

    return xi_inverse;   
}

double PolarCode::bisection(double y, double a, double b, double eps ){
    double s, diff;
    int i=0;
    
    while (!(fabs(a-b)<eps)){
    if( i > 10000 ) {
        diff = fabs( a- b );
        printf("failed to converge., y = %f, a = %f, b = %f, diff = %f(%e)\n", y, a, b, diff, diff );
        return - 100.0;
    }
    i++;

    s = (a+b)/2.0;
    if( ff(s,y) * ff(a,y)<0.0) b=s;
    else a = s;
    };
    return s;    
}

double PolarCode::ff(double x, double y)
{
  double f;
  double k0;


  k0 = 8.554;
  f =  - 0.25 * x - 0.5 * log( x ) + 0.5 * log( M_PI ) + log( 1.0 - M_PI * M_PI / (4.0 * x )  + k0 / ( x * x ) ) - y;

  return f;
}

double PolarCode::calcGamma(double u){
    double gamma;
    double z;

    double threshold0 = 0.2;
    double threshold1 = 0.7;
    double threshold2 = 10.0;

    double a = -0.4527;
    double b = 0.0218;
    double c = 0.86;
    double k0 = 8.554;


    if(u <= threshold0){
        gamma = 0.5 * pow(u,2.0) - 0.5 * pow(u,2.0) + 2.0/3.0 * pow(u,4.0);    
    }
    else{
        z = calcXi(u);
        gamma = calcXiInverse(z + log(2.0-exp(z)));
    }

    return gamma;    
}

void PolarCode::initialize_frozen_bits() {
    std::vector<double> channel_vec(_block_length);

    if(_construction_mode == 0){
        printf("construction mode is BTC\n");

        channel_vec.at(0) = exp(-pow(10.0,0.1 * _design_snr));
        for(int j = 0; j < _n; j ++ ){
        int u = ( 1 << (j+1) );

            for(int t = 0; t < (u/2); t++ ){
                double T = channel_vec.at(t);
                channel_vec.at(t) = 2.0 * T - T * T;
                channel_vec.at(t+u/2) = T * T;
            }
        }

        _channel_order_descending.resize(_block_length);
        std::size_t n_t(0);
        std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
        std::sort(  std::begin(_channel_order_descending),
                    std::end(_channel_order_descending),
                    [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] < channel_vec[_bit_rev_order.at(i2)]; } );
    }
    else if(_construction_mode == 1){

        printf("construction mode is IGA\n");
        channel_vec.at(0) = 4.0 * pow(10.0,0.1 * _design_snr);

        uint16_t J;
        double u;

        for(int i=1;i<=_n;i++){
            uint16_t J = (1 << i);

            for(int j=0;j<=(int)(J/2)-1;j++){
                u = channel_vec.at(j);
                channel_vec.at(j) = calcGamma(u);
                channel_vec.at(j + (int)J/2) = 2.0 * u;
            }
        }

        _channel_order_descending.resize(_block_length);
        std::size_t n_t(0);
        std::generate(std::begin(_channel_order_descending), std::end(_channel_order_descending), [&]{ return n_t++; });
        std::sort(  std::begin(_channel_order_descending),
                    std::end(_channel_order_descending),
                    [&](int i1, int i2) { return channel_vec[_bit_rev_order.at(i1)] > channel_vec[_bit_rev_order.at(i2)]; } );
     
    }
    else{
        exit(0);
    }

    uint16_t  effective_info_length = _info_length + _crc_size;

    for (uint16_t i = 0; i < effective_info_length; ++i) {
        _frozen_bits.at(_channel_order_descending.at(i)) = 0;
    }
    for (uint16_t i = effective_info_length; i < _block_length; ++i) {
        _frozen_bits.at(_channel_order_descending.at((i))) = 1;
    }

}

bool PolarCode::crc_check(uint8_t * info_bit_padded) {
    std::vector<uint8_t> crc_pool(_info_length+_crc_size);
    bool crc_pass = true;

    uint16_t info_crc_index = 0;
    for (uint16_t i = 0; i < _block_length; ++i) {
        if(_frozen_bits.at(i) == 0){
            crc_pool.at(info_crc_index) = info_bit_padded[info_crc_index];
            info_crc_index++;
        }
    }

    for(int i=0;i<_info_length;i++){
        if(crc_pool.at(i) != 0){
            for(int j=0;j < _crc_size+1;j++){
                crc_pool.at(i+j) = crc_pool.at(i+j)^ _crc_array.at(j);
            }
        }
    }

    for (uint16_t i = _info_length; i < _info_length + _crc_size; ++i) {
        if(crc_pool.at(i) != 0){
            crc_pass = false;
            break;
        }
    }

    return crc_pass;
}
