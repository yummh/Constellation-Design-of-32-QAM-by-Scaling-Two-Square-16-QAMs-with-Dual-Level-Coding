#include "../include/DMOFDM_AWGN_1.h" 


/////////////////////////////////////////////////////////
///////// parameter
/////////////////////////////////////////////////////////
// constellation, Nfft
int bps = 5, Nfft = 1024;
int Nbps = (int)pow(2, bps);
double sigPow = 1.0;
double snr, noise_mag;
double norms = 0;

int ch = 0;        // 0: AWGN, 1: Rayleigh fading


int bi2de(int* arr) { // Convert binary to decimal
    int num = 0;
    int j = 0;
    for (int i = bps - 1; i >= 0; i--)
        num += arr[i] * (int)pow(2, j++);
    return num;
};
int* de2bi(int* arr2, int n, int bps) { // Convert decimal to binary
    int* arr = (int*)calloc(bps, sizeof(int));
    for (int i = 0; n != 0; i++) {
        arr[i] = n % 2;
        n /= 2;
    }
    for (int i = 0; i < bps; i++)
        arr2[bps - i - 1] = arr[i];
    free(arr);
    return arr2;
};
double gaussianRandom(double average, double stdev) { // Returns a Gaussian random value of a normal distribution with mean and standard deviation of stdev
    double v1, v2, s, temp;

    do {
        v1 = 2 * ((double)rand() / RAND_MAX) - 1;      // range : -1.0 ~ 1.0 (random value)
        v2 = 2 * ((double)rand() / RAND_MAX) - 1;      // range : -1.0 ~ 1.0 (random value)
        s = v1 * v1 + v2 * v2;
    } while (s >= 1 || s == 0);

    s = sqrt((-2 * log(s)) / s);

    temp = v1 * s;
    temp = (stdev * temp) + average;

    return temp;
}
void rand_gaussian_temp(double* arr, int n, double stdev) {
    for (int i = 0; i < n; i++) {
        arr[i] = gaussianRandom(0, stdev);
    }
}
double sum_all_symbol(complex* con, complex channel, int con_select, int n_sym, complex R_Sym, double noise_mag) {
    double prob = 0;
    if (con_select == 0) {     
        for (int i = 0; i < n_sym; i++) {   // con[0~15]
            prob += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(R_Sym - con[i] * channel, 2))); // sum to probability of received symbol and each symbols(no scaling 0~15)
        }
    }
    else if (con_select == 1) { 
        for (int i = 0; i < n_sym; i++) {   // com[16~31]
            prob += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(R_Sym - con[i + 16] * channel, 2))); // sum to probability of received symbol and each symbols(scaling 16~31)
        }
    }

    return prob;
}
void LLR_index(complex* SymRx, double* received_LR, complex* channel, complex* con, double noise_mag) {
    double sum_0 = 1, sum_1 = 1;

    for (int subc = 0; subc < Nfft; subc++)
    {
        sum_0 *= sum_all_symbol(con, channel[subc], 0, (int)pow(2, bps - 1), SymRx[subc], noise_mag); // probability at all subcarriers for con_select = 0
        sum_1 *= sum_all_symbol(con, channel[subc], 1, (int)pow(2, bps - 1), SymRx[subc], noise_mag); // probability at all subcarriers for con_select = 1
        received_LR[subc] = log(sum_0) - log(sum_1);                                                  // calculate LLR values for each subcarrier

        if (received_LR[subc] > 30.0) received_LR[subc] = 30.0;
        else if (received_LR[subc] < -30.0) received_LR[subc] = -30.0;
        else received_LR[subc] = received_LR[subc];

        sum_0 = 1;  
        sum_1 = 1;
    }
}
double sym_bit_calculator(complex SymRx, complex channel, complex* con, int n_th_bit, int con_selec, double noise_mag) {
    double sum_0 = 0, sum_1 = 0;
    for (int i = 0; i < Nbps / 2; i++) {
        int j;

        if (con_selec == 0) j = i;      // if the MSB estimated bit is 0, j = i
        else j = i + 16;                // if the MSB estimated bit is 1, j = i+16

        if (n_th_bit == 0) { // 2nd bit
            if (i <= 7)      // probability that the 2nd bit is 0
                sum_0 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
            else
                sum_1 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
        }
        else if (n_th_bit == 1) { // 3rd bit
            if (i / 4 == 0 || i / 4 == 2)
                sum_0 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
            else
                sum_1 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
        }
        else if (n_th_bit == 2) { // 4th bit
            if (i % 4 == 0 || i % 4 == 1)
                sum_0 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
            else
                sum_1 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
        }
        else if (n_th_bit == 3) { // 5th bit
            if (i % 4 == 0 || i % 4 == 2)
                sum_0 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
            else
                sum_1 += exp(-(1 / (2 * noise_mag * noise_mag)) * abs(pow(SymRx - con[j] * channel, 2)));
        }
    }
    double llr = log(sum_0) - log(sum_1);

    if (llr > 30.0) llr = 30.0;
    else if (llr < -30.0) llr = -30.0;
    return llr;
}
void LLR_symbol(complex* SymRx, double* received_LR, char* est_1st_bit, complex* channel, complex* con, double noise_mag) {
    for (int subc = 0; subc < Nfft; subc++) {
        for (int bit = 0; bit < bps - 1; bit++) {
            received_LR[subc * (bps - 1) + bit] = sym_bit_calculator(SymRx[subc], channel[subc], con, bit, est_1st_bit[subc], noise_mag);
        }
    }
}

void OFDM_1(
    double* index_LR, complex* received_symbol, char* cw_1st, char* cw_2nd, double EsN0, double con_k, complex* channel
) {
    int* arr = (int*)malloc(sizeof(int) * bps);                  // to make a binary symbol into decimal
    int* bit_Tx = (int*)malloc(sizeof(int) * Nfft);              // transmitted symbols (decimal)
    complex* SymTx = (complex*)malloc(sizeof(complex) * Nfft);   // transmitted symbols (complex)
    complex* SymTx2 = (complex*)malloc(sizeof(complex) * Nfft);  
    complex* noise = (complex*)malloc(sizeof(complex) * Nfft);   // AWGN noise
    complex* SymRx = (complex*)malloc(sizeof(complex) * Nfft);   // received symbols (complex)
    complex* chan = (complex*)malloc(sizeof(complex) * Nfft);
    double* real = (double*)calloc(Nfft, sizeof(double));        // real of the noise
    double* imag = (double*)calloc(Nfft, sizeof(double));        // imag of the noise

    // constellation
    complex* con = (complex*)malloc(sizeof(complex) * Nbps);

    // 32-DaSQAM
    double* inphase = new double[Nbps] { -3.0, -3.0, -3.0, -3.0, -1.0, -1.0, -1.0, -1.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, -3.0, -3.0, -3.0, -3.0, -1.0, -1.0, -1.0, -1.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0 };
    double* quadrature = new double[Nbps] { 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0};

    for (int i = 0; i < Nbps; i++)
    {
        if (i < Nbps / 2)
            con[i] = complex(inphase[i], quadrature[i]);
        else
            con[i] = complex(inphase[i] * con_k, quadrature[i] * con_k);
    }

    delete[] inphase;
    delete[] quadrature;

    // normalization
    for (int i = 0; i < Nbps; i++)
        norms += pow(abs(con[i]), 2);

    norms = sqrt(norms / (Nbps));

    for (int i = 0; i < Nbps; i++)
        con[i] /= norms;

    for (int subc = 0; subc < Nfft; subc++) {
        for (int bit = 0; bit < bps; bit++)
        {
            if (bit == 0)
                arr[bit] = (int)cw_1st[subc];
            else
                arr[bit] = (int)cw_2nd[subc * (bps - 1) + (bit - 1)];
            if (bit == bps - 1)
                bit_Tx[subc] = bi2de(arr);
        }
    }

    // mapper
    for (int subc = 0; subc < Nfft; subc++)
        SymTx[subc] = con[bit_Tx[subc]];

    rand_gaussian_temp(real, Nfft, 1 / sqrt(2));
    rand_gaussian_temp(imag, Nfft, 1 / sqrt(2));

    // channel
    for (int i = 0; i < Nfft; i++)
    {
        if (ch == 0) 
        {
            chan[i] = complex(1.0, 0.0);
            channel[i] = chan[i];
        }
        else if (ch == 1)
        {
            chan[i] = complex(real[i], imag[i]);
            SymTx2[i] = SymTx[i] * chan[i];
            channel[i] = chan[i];
        }
    }

    sigPow = 1.0;

    // add AWGN noise
    double noise_mag = sqrt(pow(10, (-EsN0 / 10)) * sigPow / 2);

    // noise
    rand_gaussian_temp(real, Nfft, noise_mag);
    rand_gaussian_temp(imag, Nfft, noise_mag);

    for (int i = 0; i < Nfft; i++)
        noise[i] = complex(real[i], imag[i]);

    for (int i = 0; i < Nfft; i++) {
        if (ch == 0) {
            SymRx[i] = SymTx[i] + noise[i];
            received_symbol[i] = SymRx[i];
        }
        else if (ch == 1) {
            SymRx[i] = SymTx2[i] + noise[i];
            received_symbol[i] = SymRx[i];
        }
    }
        
    LLR_index(SymRx, index_LR, chan, con, noise_mag);

    free(bit_Tx);
    free(SymTx);
    free(SymRx);
    free(noise);
    free(imag);
    free(real);
    free(chan);
    free(arr);
    free(con);
}

void OFDM_2(
    double* symbol_LR, complex* received_symbol, char* est_1st_bit, double EsN0, double con_k, complex* channel
) {
    double* real = (double*)calloc(Nfft, sizeof(double));      
    double* imag = (double*)calloc(Nfft, sizeof(double));     

    // constellation
    complex* con = (complex*)malloc(sizeof(complex) * Nbps);

    // 32-DaSQAM
    double* inphase = new double[Nbps] { -3.0, -3.0, -3.0, -3.0, -1.0, -1.0, -1.0, -1.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, -3.0, -3.0, -3.0, -3.0, -1.0, -1.0, -1.0, -1.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0 };
    double* quadrature = new double[Nbps] { 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0, 3.0, 1.0, -3.0, -1.0};

    for (int i = 0; i < Nbps; i++)
    {
        if (i < Nbps / 2)
            con[i] = complex(inphase[i], quadrature[i]);
        else
            con[i] = complex(inphase[i] * con_k, quadrature[i] * con_k);
    }

    delete[] inphase;
    delete[] quadrature;

    // normalization
    for (int i = 0; i < Nbps; i++)
        norms += pow(abs(con[i]), 2);

    norms = sqrt(norms / (Nbps));

    for (int i = 0; i < Nbps; i++)
        con[i] /= norms;

    sigPow = 1.0;

    // add AWGN noise
    double noise_mag = sqrt(pow(10, (-EsN0 / 10)) * sigPow / 2);

    LLR_symbol(received_symbol, symbol_LR, est_1st_bit, channel, con, noise_mag);

    free(received_symbol);
    free(channel);
    free(con);
}
