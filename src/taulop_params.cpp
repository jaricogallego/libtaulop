//
//  taulop_params.cpp
//  TauLopCost
//
//  Created by jarico on 15/5/16.
//  Copyright © 2016 Juan A. Rico. All rights reserved.
//

#include "taulop_params.hpp"
#include "config.h"

#include <cmath>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
using namespace std;


// Initialization of tatic variables
TauLopParam* TauLopParam::single = nullptr;
vector<string> TauLopParam::networks;



// PRIVATE methods

//private constructor
TauLopParam::TauLopParam() {

    // Create communication channels information
    int num = 0;
    for (auto it = TauLopParam::networks.begin(); it != TauLopParam::networks.end(); ++it, num++) {
        this->channel.push_back(new TaulopParamChannel (*it, num));
    }
    
    // Take data from the first channel (we assume that all channels have the next same values)
    this->max_tau = this->channel[0]->getNumTau(); // TBD
    this->max_idx = this->channel[0]->getNumM(); // TBD
    this->sizes   = this->channel[0]->getSizes();
    
    // Create the matrices to store the data of P2P transmissions
    for (int channel_nr = 0; channel_nr < this->channel.size(); channel_nr++) {
        double **chn = new double * [this->channel[channel_nr]->getNumTau()];
        for (int i = 0; i < this->max_tau; i++) {
            chn[i] = new double [this->max_idx];
        }
        this->p2p.push_back(chn);
    }
        
    // Set P2P transmission times values
    this->setP2P();
    
#if TLOP_DEBUG == 1
    show();
#endif
}



void TauLopParam::setP2P  () {
    
    int chn_nr = 0;
    
    for (auto it = TauLopParam::networks.begin(); it != TauLopParam::networks.end(); ++it, chn_nr++) {
        
        double **P2P = this->p2p[chn_nr];

        for (int tau = 0; tau < this->max_tau; tau++) {
            
            for (int idx = 0; idx < this->max_idx; idx++) {
                
                // Cost depends on the channel type (name)
                // TBD: this cost formulation should be provided by the user.
                
                double T = 0.0;
                string chn_type = *it;
                
                if (chn_type == SHM_NET) {
                    
                    double o0 = this->channel[chn_nr]->getO(idx);
                    double L0 = this->channel[chn_nr]->getL(idx, tau);
                    
                    T = o0 + (2 * L0);
                    
                } else if (chn_type == IB_NET) {
                    
                    double o1 = this->channel[chn_nr]->getO(idx);
                    double L1 = this->channel[chn_nr]->getL(idx, tau);

                    T = o1 + L1;
                    //T = L1;  // RICO: PROB. MENSAJES PEQUEÑOS
                    
                } else if (chn_type == TCP_NET) {
                    
                    double o1 = this->channel[chn_nr]->getO(idx);
                    double L0 = this->channel[0]->getL(idx, tau);
                    double L1 = this->channel[chn_nr]->getL(idx, tau);
                    
                    T = o1 + (2 * L0) + L1;  // RICO: PROB. MENSAJES PEQUEÑOS
                    //T = (2 * L0) + L1;
                    
                } else {
                    cerr << "ERROR: channel name not known: " << chn_type << endl;
                }                
                
                P2P[tau][idx] = T;
            }
        }

    }
}




// PUBLIC interface

TauLopParam::~TauLopParam () {
    
    TauLopParam::single = nullptr;
    
    for (int chn_nr; chn_nr < this->networks.size(); chn_nr++) {
        delete this->channel[chn_nr];
    }
    
    for (auto it = this->p2p.begin(); it != this->p2p.end(); ++it) {
        
        double **p2p;
        for (int i = 0; i < this->max_tau; i++) {
            delete [] p2p[i];
        }
        delete [] p2p;
    }
}


void TauLopParam::setInstance(vector<string> networks) {
    
    if (!TauLopParam::single) {
        
        TauLopParam::networks = networks;
        TauLopParam::single   = new TauLopParam();
        
    } else {
        
        cerr << "ERROR: network parameters already loaded from: " << endl;
        
        for (auto i = TauLopParam::networks.begin(); i != TauLopParam::networks.end(); ++i) {
            cout << *i << endl;
        }
        
    }
}


TauLopParam* TauLopParam::getInstance() {
    
    if (!TauLopParam::single) {
        TauLopParam::single = new TauLopParam();
    }
    
    return single;
}


double TauLopParam::getTime (long m, int tau, int chn) {
    
    double   t = 0.0;
    double **p2p;
    int      idx;
    

    if (m <= 0) {
        cerr << "ERROR: message size must be greater than 0: " << m << endl;
        return 0.0;
    }

    if (chn >= this->networks.size()) {
        cerr << "ERROR: unknown channel: " << chn << endl;
        return -1;
    }

    if (tau > this->max_tau) {
        cerr << "ERROR: value of tau is too high: " << tau << " (" << this->max_tau << ")" << endl;
        tau = max_tau;
        //return -1;
    }
    
    
    p2p = this->p2p[chn];
    
    for (idx = 0; idx < this->max_idx; idx++) {
        if (this->sizes[idx] >= m) break;
    }
    
    if (this->sizes[idx] == m) {
        t = p2p[tau-1][idx];
    }
    
    else if (idx == 0) {
        t = p2p[tau-1][0];
    }
    
    else if (idx == this->max_idx) {
        t = p2p[tau-1][idx-1] + getTime(m - this->sizes[this->max_idx-1], tau, chn);
    }
    
    else {
        
        int idx_up = idx;
        int idx_dw = idx - 1;
        
        double t_up = p2p[tau-1][idx_up];
        double t_dw = p2p[tau-1][idx_dw];
        
        long n_up = this->sizes[idx_up];
        long n_dw = this->sizes[idx_dw];
        
        t = t_dw + ((m - n_dw) * (t_up - t_dw)) / (n_up - n_dw);
        
        // Security rule: if T < 0, t_up < t_dw, then use t_dw.
        if (t < 0) {
            t = t_dw;
        }
    }
    
    return t;
}


long TauLopParam::getBytes (double t, int tau, int chn) {
    
    long     m = 0;
    int      idx;
    double **p2p;
    

    if (t <= 0.0) {
        cerr << "ERROR: time must be greater than 0.0: " << t << endl;
        return 0;
    }
    
    if (chn >= this->networks.size()) {
        cerr << "ERROR: unknown channel: " << chn << endl;
        return -1;
    }
    
    if (tau > this->max_tau) {
        cerr << "ERROR: value of tau is too high: " << tau << " (" << this->max_tau << ")" << endl;
        tau = max_tau;
        //return -1;
    }

    
    p2p = this->p2p[chn];
    
    for (idx = 0; idx < this->max_idx; idx++) {
        if (p2p[tau-1][idx] > t) break;
    }
    
    if (idx == 0) {
        return 1; // Time short. Only 1 block.
    }
    
    if (idx == this->max_idx) {
        
        if (t == p2p[tau-1][this->max_idx-1]) {
            
            m = this->sizes[this->max_idx-1];
            
        } else {
        
            double t_max = p2p[tau-1][this->max_idx-1];
            m = this->sizes[this->max_idx-1] + getBytes(t - t_max, tau, chn);
            
        }
        
    } else { // Time is in an interval
        
        long b_up = this->sizes[idx];
        double t_up = p2p[tau-1][idx];
        
        long b_dw = this->sizes[idx-1];
        double t_dw = p2p[tau-1][idx-1];
                
        m = b_dw + ceil((t - t_dw) * double(b_up - b_dw) / (t_up - t_dw));
        
        if (m < 1) {
            m = 1;
        }
    }
    
    return m;
}


void TauLopParam::show () {
    
    int chn_nr = 0;
    
    for (auto it = this->networks.begin(); it != this->networks.end(); ++it, chn_nr++) {
        
        cout << "Network:  " << *it << endl;
        double **P2P = this->p2p[chn_nr];
        
        for (int i = 0; i < this->max_idx; i++) {
            
            cout << this->sizes[i] << ")  ";
            
            for (int j = 0; j < this->max_tau; j++) {
                cout << fixed << setprecision(9) << P2P[j][i] << " \t ";
            }
            cout << endl;
        }
        cout << endl;
        
    }
}




