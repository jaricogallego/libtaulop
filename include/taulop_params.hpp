//
//  taulop_params.hpp
//  TauLopCost
//
//  Created by jarico on 15/5/16.
//  Copyright © 2016 Juan A. Rico. All rights reserved.
//

#ifndef taulop_params_hpp
#define taulop_params_hpp

#include "taulop_config.h"
#include "taulop_params_channel.hpp"

#include <iostream>
#include <vector>
using namespace std;


// TBD: Comment the attributes and methods.


// Singleton pattern
class TauLopParam {
    
private:
        
    vector<TaulopParamChannel *> channel;
    
    // P2P cost in both channels. (TBD: maybe it must be a class for any alg.)
    vector<double **> p2p;
    
    vector<long> sizes;
    
    int  max_idx;
    int  max_tau;
    
    static bool instanceFlag;
    static TauLopParam *single;
    
    static vector<string> channel_names;

    
    TauLopParam(); //private constructor
    
    void  setP2P ();
    
    
public:
    
    static void setInstance(vector<string> channel_names);
    static TauLopParam* getInstance();
    
    ~TauLopParam();
    
    double  getTime  (long   m, int tau, int chn);
    long    getBytes (double t, int tau, int chn);
    
    void show ();
};



#endif /* taulop_params_hpp */
