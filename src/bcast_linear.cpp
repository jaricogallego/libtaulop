//
//  bcast_binomial.cpp
//  TauLopCost
//
//  Created by jarico on 17/Nov/16.
//  Copyright © 2016 Juan A. Rico. All rights reserved.
//

#include "bcast_linear.hpp"

#include "transmission.hpp"
#include "collective.hpp"
#include "communicator.hpp"
#include "taulop_concurrent.hpp"
#include "taulop_sequence.hpp"

#include <cmath>
#include <iostream>
using namespace std;




BcastLinear::BcastLinear () {
    
}

BcastLinear::~BcastLinear () {
    
    
}


double BcastLinear::evaluate (Communicator *comm, int *size, int root) {
    
    TauLopConcurrent *conc;
    TauLopSequence   *seq;
    Transmission     *c;
    Process          *p_src, *p_dst;
    
    int P = comm->getSize();
    
    double tm = 0.0;
    
    conc = new TauLopConcurrent ();
    seq  = new TauLopSequence ();
    
    for (int rank = 0; rank < P; rank++) {
        
        if (rank == root)  continue;

        p_src = new Process (root, comm->getNode(root));
        p_dst = new Process (rank, comm->getNode(rank));
        
        int tau = 1;
        
        c = new Transmission(p_src, p_dst, *size, tau);
        seq->add(c);
    }
    
    conc->add(seq);
    
    TauLopCost *t = new TauLopCost();
    
#if TLOP_DEBUG == 1
    cout << " ----  Root " << root << endl;
    conc->show();
#endif
    
    conc->evaluate(t);
    
#if TLOP_DEBUG == 1
    cout << "  --------  Cost:  " << t->getTime() << endl;
    t->show();
#endif
    
    tm = t->getTime();
    
    delete t;
    delete conc;
    
    
    return tm;
}
