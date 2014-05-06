#include <fstream>
#include <ctime>
#include <assert.h>
#include "common.h"
#include "GCrod.h"
#include "Fanglespring2.h"

int main()
{
    //! [Parameter]
    const double temperature = 298;
    const double viscosity = 0.0008904;
    const double timeStep = 5e-12;
    const unsigned long steps = 55000000;
    const double rodLength = 4e-09;
    const double beadRadius = 1.6e-09;
    const double persistence = 50e-9;
    const unsigned beadNumber = 11; 
    const double stiffness = M_KB*temperature*persistence/rodLength;
    const unsigned long skipSteps = 0;
    const unsigned interval = 100;
    const double sf = 20;
    const double stretch = M_KB*temperature*sf/persistence;
    //! [Parameter end]

    //! [wlc chain]
    VectorXV3d wlc(beadNumber);
    //initialize ...
    wlc[0] = Vector3d(0,0,0);
    for (unsigned i=1; i < beadNumber; ++i)
        wlc[i] = wlc[i-1] + Vector3d(0,0,rodLength);
    //! [wlc chain end]
    //! [Output]
    unsigned long count = 0;
    
    std::ofstream fscale("scale.st",std::ios::app);
    long double scale = 0.0;
    //! [Output end]
    LINCS lincs(temperature,viscosity,timeStep);
    lincs.setSeed(1);
    for(unsigned long istep=0;istep<steps;++istep)
    {
        //! [wlc chain lincs]
        lincs.swapBeadPositions(wlc);
        lincs.resizeGC(beadNumber-1);
        lincs.setBeadRadius(beadRadius);
        lincs.reset();
        //forces,constraints,...
        //rod constraints.
        for (unsigned i=0; i<beadNumber-1; ++i)
            rodLink(rodLength,i,i+1,lincs);        
        //angle force
        for (unsigned i=1; i<beadNumber-1; ++i)       
            Fanglespring2_pi(stiffness,i-1,i,i+1,lincs);
        
        //stretching force
        lincs.plusF(0           , -stretch*Vector3d::UnitZ());
        lincs.plusF(beadNumber-1,  stretch*Vector3d::UnitZ()); 
        
        lincs();
        lincs.swapBeadPositions(wlc);
        //! [wlc chain lincs end]
        //! [Output lincs]
//        if(istep > skipSteps && istep % interval == 0)
//        {
           ++count;
           
           scale += wlc[beadNumber-1].z()-wlc[0].z();
//            if(istep%(interval*100) == 0)
//            {
            std::cout << wlc[beadNumber-1][2]-wlc[0][2] << std::endl;
           std::cout << scale/count << std::endl;
//            }
//        }
        //! [Output lincs end]
    }
    //! [Output out]
    
    fscale << scale/count << std::endl;
    //! [Output out end]
}

