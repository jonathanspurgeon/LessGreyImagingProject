#include "network.h"
#include "tools.h"

#include <boost/math/special_functions/erf.hpp>

using namespace std;

void network::pulsateVessel()
{
    initialiseSimulation();

    cout<<"Starting Pulsating... "<<endl;

    double startTime,endTime;
    startTime=tools::getCPUTime();

    double timeSoFar=0;

    //post-processing
    if(videoRecording)
        record=true;

    timeStep=0.001;
    double xCentre=xEdgeLength/2.;
    double yCentre=yEdgeLength/2.;
    double zCentre=zEdgeLength/2.;

    map<node*,double> nodalPositionX;
    map<node*,double> nodalPositionY;
    map<node*,double> nodalPositionZ;
    map<node*,double> nodalRadius;
    map<pore*,double> poreRadius;
    for (int i = 0; i < totalNodes; ++i) {
        node* n=getNode(i);
        if(!n->getClosed())
        {
            nodalPositionX[n]=n->getXCoordinate();
            nodalPositionY[n]=n->getYCoordinate();
            nodalPositionZ[n]=n->getZCoordinate();
            nodalRadius[n]=n->getRadius();
        }
    }

    for (int i = 0; i < totalPores; ++i) {
        pore* n=getPore(i);
        if(!n->getClosed())
        {
            poreRadius[n]=n->getRadius();
        }
    }

    while(timeSoFar<simulationTime)
    {
        timeSoFar+=timeStep;

        double a=0.2;
        double d=0.7;
        double h=3;
        double s=0.05;
        double w=0.02;
        double L=1.5*d;
        double t=timeSoFar-ceil(timeSoFar/L-0.5)*L;

        for (int i = 0; i < totalNodes; ++i) {
            node* n=getNode(i);
            if(!n->getClosed())
            {
                double pulse=(a*(exp(-pow(t+d,2)/(2*w))+exp(-pow(t-d,2)/(2*w)))+(h-abs(t/s)-t)*exp(-pow(7*t,2)/2))*0.04+1;
                double scaleFactor=1-sqrt(pow(nodalPositionX[n]-xCentre,2)+pow(nodalPositionY[n]-yCentre,2)+pow(nodalPositionZ[n]-zCentre,2))/xEdgeLength;
                n->setXCoordinate(xCentre+(nodalPositionX[n]-xCentre)*pulse*scaleFactor);
                n->setYCoordinate(yCentre+(nodalPositionY[n]-yCentre)*pulse*scaleFactor);
                n->setZCoordinate(zCentre+(nodalPositionZ[n]-zCentre)*pulse*scaleFactor);
                double scaleR=(a*(exp(-pow(t+d,2)/(2*w))+exp(-pow(t-d,2)/(2*w)))+(h-abs(t/s)-t)*exp(-pow(7*t,2)/2))*0.1+1;
                n->setRadius(nodalRadius[n]/scaleR);
            }
        }

        for (int i = 0; i < totalPores; ++i) {
            pore* n=getPore(i);
            if(!n->getClosed())
            {
                double scaleR=(a*(exp(-pow(t+d,2)/(2*w))+exp(-pow(t-d,2)/(2*w)))+(h-abs(t/s)-t)*exp(-pow(7*t,2)/2))*0.1+1;
                n->setRadius(poreRadius[n]/scaleR);
            }
        }

        emitPlotSignal();
        //Thread Management
        if(cancel)break;
    }

    //post-processing
    if(videoRecording)
    {
        record=false;
        extractVideo();
    }

    cout<<"Simulation Time "<<timeSoFar<<endl;

    endTime=tools::getCPUTime();
    cout<<"Processing Time: "<<endTime-startTime<<" s"<<endl;
}

void network::runCoupledCell()
{
    initialiseSimulation();

    cout<<"Starting Coupled Cell... "<<endl;

    double startTime,endTime;
    startTime=tools::getCPUTime();

    double timeSoFar=0;
    double timeToRemodel=0;

    //post-processing
    if(videoRecording)
        record=true;

    loop();
    end();

    cancel=false;

    setupTissueProperties();
    setupTAFDistribution();
    setupFNDistribution();
    assignInitialBloodViscosities();

    set<int> sproutTips;
    initialiseSroutTips(sproutTips);

    calculateTimeStepForAngio();

    while(timeSoFar<simulationTime)
    {
        timeSoFar+=timeStep;
        timeToRemodel+=timeStep;

        updateSproutTipPositions(sproutTips);

        setBranching(sproutTips);
        if(branchingWSS)
            setBranchingWSS(sproutTips);

        emitPlotSignal();

        //Thread Management
        if(cancel)break;
    }
    cancel=false;

    loop();
    end();

    //post-processing
    if(videoRecording)
    {
        record=false;
        extractVideo();
    }

    cout<<"Simulation Time "<<timeSoFar<<endl;

    endTime=tools::getCPUTime();
    cout<<"Processing Time: "<<endTime-startTime<<" s"<<endl;
}

void network::runStaticModelAdaptation()
{
    initialiseSimulation();

    cout<<"Starting Angiogenesis On-Lattice... "<<endl;

    double startTime,endTime;
    startTime=tools::getCPUTime();

    double timeSoFar=0;

    timeStep=0.005;

    //post-processing
    if(videoRecording)
        record=true;

    for(int i=0;i<totalPores;++i)
    {
        pore*p=getPore(i);
        p->setRadius(6e-6);
    }

    while(timeSoFar<simulationTime)
    {
        double remodel1=true;

        while(remodel1)
        {
            double ff=0;

            remodel1=phaseSeparation?solvePressureInAngioModelWthPhaseSeparation():
                                     solvePressureInAngioModel();
            while(ff<0.1)
            {
                double additionalTime=phaseSeparation?runHaematocritFlowWithPhaseSeparation():
                                                        runHaematocritFlow();
                ff+=additionalTime;

                //Thread management
                if(cancel)
                    break;
            }

            timeSoFar+=ff;

            if(shuntPrevention)
            {
                calculateConvectedStimuli();
                calculateConductedStimuli();
            }

            double remodel2=true;
            while(remodel2)
            {
                remodel2=recalculateRadii();

                //Thread management
                if(cancel)
                    break;
            }

            remodel1=phaseSeparation?solvePressureInAngioModelWthPhaseSeparation():
                                     solvePressureInAngioModel();

            //Thread Management
            if(cancel)break;
        }

        //Thread Management
        if(cancel)break;
    }

    //post-processing
    if(videoRecording)
    {
        record=false;
        extractVideo();
    }

    cout<<"Simulation Time "<<timeSoFar<<endl;

    endTime=tools::getCPUTime();
    cout<<"Processing Time: "<<endTime-startTime<<" s"<<endl;
}

void network::runTwoLayerDiffusion()
{
    cout<<"Starting Diffusion in Two Layers... "<<endl;

    double startTime,endTime;
    startTime=tools::getCPUTime();

    ofstream file("Results/releaseProfile.txt");
    ofstream file2("Results/AnalyticalResleaseProfile.txt");
    file<<"t EG_Sim DM_Sim"<<endl;
    file2<<"t EG_Analytical DM_Analytical"<<endl;

    double timeSoFar=0;
    cancel=false;

    double hx=xEdgeLength/meshSizeX;
    double coefX=1/pow(hx,2);

    double D1=DT;
    double D2=PVT;

    double initialDrugMass(0);
    double initialOxygenMass(0);

    for(int i=0;i<totalBlocks;++i)
    {
        block* b=getBlock(i);
        if(!b->getClosed())
        {
            if(b->getXCoordinate()<48e-6)
            {
                b->setConcentration(1);
                b->setOxygenConcentration(0);
            }
            else if (b->getXCoordinate()<90e-6)
            {
                b->setConcentration(0);
                b->setOxygenConcentration(1);
            }
            else
            {
                b->setConcentration(0);
                b->setOxygenConcentration(0);
            }
            initialDrugMass+=b->getConcentration()*b->getVolume();
            initialOxygenMass+=b->getOxygenConcentration()*b->getVolume();
        }
    }

    emitPlotSignal();

    //post-processing
    if(videoRecording)
        record=true;

    timeStep=1e50;
    for(int i=0;i<totalBlocks;++i)
    {
        block* n=getBlock(i);
        if(!n->getClosed())
        {

            double step=1./(2*max(DT,PVT)*coefX);
            if(step<timeStep)
                timeStep=step;
        }
    }

    cout<<"Time step with diffusion: "<<timeStep<<endl;

    double deltaT=timeStep/2;

    while(timeSoFar<simulationTime)
    {

        map<block*,double> blockConcentration;
        map<block*,double> blockConcentration2;

        for(int i=0;i<totalBlocks;++i)
        {
            block* n=getBlock(i);
            if(!n->getClosed() && n->getXCoordinate()<90e-6)
            {
                double ii=n->getX();
                double jj=n->getY();
                double kk=n->getZ();

                block *nW,*nE;
                nW=getBlock(ii-1,jj,kk);
                nE=getBlock(ii+1,jj,kk);

                double newConcentration(0);
                double newConcentration2(0);

                /////Drug

                newConcentration+=n->getConcentration()*(1-deltaT*(2*D1*coefX));
                if(nW!=0) if(!nW->getClosed()) newConcentration+=deltaT*D1*coefX*nW->getConcentration();
                if(nE!=0 && nE->getXCoordinate()<90e-6) if(!nE->getClosed()) newConcentration+=deltaT*D1*coefX*nE->getConcentration();

                //boundry
                if(nW==0) newConcentration+=deltaT*D1*coefX*n->getConcentration();

                /////Oxygen

                newConcentration2+=n->getOxygenConcentration()*(1-deltaT*(2*D2*coefX));
                if(nW!=0) if(!nW->getClosed()) newConcentration2+=deltaT*D2*coefX*nW->getOxygenConcentration();
                if(nE!=0 && nE->getXCoordinate()<90e-6) if(!nE->getClosed()) newConcentration2+=deltaT*D2*coefX*nE->getOxygenConcentration();

                //boundry conditions
                if(nW==0) newConcentration2+=deltaT*D2*coefX*n->getOxygenConcentration();

                blockConcentration[n]=newConcentration;
                blockConcentration2[n]=newConcentration2;
            }
        }

        double drugMass(0);
        double oxygenMass(0);

        for(map<block*,double>::iterator iterator = blockConcentration.begin(); iterator != blockConcentration.end(); iterator++)
        {
            block* b=iterator->first;
            double concentration=iterator->second;

            b->setConcentration(concentration);
            drugMass+=b->getConcentration()*b->getVolume();

            if(concentration<0 || concentration>1.01) {cout<<"block concentration1 out of range: "<<b->getConcentration()<<endl;cancel=true;}
        }

        for(map<block*,double>::iterator iterator = blockConcentration2.begin(); iterator != blockConcentration2.end(); iterator++)
        {
            block* b=iterator->first;
            double concentration=iterator->second;

            b->setOxygenConcentration(concentration);
            oxygenMass+=b->getOxygenConcentration()*b->getVolume();

            if(concentration<0 || concentration>1.01) {cout<<"block concentration2 out of range: "<<b->getOxygenConcentration()<<endl;cancel=true;}
        }

        timeSoFar+=deltaT;

        file<<timeSoFar/(24*7*3600)<<" "<<(1-drugMass/initialDrugMass)*20.46<<" "<<(1-oxygenMass/initialOxygenMass)*68.2<<endl;

        double drug1Release(0);
        double drug2Release(0);
        //Analytical solution
        for(int i=0;i<100;++i)
        {
            drug1Release+=8*90e-6*pow(-1,i)/(48e-6*pow(tools::pi(),2)*pow(2*i+1,2))*sin(48e-6*tools::pi()*(2*i+1)/(2*90e-6))*exp(-DT*pow(tools::pi(),2)*pow(2*i+1,2)*timeSoFar/(4*pow(90e-6,2)));
            drug2Release+=8*90e-6/(42e-6*pow(tools::pi(),2)*pow(2*i+1,2))*(1-pow(-1,i)*sin(48e-6*tools::pi()*(2*i+1)/(2*90e-6)))*exp(-PVT*pow(tools::pi(),2)*pow(2*i+1,2)*timeSoFar/(4*pow(90e-6,2)));
        }
        file2<<timeSoFar/(24*7*3600)<<" "<<(1-drug1Release)*20.46<<" "<<(1-drug2Release)*68.2<<endl;

        emitPlotSignal();

        //Thread Management
        if(cancel)break;
    }

    //post-processing
    if(videoRecording)
    {
        record=false;
        extractVideo();
    }

    cout<<"Simulation Time "<<timeSoFar<<endl;

    endTime=tools::getCPUTime();
    cout<<"Processing Time: "<<endTime-startTime<<" s"<<endl;
}

void network::runLayerDissolution()
{
    cout<<"Starting Dissolution In two layers... "<<endl;

    double startTime,endTime;
    startTime=tools::getCPUTime();

    double timeSoFar=0;
    cancel=false;

    double theta=-0.41;
    double s=0;

    double hx=xEdgeLength/meshSizeX;
    double coefX=1/pow(hx,2);

    for(int i=0;i<totalBlocks;++i)
    {
        block* b=getBlock(i);
        if(!b->getClosed())
        {
            if(b->getXCoordinate()<45e-6)
            {
                b->setConcentration(1);
                b->setConnectedToVessel(true);
            }

            else
            {
                b->setConcentration(0);
                b->setConnectedToVessel(false);
            }

        }
    }

    emitPlotSignal();

    //post-processing
    if(videoRecording)
        record=true;

    timeStep=1e50;
    for(int i=0;i<totalBlocks;++i)
    {
        block* n=getBlock(i);
        if(!n->getClosed())
        {

            double step=1./(2*(DT)*coefX);
            if(step<timeStep)
                timeStep=step;
        }
    }

    double deltaT=timeStep/2;


    while(timeSoFar<simulationTime)
    {

        timeSoFar+=deltaT;
        s=theta*sqrt(timeSoFar);
        if(s>=-1)

        {
            for(int i=0;i<totalBlocks;++i)
            {
                block* n=getBlock(i);
                if(!n->getClosed() && n->getXCoordinate()<(1-abs(s))*45e-6)
                {
                    n->setConnectedToVessel(true);
                }
                else if(!n->getClosed() && n->getXCoordinate()>=(1-abs(s))*45e-6)
                {
                    double x=n->getXCoordinate()/45e-6-1;
                    n->setConnectedToVessel(false);
                    n->setConcentration((boost::math::erf(x/(2*sqrt(timeSoFar)))-1)/(boost::math::erf(theta/2)-1));
                }
            }
        }
        else
        {
            map<block*,double> blockConcentration;

            for(int i=0;i<totalBlocks;++i)
            {
                block* n=getBlock(i);
                if(!n->getClosed())
                {
                    double ii=n->getX();
                    double jj=n->getY();
                    double kk=n->getZ();

                    block *nW,*nE;
                    nW=getBlock(ii-1,jj,kk);
                    nE=getBlock(ii+1,jj,kk);

                    double newConcentration(0);

                    newConcentration+=n->getConcentration()*(1-deltaT*(2*DT*coefX));
                    if(nW!=0) if(!nW->getClosed()) newConcentration+=deltaT*DT*coefX*nW->getConcentration();
                    if(nE!=0) if(!nE->getClosed()) newConcentration+=deltaT*DT*coefX*nE->getConcentration();
                    //boundry
                    if(nW==0) newConcentration+=deltaT*DT*coefX*n->getConcentration();
                    if(nE==0) newConcentration+=deltaT*DT*coefX*n->getConcentration();

                    blockConcentration[n]=newConcentration;
                }
            }

            for(map<block*,double>::iterator iterator = blockConcentration.begin(); iterator != blockConcentration.end(); iterator++)
            {
                block* b=iterator->first;
                double concentration=iterator->second;

                b->setConcentration(concentration);
                b->setConnectedToVessel(false);

                if(concentration<0 || concentration>1.01) {cout<<"block concentration1 out of range: "<<b->getConcentration()<<endl;cancel=true;}
            }
        }

        emitPlotSignal();

        //Thread Management
        if(cancel)break;
    }

    //post-processing
    if(videoRecording)
    {
        record=false;
        extractVideo();
    }

    cout<<"Simulation Time "<<timeSoFar<<endl;

    endTime=tools::getCPUTime();
    cout<<"Processing Time: "<<endTime-startTime<<" s"<<endl;
}


//    ofstream file2("heatmap.txt");

//    double unitLength=20*length/200;

//    vector<vector<vector<bool> > > isVessel;

//    isVessel.resize(200);
//    for (int i = 0; i < 200; ++i) {
//    isVessel[i].resize(200);

//    for (int j = 0; j < 200; ++j)
//      isVessel[i][j].resize(200);
//    }


//    for (int i = 0; i < 200; ++i)
//        for (int j = 0; j < 200; ++j)
//            for (int k = 0; k < 200; ++k)
//            {isVessel[i][j][k]=false;}

//    for (int i = 0; i < 20; ++i)
//        for (int j = 0; j < 20; ++j)
//            for (int k = 0; k < 20; ++k)
//            {
//                pore* pX=getPoreXout(i,j,k);
//                if(pX->getPhaseFlag()=='o')
//                {
//                    int thickness=floor(pX->getRadius()/unitLength);
//                    for (int ii = 10*i; ii < 10*i+10; ++ii)
//                        for (int jj = 10*j-thickness; jj <= 10*j+thickness; ++jj)
//                            for (int kk = 10*k-thickness; kk <= 10*k+thickness; ++kk)
//                            {
//                                if(jj>=0 && kk>=0 && jj<200 && kk<200)
//                                    isVessel[ii][jj][kk]=true;
//                            }
//                }

//                pore* pY=getPoreYout(i,j,k);
//                if(pY->getPhaseFlag()=='o')
//                {
//                    int thickness=floor(pY->getRadius()/unitLength);
//                    for (int ii = 10*i-thickness; ii <= 10*i+thickness; ++ii)
//                        for (int jj = 10*j; jj < 10*j+10; ++jj)
//                            for (int kk = 10*k-thickness; kk <= 10*k+thickness; ++kk)
//                            {
//                                if(ii>=0 && kk>=0 && ii<200 && kk<200)
//                                    isVessel[ii][jj][kk]=true;
//                            }
//                }

//                pore* pZ=getPoreZout(i,j,k);
//                if(pZ->getPhaseFlag()=='o')
//                {
//                    int thickness=floor(pZ->getRadius()/unitLength);
//                    cout<<thickness<<" ";
//                    for (int ii = 10*i-thickness; ii <= 10*i+thickness; ++ii)
//                        for (int jj = 10*j-thickness; jj <= 10*j+thickness; ++jj)
//                            for (int kk = 10*k; kk < 10*k+10; ++kk)
//                            {
//                                if(jj>=0 && ii>=0 && jj<200 && ii<200)
//                                    isVessel[ii][jj][kk]=true;
//                            }
//                }
//            }

//    for (int i = 0; i < 200; ++i)
//        for (int j = 0; j < 200; ++j)
//            for (int k = 0; k < 200; ++k)
//            {file2<<i*unitLength<<" "<<j*unitLength<<" "<<k*unitLength<<" "<<isVessel[i][j][k]<<endl;}




//        ifstream file("heatmap.txt");
//        ofstream file2("heatmapReduced.txt");

//        double value1,value2,value3;
//        int value;
//        while(file>>value1)
//        {
//            file>>value2;
//            file>>value3;
//            file>>value;
//            if(value==1)
//                file2<<value1<<" "<<value2<<" "<<value3<<" "<<value<<endl;
//        }
