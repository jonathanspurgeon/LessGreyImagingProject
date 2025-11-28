#ifndef NETWORK_H
#define NETWORK_H

#include "node.h"
#include "pore.h"
#include "cluster.h"
#include "block.h"
#include "particle.h"
#include "PDE.h"
#include "Cell.h"
#include "Fibre.h"

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <utility>
#include <fstream>
#include <cmath>

#include <boost/random/mersenne_twister.hpp>

#include <QObject>

using namespace std;

class network : public QObject
{
    Q_OBJECT
public:
    explicit network(QObject *parent = 0);
    ~network();
    void destroy();
    void reset();
    void setupModel();

    ////Regular Model
    void setupRegularModel();
    void createNodes();
    void createPores();
    void setNeighboors(); 
    void setNeighboorsForEntireNetwork();
    void applyCoordinationNumber();
    void assignRadii();
    void assignLengths();
    void distortNetwork();
    void assignShapeFactors();
    void assignShapeFactorConstants();
    void assignVolumes();
    void assignConductivities();
    void setActiveElements();
    void assignViscosities();
    void updateRanking();
    void setNeighboorsForGenericModel();
    void cleanGenericNetwork();

    ////Tissue
    void generateTissue();
    void tissueNetworkCollisionAnalysis();
    void tissueNetworkCollisionAnalysisRegular();
    void setupTissueProperties();

    ////Simulations
    void runSimulation();
    void runDrugFlowWithoutDiffusion();
    void runDrugFlowWithDiffusion();
    void runParticleFlow();
    void runAngiogenesisOnLattice();
    void runStaticModelAdaptation();
    void runRetinaModel();
    void runCellProliferation();

    ////Mouse Network
    void setupMouseNetwork();
    void loadMouseNetworkFromFiles();

    ////Artificial Network
    void setupArtificialNetwork();
    void buildArtificialNetwork();
    void buildArtificialNetwork2();
    void buildArtificialNetwork3();
    void buildArtificialNetwork4();
    void buildArtificialNetwork5();
    void buildArtificialNetwork6();
    void buildArtificialNetwork7();

    ////Angiogenesis Network
    void createParentVessel();
    void createParentVessel3D();
    void setupTAFDistribution();
    void setupFNDistribution();
    void initialiseSroutTips(std::set<int> &);
    void calculateTimeStepForAngio();
    void updateChemicalConcentrations();
    void updateSproutTipPositions(std::set<int> &);
    void setBranching(std::set<int> &);
    void setBranchingWSS(std::set<int> &);
    void assignInitialBloodViscosities();
    void assignBloodViscosities();
    void calculateConvectedStimuli();
    void calculateConductedStimuli();
    void computeConvectiveStimuliRecursive(node *n, double stimulus);
    void computeConductiveStimuliRecursive(node *n, double stimulus);
    bool solvePressureInAngioModel();
    bool solvePressureInAngioModelWthPhaseSeparation();
    double runHaematocritFlow();
    double runHaematocritFlowWithPhaseSeparation();
    void remodelVasculature();
    bool recalculateRadii(double time=0);
    double Chi_func(double);
    std::vector<double> generateEndothelialCellProbabilities(node *);
    node *addNode(node *, int, int, int);
    pore *addVessel(node* ,node*, int);

    ///Retina Model
    void createRetinaParentVessels();

    ///Pure Diffusion/Convection
    void runTwoLayerDiffusion();
    void runLayerDissolution();

    ///Pulsating
    void pulsateVessel();

    ///Coupled Cell
    void runCoupledCell();

    ////Solvers and Permeabilities
    void solvePressures();
    void solvePressuresForRegularModel();
    double updateFlows();
    void calculatePermeabilityAndPorosity();


    ////Misc

    //tools
    double getOutletFlow();

    //Video Recording
    void extractVideo();

    //randomness
    int uniform_int(int a=0, int b=1);
    double uniform_real(double a=0, double b=1);
    double rayleigh(double, double, double);
    double triangular(double, double, double);
    double normal(double,double,double,double);
    double weibull(double,double,double,double);

    //data extraction
    void extractDrugFlowResults(double, double, double &, int &, bool forceExtraction=false);
    void extractParticleFlowResults(double, double, double &, int &, bool forceExtraction=false);

    //initialisation
    void initialiseSimulation();

    //flow rate
    void setConstantFlowRateAker();
    void setConstantFlowRateSecant(std::set<pore*> &);
    void massConservationCheck();

    ////clustering
    int hkFind(int, std::vector<int>&);
    int hkUnion(std::vector<int>&,std::vector<int>&);
    int hkMakeSet(std::vector<int>&);
    //Regular
    void clusterPores(cluster*(pore::*)(void) const,void(pore::*)(cluster*),char(pore::*)(void) const,char,std::vector<cluster *> &);
    void clusterEverything();

    ////loading data
    void loadData();
    void loadNetworkData();
    void loadTwoPhaseData();

    ///// Plotting
    void emitPlotSignal();

    ////Access to pores/nodes/elements
    pore *getPoreX(int,int,int) const;
    pore *getPoreY(int,int,int) const;
    pore *getPoreZ(int,int,int) const;
    pore *getPoreXout(int,int,int) const;
    pore *getPoreYout(int,int,int) const;
    pore *getPoreZout(int,int,int) const;
    pore *getPore(int) const;
    node *getNode(int,int,int) const;
    node *getNode(int) const;
    particle *getParticle(int) const;
    block *getBlock(int,int,int) const;
    block *getBlock(int) const;

    int getTotalPores() const;
    int getTotalNodes() const;
    int getTotalBlocks() const;
    int getTotalOpenedPores() const;

    ////Getters/Setters

    //ThreadManagement
    bool getReady() const;
    void setCancel(bool value);
    bool getSimulationRunning() const;
    void setSimulationRunning(bool value);

    bool getRecord() const;
    bool getVideoRecording() const;
    bool getStartRecording() const;
    void setStartRecording(bool value);

    //Getters for network attributes

    int getNetworkSource() const;
    int getTotalOpenedNodes() const;
    int getTotalParticles() const;
    double getAbsolutePermeability() const;
    double getPorosity() const;

    double getXEdgeLength() const;
    double getYEdgeLength() const;
    double getZEdgeLength() const;

    int getNx() const;
    void setNx(int value);

    int getNy() const;
    void setNy(int value);

    int getNz() const;
    void setNz(int value);


    /////Cell Proliferation

    struct Box {
      /// @brief cells belonging to box
      std::vector<Cell> cells;
      /// @brief (pointer to) fibres belonging to box
      std::vector<Fibre*> p_fibres;
      // @ brief mesh elements contained in the box
      std::vector<int>  v_triangles;
    };

    /// @brief input parameters
    Param params;

    std::string input_file_name;

    ///@brief the list of fibres in the model
    std::vector<Fibre> fibres;

    /// @brief finite element solver
    PDE oxy_diff;


    /// @brief boxes in the domain
    Box ***boxes_A,***boxes_new_A;

    int maxx,maxy,maxz;
    int minx,miny,minz;
    int new_maxx,new_maxy,new_maxz;
    int new_minx,new_miny,new_minz;

    /// @brief physical size of boxes
    double box_sizex, box_sizey, box_sizez;

    /// @brief total number of cells
    unsigned int total_no_of_cells;
    /// @brief max cell number allowed
    unsigned int max_cell;


    // cell counters (per type)
    std::vector<unsigned int> totNorm,totHypo,totDead;
    // counter per phenotype at each time step
    double phenotype1_count, phenotype2_count;


    // NEW 17/0/17
    double Area_forcex;
    double Area_forcey;
    double Area_forcez;

    // ------------------------
    ///@brief maximum number of cell allowed in a box
    ///@todo to be removed when allocating box element dynamically
    unsigned int max_cell_in_box;

    /// @brief current time step
    unsigned int reloj;

    /// @brief interaction distance
    double epsilon;
    /// @brief cell initial radius [micron]
    double Ro;
    double max_radius_cell;

    unsigned int movez;
    // further model parameter (put in param class)
    ///@brief variance of random motion
    double vf;

    ///@brief displacement of cells after birth
    double birth_step;

    /// @brief energy bound for birth
    double birth_energy;

    // to check
    int initial_cells_in_box;
    int cells_counter;
    unsigned int daughter1;
    int dauther_birth;

    // ****************
    // methods
    /// @brief initialize models (using input file)
    void init(std::string f);

    /// @brief display no. of cell for each box
    void count_cells_per_box();
    /// @brief store the number of cells of different status
    void count_cells_per_type();

    // ******************
    // random generators
    /// @brief random between 0 and 1
    double aleatorio();
    /// @brief random between a and b
    double aleatorio(const double a, const double b);
    /// @brief sample from Normal distribution N(m,s)
    double box_muller(const double m, const double s);
    /// @brief sample from Normal distribution N(m,s) + bounded
    double box_muller(const double m, const double s,
              const double M, const double k);

    // ******************
    // cell functions
    /// @brief distance between 2 cells
    double DISTANCE(const Cell& c1, const Cell& c2);
    /// @brief distance between a cell and a fiber
    double DISTANCE(const Cell& cell, const Fibre& fibre);

    /// @brief change the status of cell according to O2 concentration
    void oxy_in_cell(Cell& cell);
    /// @brief search a cell in a given box (-1: cell not found)
    int search_in_box(const vector<int>& box_number,
              const unsigned int cell_name);// In solve_system
    // @brief revert cell phenotype hypoxic->normoxic (stochastic)
    void reverse_phenotype(Cell& cell);



    // ******************
    // box functions
    /// @brief allocate memories for position and boxes
    void allocate_compare_box();
    /// @brief cleans the elements of the box and update the new cells
    void update_box();
    /// @brief updates the maximum value in boxes
    void update_maximum();

    // ******** vessel related
    /// @brief set initial vessels in the system
    void set_ic_vessels();
    /// @brief determines vessels in contact with cells
    void compute_cell_vessels_contact(const int u,
              const int v,
              const int w);
    /// @brief distance between a cell and a pore
    vector<double> DISTANCE(const Cell& cell, pore* vessel);
    /// @brief point on vessel with shortest distance to cell
    float vessel_point[3];
    // ****************************************


    // ******************
    /// @brief set initial fibres in the system
    void set_ic_fibres();

    ///@brief assign different ID to different parts of the domain
    int get_sub_domain_id(double x,double y,double z);
    /// @brief set initial cells in the system
    void set_ic_cells();

    /// @brief caculates the new cells in the system
    void  cell_birth(Cell& cell);


    /// @brief compute new cell position
    void movement(const Cell& cell,
          const int u, const int v, const int w,
          const unsigned int cont_cell);
    void fibre_induced_movement(Cell& cell); // - CICELY NEW

    /// @brief calculate contact forces for cells in a given box
    void contact_forces(const int u,
                const int v,
                const int w,
                const unsigned int cells_number);

     /// @brief determines fibres in contact with cells
    void compute_cell_fibres_contact(const int u,
                const int v,
                const int w);

    /// @brief compute force acting on the cell
    void hertz(Cell& cell);

    ///@brief cell velocity correction due to surrounding fibres
    void cell_fibres_interaction(Cell& cell);

    ///@brief set the fibre_exists = 0 if the fibre is damaged
    void fibre_degradation(Cell& cell, Fibre& fibre, double cell_fibre_min_dist);

    // cell-fem coupling
    /// @brief set tetra per box
    /// @todo check this function
    void setElementsInBox(const Mesh& m);
    /// @brief check no. of element assigned to each box
    void checkElementsInBoxes();

    void compare_elements(int u, int v, int w,
              unsigned int cells_number,Mesh& _mesh);


    /// @brief time loop of coupled model
    void loop();

    /// @brief print all cells
    void printAllCells();

    /**
       @brief write cell list in a .vtk file
       @param[in] onlyCoord = 1 : write only the coordinates and radius
    */
    void writeVtk(std::string f,unsigned int onlyCoord=0);
    void writeBoxesVtk(string f);
    void writeFibresVtk(string f);
    void writeVesselsVtk(string f);

    /// @brief perform final operations
    void end();

signals:
    void plot();

private:
    ////////////// Network Attributes //////////////
    int networkSource;
    //Basic
    int Nx;
    int Ny;
    int Nz;
    std::vector<std::vector<std::vector<node*> > > tableOfNodes;
    std::vector<std::vector<std::vector<pore*> > > tableOfPoresX;
    std::vector<std::vector<std::vector<pore*> > > tableOfPoresY;
    std::vector<std::vector<std::vector<pore*> > > tableOfPoresZ;
    std::vector<pore*> tableOfAllPores;
    std::vector<node*> tableOfAllNodes;
    std::vector<particle*> tableOfParticles;
    std::vector<block*> tableOfAllBlocks;

    int totalPores;
    int totalOpenedPores;
    int totalNodes;
    int totalOpenedNodes;
    int totalBlocks;
    int totalParticles;


    double totalPoresVolume;
    double totalNodesVolume;
    double coordinationNumber;
    double minRadius;
    double maxRadius;
    double minNodeRadius;
    double maxNodeRadius;
    int radiusDistribution;
    double length;
    double degreeOfDistortion;
    double aspectRatio;
    double shapeFactor;

    double poreVolumeConstant;
    double poreVolumeExponent;
    double poreConductivityConstant;
    double poreConductivityExponent;

    double rayleighParameter;
    double triangularParameter;
    double normalMuParameter;
    double normalSigmaParameter;

    //Extracted Network
    double xEdgeLength;
    double yEdgeLength;
    double zEdgeLength;
    int maxConnectionNumber;

    double pressureIn;
    double pressureOut;
    double flow;
    double absolutePermeability;
    double porosity;

    //seed
    int seed;

    //solver
    int solverChoice;

    //perm Calc
    bool absolutePermeabilityCalculation;

    ///////////// Flow ////////////////////////
    double flowRate;
    double deltaP;
    double timeStep;
    double simulationTime;
    double recoveryTime;


    ////////////// Misc Attributes //////////////
    bool record;
    bool videoRecording;
    bool extractData;
    double extractionTimestep;  

    ////////////// fluids properties //////////////
    double plasmaViscosity;

    ////////////// Simulations //////////////
    bool drugFlowWithDiffusion;
    bool drugFlowWithoutDiffusion;
    bool particleFlow;
    bool angiogenesisTumour;
    bool angiogenesisRetina;
    bool cellProliferationSimulation;

    ////////////// Tissue Data //////////////
    bool buildTissue;
    double meshSizeX;
    double meshSizeY;
    double meshSizeZ;
    double PVT;
    double DT;
    double sigma;
    bool tissueHomo;
    bool tissueRandom;
    bool tissueCircular;
    double tissueCircularR;
    double tissueCircularX;
    double tissueCircularY;
    double tissueCircularZ;
    bool closedBoundaries;

    ////////////// Drug Data //////////////
    bool bolusInjection;
    double bolusDuration;
    bool AIFInjection;

    ////////////// Artificial Network //////////////
    bool injectParticles;
    double injectionInterval;

    ////////////// Angiogenesis On-Lattice data //////////////

    double h_sq;
    double angio_D;
    double angio_Chi;
    double angio_Delta;
    double angio_Rho;
    double angio_Eta;
    double angio_Beta;
    double angio_Gamma;
    double angio_Alpha;
    double angio_Epsilon;
    double angio_Mu;
    double angio_Psi;
    double angio_TauMax;
    bool circularTumour;
    bool linearTumour;
    bool updateBlockAttributes;
    bool phaseSeparation;
    bool branchingWSS;
    bool shuntPrevention;
    int initialTipsNumber;
    double Kp,Km,Ks,Kc;
    double Qref, QHDref, tauRef, J0;
    double decayConv, decayCond;

    ////////////// Clustering Attributes //////////////

    std::vector<cluster*> existClusters;


    ////////// Thread Management ///////////////
    bool cancel;
    bool ready;
    bool simulationRunning;

    ////////// Random generator ////////////////
    boost::random::mt19937 gen;
};

#endif // NETWORK_H
