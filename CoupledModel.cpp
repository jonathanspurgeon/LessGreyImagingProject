#include "CoupledModel.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <time.h>

/* ****************************************************************************/
using namespace std;

double PIG=3.1415926535897932384626433832795;

/*ALEATORIO: Generates random numbers between 0 and 1               */
double CoupledModel::aleatorio()
{
  double random_value = rand();
  return (random_value/RAND_MAX);
}

/*ALEATORIO_moves: Generates random numbers between a and b        */
double CoupledModel::aleatorio(const double a, const double b)
{
  double random_value = 0.;
 
  if (a<=b) {
    random_value = b - (b-a)*this->aleatorio();
  } else {
    random_value = 1. - 2.*this->aleatorio();
  }
  return(random_value);
}

/* *****************************************************************************
  Routine for generating gaussian random numbers N(m,s)
  ***************************************************************************** */
double CoupledModel::box_muller(const double m, const double s)	
{				        
  double w, y1;
  double x1, x2;
  static double y2;
  static int use_last = 0;
  if (use_last)	{	        /* use value from previous call */
    y1 = y2;
    use_last = 0;
  } else {
    do {
      x1 = 2.0 * this->aleatorio() - 1.0;
      x2 = 2.0 * this->aleatorio() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    
    w = sqrt( (-2.0 * log( w ) ) / w );
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
  
  return( m + y1 * s );
}

/** ************************************************************************ 
 Routine for generating gaussian random numbers N(m,s) + bound
 *****************************************/
double CoupledModel::box_muller(const double m, const double s, 
				const double MAX, const double MIN)	
{	
  double val = this->box_muller(m, s);
			        
  if (val < MIN)
    return MIN;
  else if (val > MAX)
    return MAX;
  else  
    return val;
}

/** ************************************************************************ 
 Routine for determining distance between cells
****************************************************************************/
double CoupledModel::DISTANCE(const Cell& c1, const Cell& c2)
{
  double dist = 
    pow(c1.position[0] - c2.position[0],2)+
    pow(c1.position[1] - c2.position[1],2)+
    pow(c1.position[2] - c2.position[2],2);

  return(sqrt(dist));
}
/****************************************************************************/
//Routine for determining distance between cell and fibre
/**
   @brief 
   The distance between cell and fibre is computed using the
   scalar product Y = <c_to_f,fibre_vector> where
   c_to_f = cell_center - fibre_start
   fibre_vector = fibre.length*fibre.direction
   -> if Y<0, then the cell is closer to fibre.start
   -> if Y>fibre.length^2, then the cell is closer to fibre.end
   -> otherwise, we take the projection of cell center onto the fibre
 */
double CoupledModel::DISTANCE(const Cell& cell, const Fibre& fibre)
{
  double c_to_f[3]; // cell-to-fibre vector
  double c_to_f_length_2; // length of c_to_f squared
  double c_to_f_dot_v_f; // scalar product c_to_f * (f_end-f_start)
  
  c_to_f_length_2 = 0.;
  c_to_f_dot_v_f = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_to_f[i] = cell.position[i]-fibre.start[i];
    c_to_f_length_2 = c_to_f_length_2 + c_to_f[i]*c_to_f[i];
    c_to_f_dot_v_f = c_to_f_dot_v_f + c_to_f[i]*fibre.length*fibre.direction[i];
  }
  
  if (c_to_f_dot_v_f < 0.) {
    // in this case the cell is closest to the fibre start point
    //cout << "the cell is closest to the start of the fibre" << endl;
    return sqrt(c_to_f_length_2);
    
    } else if (c_to_f_dot_v_f > fibre.length*fibre.length) {
    // in this case the cell is closest to the fibre end point
    //cout << "the cell is closest to the end of the fibre" << endl;
    double distance = 0.;
    for (unsigned int i=0; i<3; i++) {
      distance = distance +
	(cell.position[i]-(fibre.start[i]+fibre.length*fibre.direction[i]))*
	(cell.position[i]-(fibre.start[i]+fibre.length*fibre.direction[i]));
    }
    return sqrt(distance);
    
    } else {
    // in this case the cell is closest to a point along the fibre
    //cout << "the cell is closest a point along the fibre" << endl;
    double c_to_f_length_cos_alpha_2 = c_to_f_dot_v_f*c_to_f_dot_v_f/
      (fibre.length*fibre.length);
    return sqrt( c_to_f_length_2 - c_to_f_length_cos_alpha_2);
    }
  
}
/****************************************************************************/

/****************************************************************************/
//Routine for determining distance between cell and vessel
/**
   @brief 
   The distance between cell and vessel is computed using the
   scalar product Y = <c_to_v,vessel_vector> where
   c_to_v = cell_center - vessel_start
   vessel_vector = vessel.length*vessel.direction
   -> if Y<0, then the cell is closer to vessel.start
   -> if Y>fibre.length^2, then the cell is closer to vessel.end
   -> otherwise, we take the projection of cell center onto the fibre
 */
double CoupledModel::DISTANCE(const Cell& cell, const Vessel& vessel)
{
  double c_to_v[3]; // cell-to-vessel vector
  double c_to_v_length_2; // length of c_to_v squared
  double c_to_v_dot_v_v; // scalar product c_to_v * (v_end-v_start)
  
  c_to_v_length_2 = 0.;
  c_to_v_dot_v_v = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_to_v[i] = cell.position[i]-vessel.ves_start[i];
    c_to_v_length_2 = c_to_v_length_2 + c_to_v[i]*c_to_v[i];
    c_to_v_dot_v_v = c_to_v_dot_v_v + c_to_v[i]*vessel.ves_length*vessel.ves_direction[i];
  }
  
  if (c_to_v_dot_v_v < 0.) {
    return sqrt(c_to_v_length_2);    
    } else if (c_to_v_dot_v_v > vessel.ves_length*vessel.ves_length) {
    double distance = 0.;
    for (unsigned int i=0; i<3; i++) {
      distance = distance +
	(cell.position[i]-(vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i]))*
	(cell.position[i]-(vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i]));
    }
    return sqrt(distance);
  } else {
    double c_to_v_length_cos_alpha_2 = c_to_v_dot_v_v*c_to_v_dot_v_v/
      (vessel.ves_length*vessel.ves_length);
    return sqrt( c_to_v_length_2 - c_to_v_length_cos_alpha_2);
    }
  
}
/****************************************************************************/

/****************************************************************************/
//Routine for determining point on vessel which provides minimum distance
void CoupledModel::VESSELPOINT(const Cell& cell, const Vessel& vessel)
{
  double c_v[3]; // cell-to-vessel vector
  double c_v_length_2; // length of c_to_v squared
  double c_v_dot_v_v; // scalar product c_to_v * (v_end-v_start)
  
  c_v_length_2 = 0.;
  c_v_dot_v_v = 0.;
  for (unsigned int i=0; i<3; i++) {
    c_v[i] = cell.position[i]-vessel.ves_start[i];
    c_v_length_2 = c_v_length_2 + c_v[i]*c_v[i];
    c_v_dot_v_v = c_v_dot_v_v + c_v[i]*vessel.ves_length*vessel.ves_direction[i];
  }
  
  if (c_v_dot_v_v < 0.) {
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i];
    }
    //cout << " case 1 " << vessel_point[0] << " "  << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } else if (c_v_dot_v_v > vessel.ves_length*vessel.ves_length) {
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i]+vessel.ves_length*vessel.ves_direction[i];
    }
    //cout << " case 2 " << vessel_point[0] << " " << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } else {
    double dist_to_vessel_point = c_v_dot_v_v/(vessel.ves_length);
    for (unsigned int i=0; i<3; i++) {
      this->vessel_point[i] = vessel.ves_start[i]+dist_to_vessel_point*vessel.ves_direction[i];
    }
    //cout << " case 3 " << vessel_point[0] << " "  << vessel_point[1] << " " <<  vessel_point[2] << endl;
  } 
}
/****************************************************************************/

/* initialize the class (read input file, set parameters) */
void CoupledModel::init(string f)
{
  // read input parameters
  this->input_file_name = f;
  //read an input file with GetPot
  params.readFile(f);

  // initialize simulation
  // set box size
  this->box_sizex=params.lattice_length_x/(params.boxesx+0.0);
  this->box_sizey=params.lattice_length_y/(params.boxesy+0.0);
  this->box_sizez=params.lattice_length_z/(params.boxesz+0.0);

  // random seed
  srand(static_cast<unsigned>(time(NULL)));

  this->vf = params.variance_motion;

  //allocate memory
  this->allocate_compare_box();
  
  // initial values for range of occupied boxes
  this->maxx=0;
  this->maxy=0;
  this->maxz=0;
  this->minx=10000;
  this->miny=10000;
  this->minz=10000;
  
  this->new_maxx=0;
  this->new_maxy=0;
  this->new_maxz=0;  
  this->new_minx=10000;
  this->new_miny=10000;
  if (params.dimension == 3)  this->new_minz=10000;
  else this->new_minz=0;
  
  /// @brief: max cell per box: depends on cell volume and box size
  this->epsilon = 10;
  this->Ro = this->epsilon/2.;
  this->max_radius_cell = this->Ro;
  double cell_estimated_volume = 4./3.*PIG*pow(this->max_radius_cell,3.);
  double box_volume =  this->box_sizex*this->box_sizey*this->box_sizez;
  if (params.dimension == 2) {
    cell_estimated_volume = 2.*PIG*pow(this->max_radius_cell,2.);
    box_volume =  this->box_sizex*this->box_sizey;
  }
  double cell_compressibility = 6.; // compressibility of cells
  this->max_cell_in_box = box_volume/cell_estimated_volume*cell_compressibility;
  if (params.dimension == 2) this->max_cell_in_box *= 2;
  cout << " === max. number of cell per box allowed: " <<  this->max_cell_in_box << endl;
  this->max_cell = params.max_cell;

  cout << " set initial conditions for fibres" << endl;
  this->set_ic_fibres();
  //this->place_fibres_in_boxes();

  cout << " set initial conditions for vessels" << endl;
  this->set_ic_vessels();

  cout << " set initial conditions for cells" << endl;
  // set IC for cells
  this->set_ic_cells();

  // updates the maximum value in boxes
  this->update_maximum(); 
  this->update_box();
  

  // print initial configuration on screen
  if ( this->params.verbose > 0) {
    cout << " === number of boxes (" 
	 << params.dimension << "D): " 
	 << params.boxesx << " x "
	 << params.boxesy << " x "
	 << params.boxesz << endl;

    cout << " === initial configuration: " << endl;
    for(unsigned int k=0; k<params.boxesx;k++) {
      for(unsigned int l=0; l<params.boxesy;l++) {
	for(unsigned int n=0; n<params.boxesz;n++) {
	  if (this->boxes_A[k][l][n].cells.size()>0) {
	    cout << " -> box " << k << "," << l << "," << n
		 << ": " << this->boxes_A[k][l][n].cells.size()
		 << " cell(s) " << endl;
	    cout << "type: ";
	  for (unsigned int j=0; j<this->boxes_A[k][l][n].cells.size(); j++) {
	    cout << this->boxes_A[k][l][n].cells[j].type << " ";
	  }
	  cout << endl;
	  }
	}
      }
    }	
    cout << endl;
  }

  ///@todo check the birth step: too small? depends on dimension?
  this->birth_step=.9/pow(2,1./3.) * this->Ro; // cell displacement due to newbirth
  
  if (params.dimension==2){
    this->movez = 0;
  } else if (params.dimension==3) {
    this->movez = 1;
  } else {
    cout << " ** ERROR: dimension = " << params.dimension
	 << " not supported. " << endl;
    exit(1);
  }

  // to check
  initial_cells_in_box = 0;
  cells_counter = 0;
  daughter1 = 1;

}

/*****************************************************************************/
/* allocate_compare_box : allocates memory for name, distancia and posiciones*/
/*****************************************************************************/
void CoupledModel::allocate_compare_box()
{

  // dynamic allocated memory for cells in each box
  this->boxes_A = new Box** [params.boxesx];
  this->boxes_new_A = new Box** [params.boxesx];
  
  for(unsigned int i=0; i<params.boxesx; i++) {
    this->boxes_A[i] = new Box* [params.boxesy]; 
    this->boxes_new_A[i] = new Box* [params.boxesy]; 
  } 
  for(unsigned int i=0; i<params.boxesx; i++) {
    for(unsigned int j=0; j<params.boxesy; j++) {
      this->boxes_A[i][j] = new Box [params.boxesz]; 
      this->boxes_new_A[i][j] = new Box [params.boxesz]; 
    } 
  }

  for(unsigned int i=0; i<params.boxesx; i++) {
    for(unsigned int j=0; j<params.boxesy; j++) {
      for(unsigned int k=0; k<params.boxesz; k++) {
	this->boxes_A[i][j][k].cells.resize(0);
	this->boxes_new_A[i][j][k].cells.resize(0);
	this->boxes_A[i][j][k].v_triangles.resize(0);
	this->boxes_new_A[i][j][k].v_triangles.resize(0);
	this->boxes_A[i][j][k].p_fibres.resize(0);
	this->boxes_new_A[i][j][k].p_fibres.resize(0);
      }
    } 
  }

}
/********************************************************************/
// place initial cells in the system
void CoupledModel::set_ic_cells()
{
  Cell cell;
  int newbox[3];
  this->total_no_of_cells = params.n_initial_cells;
  this->phenotype1_count = 0.0;
  this->phenotype2_count = 0.0;
  //Initialise the boxes cells number:
  for(unsigned int l=0; l < this->total_no_of_cells ; l++) {

    //Position
    cell.position[0]=params.ic_cell_x[l];
    cell.position[1]=params.ic_cell_y[l];
    cell.position[2]=params.ic_cell_z[l];

    cell.position_old[0]=params.ic_cell_x[l];
    cell.position_old[1]=params.ic_cell_y[l];
    cell.position_old[2]=params.ic_cell_z[l];

    //velocity
    cell.polarity[0]=params.polarity_x[l];
    cell.polarity[1]=params.polarity_y[l];
    cell.polarity[2]=params.polarity_z[l];
    
    for (unsigned int j=0; j<3; j++){
      cell.vel[j]=0.;
    }

    //Cell characteristics
    cell.name=l;  
    cell.contacts=0;
    cell.mother_name = -1;
    cell.birthday = this->reloj;
    // initialising the cell contact area
    //cell.contact_area_old=0.0;
    //cell.variation_area=0.0;
    
    cell.type=1;
    cell.radius=this->Ro;

    cell.phenotype=params.threshold_hypo;
    cell.phenotype_counter = 0;

    cell.polarised = 0; // CICELY NEW

    // oxygen concentration
    cell.O2 = 100.; // !! TODO: this should be given by the diffusion solver !!
    cell.dxO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dyO2 = 0.; // !! TODO: this should be given by the diffusion solver !!
    cell.dzO2 = 0.; // !! TODO: this should be given by the diffusion solver !!

    // -----
    cell.type = 1;
    cell.hypoxic_count = 0;

    // set the interaction phenotype
    cell.interaction_phenotype = params.ic_phenotype[l];

    // set whether the cell is a follower or leader
    cell.follow_lead = params.ic_follower_leader[l];

    // CICELY NEW - Change the adhesion_value 
    cell.adhesion=params.adhesion_value[cell.interaction_phenotype];
    
    if (cell.interaction_phenotype==0)
    {
      this->phenotype1_count = this->phenotype1_count+1;
    }
    else
    {
      this->phenotype2_count = this->phenotype2_count+1;
    }
    
    //Cell box  
    int u=(int)(floor(cell.position[0]/this->box_sizex));
    int v=(int)(floor(cell.position[1]/this->box_sizey));
    int w=(int)(floor(cell.position[2]/this->box_sizez));
    
    cell.box[0]=u;
    cell.box[1]=v;
    cell.box[2]=w;

    cell.new_box[0]=u;
    cell.new_box[1]=v;
    cell.new_box[2]=w;
    cout << " cell " << l << " placed in box: "
	 << u << " " << v << " " << w << endl;
    cout << " position: " << cell.position[0] << " " 
	 << cell.position[1] << " " << cell.position[2] <<
      " phenotype " << cell.interaction_phenotype << endl;


    this->boxes_A[u][v][w].cells.push_back(cell);
    this->boxes_new_A[u][v][w].cells.push_back(cell);

    // box of the new cell
    newbox[0]=u; 
    newbox[1]=v;
    newbox[2]=w;
   

    // compute the extrema of occupied boxes
    if(this->maxx<u && u<(int) params.boxesx) {
      this->new_maxx=u;	
    }
    if(this->maxy<v && v<(int) params.boxesy) {
      this->new_maxy=v;	
    }
    
    if(this->maxz<w && w<(int) params.boxesz) {
      this->new_maxz=w;	
    }	
    
    if(this->minx>u && u>0) {
      this->new_minx=u;	
    }
    
    if(this->miny>v && v>0) {
      this->new_miny=v;	
    }
    
    if(this->minz>w && w>0) {
      this->new_minz=w;	
    }	
    
    /// @todo can we set newbox=[u,v,w]?
    /// @todo what is the differnce between maxx and new_maxx?
    if(this->maxx<newbox[0] && newbox[0]<(int) params.boxesx) {
      this->maxx=newbox[0];	
    }
    if(this->maxy<newbox[1]&& newbox[1]<(int) params.boxesy) {
      this->maxy=newbox[1];
    }
    if(this->maxz<newbox[2]&& newbox[2]<(int) params.boxesz) {
      this->maxz=newbox[2];
    }	
    
    if(this->minx>newbox[0]&& newbox[0]>=0 ) {
      this->minx=newbox[0];	
    }
    if(this->miny>newbox[1]&& newbox[1]>=0 ) {
      this->miny=newbox[1]	;
    }
    if(this->minz>newbox[2]&& newbox[2]>=0 ) {
      this->minz=newbox[2];	
    }	

  }

  cout << " phenotype 1 - " << this->phenotype1_count << " cells " << endl;
  cout << " phenotype 2 - " << this->phenotype2_count << " cells " << endl;
  
  //cout << " -- occupied box with smallest coordinate (min box):";
  //cout << this->minx << " " << this->miny << " " << this->minz << endl;
  //cout << " -- occupied box with largest coordinate (max box):";
  //cout << this->maxx << " " << this->maxy << " " << this->maxz << endl;
 
}

/****************************************************************************/
// place fibres in the system 
void CoupledModel::set_ic_fibres()
{
  Fibre fibre;
  double fibre_x, fibre_y, fibre_z;

  // we divide the domain in subdomains
 
  for (unsigned int sub_domain=0; sub_domain<params.n_sub_domain; sub_domain++){

  for (int ll=0; ll < this->params.n_initial_fibres[sub_domain] ; ll++){
    fibre.fibre_name=ll;
    fibre.fibre_exists=1.;
    
    // for each fibre randomly position one end point within the domain
    /// @brief x co-ordinate for each fibre
    fibre.start[0] = aleatorio(params.x_start[sub_domain],params.x_end[sub_domain]);
    /// @brief y co-ordinate for each fibre
    fibre.start[1] = aleatorio(0.0,params.lattice_length_y); 
    /// @brief z co-ordinate for each fibre
    if (params.dimension == 3) {
      fibre.start[2] = aleatorio(0.0,params.lattice_length_z); 
    } else {
      fibre.start[2] = 0.;
    }
  
    
    // for each fibre randomly assign angles for direction
    if (params.fibre_orientation_distribution[sub_domain]==0){
      
      fibre.phi = aleatorio(0.0,2.*PIG);
      
      if (params.dimension == 3) fibre.theta = aleatorio(0.,PIG);
      else fibre.theta = PIG/2.;
      
    }
    
    if (params.fibre_orientation_distribution[sub_domain]==1){
      // for each fibre randomly assign angles for direction
      /// @brief create a normally distributed angle
      fibre.phi = box_muller(params.fibre_orientation_mean_phi[sub_domain],
			       params.fibre_orientation_variance_phi[sub_domain]); 
      fibre.theta = box_muller(params.fibre_orientation_mean_theta[sub_domain],
			     params.fibre_orientation_variance_theta[sub_domain]);
    }
    
    /// @brief for each fibre randomly assign length
    fibre.length = box_muller(params.fibre_length_mean[sub_domain],
			      params.fibre_length_variance[sub_domain]);
    
    /// @brief for each fibre randomly align throughout the domain
    fibre_x = fibre.length*sin(fibre.theta)*cos(fibre.phi);
    fibre_y = fibre.length*sin(fibre.theta)*sin(fibre.phi);
    fibre_z = fibre.length*cos(fibre.theta);
    /// @brief x direction for each fibre
    fibre.direction[0] = fibre_x/(sqrt(pow(fibre_x,2)+pow(fibre_y,2)+pow(fibre_z,2)));
    /// @brief y direction for each fibre
    fibre.direction[1] = fibre_y/(sqrt(pow(fibre_x,2)+pow(fibre_y,2)+pow(fibre_z,2)));
    /// @brief z direction for each fibre
    fibre.direction[2] = fibre_z/(sqrt(pow(fibre_x,2)+pow(fibre_y,2)+pow(fibre_z,2)));
    if (params.dimension == 2) {
      fibre_z = 0.;
      fibre.direction[2] = 0.;
    }
    /*
   /// @brief Determine whether end point of fibre lies within Domain - if not change the direction
   /// @todo is this an acceptable thing to do? Does it affect the distribution?
   double latice_length[3];
   latice_length[0] = params.lattice_length_x;
   latice_length[1] = params.lattice_length_y;
   latice_length[2] = params.lattice_length_z;
   //cout << " " << latice_length[0] << " " << latice_length[1] << " " << latice_length[2] << endl;
   for (unsigned int j=0; j<=2; j++)
   {
   while (fibre.start[j]+fibre.length*fibre.direction[j]>latice_length[j] || fibre.start[j]+fibre.length*fibre.direction[j]<0)
   {
   //fibre.direction[j]=-fibre.direction[j];
   }
   }*/
      
    //cout << " fibre " << ll << " placed at: " << fibre.start[0] << "," << fibre.start[1] << "," << fibre.start[2] << endl;
    //cout << " length: " << fibre.length << " direction " << fibre.direction[0] << "," << fibre.direction[1] << "," << fibre.direction[2]  << endl;
    
    // This part of the code checks the position of initial fibres in relation to initial cells
    // if any fibre is too close to an initial cell then new fibre parameters are calculated
    /*if (params.fibre_make_gap==1){
      this->total_no_of_cells = params.n_initial_cells;
      double f_centre[3];
      
      for (unsigned int i=0; i<=2; i++){
	f_centre[i] = fibre.start[i]+0.5*fibre.length*fibre.direction[i];
      }
      
      double c_to_f_centre;
      double cell_pos[3];
      
      for (unsigned int l=0; l < this->total_no_of_cells ; l++){
	c_to_f_centre = 0;
	cell_pos[0] = params.ic_cell_x[l];
	cell_pos[1] = params.ic_cell_y[l];
	cell_pos[2] = params.ic_cell_z[l];
	
	for (unsigned int jj=0; jj<=2; jj++){
	  c_to_f_centre = c_to_f_centre+pow(cell_pos[jj]-f_centre[jj],2);
	}
	
	double cell_fibre_dist = sqrt(c_to_f_centre);	
	double ref_dist = this->Ro + fibre.length;
	//cout << ref_dist << endl;
	
	if (cell_fibre_dist < ref_dist){
	  //cout << "fibre parameters need updating" << endl;
	  // for each fibre randomly position one end point within the domain
	  /// @brief x co-ordinate for each fibre
	  fibre.start[0] = aleatorio(0,params.lattice_length_x); 
	  /// @brief y co-ordinate for each fibre
	  fibre.start[1] = aleatorio(0,params.lattice_length_y); 
	  /// @brief z co-ordinate for each fibre
	  fibre.start[2] = aleatorio(0,params.lattice_length_z); 
	  
	  if (params.fibre_orientation_distribution==0){
	    // for each fibre randomly assign angles for direction
	    /// @brief create random angle between 0 and 2pi
	    fibre.theta = aleatorio(0,2*PIG); 
	    fibre.phi = aleatorio(0,2*PIG); 
	  }
	  
	  if (params.fibre_orientation_distribution==1){
	    // for each fibre randomly assign angles for direction
	    /// @brief create a normally distributed angle
	    fibre.theta = box_muller(params.fibre_orientation_mean_phi,params.fibre_orientation_variance_phi); 
	    fibre.phi = box_muller(params.fibre_orientation_mean_theta,params.fibre_orientation_variance_theta);
	  }
	  
	  /// @brief for each fibre randomly assign length
	  fibre.length = box_muller(params.fibre_length_mean,params.fibre_length_variance);
	  
	  /// @brief for each fibre randomly align throughout the domain
	  fibre_x = fibre.length*sin(fibre.theta)*cos(fibre.phi);
	  fibre_y = fibre.length*sin(fibre.theta)*sin(fibre.phi);
	  fibre_z = fibre.length*cos(fibre.theta);
	  /// @brief x direction for each fibre
	  fibre.direction[0] = fibre_x/(sqrt(pow(fibre_x,2)+pow(fibre_y,2)+pow(fibre_z,2)));
	  /// @brief y direction for each fibre
	  fibre.direction[1]= fibre_y/(sqrt(pow(fibre_x,2)+pow(fibre_y,2)+pow(fibre_z,2)));
	  /// @brief z direction for each fibre
	  fibre.direction[2]= fibre_z/(sqrt(pow(fibre_x,2)+pow(fibre_y,2)+pow(fibre_z,2)));
	  
	  /// @brief Determine whether end point of fibre lies within Domain - if not change the direction
	  /// @todo is this an acceptable thing to do? Does it affect the distribution?
	  double latice_length[3];
	  latice_length[0] = params.lattice_length_x;
	  latice_length[1] = params.lattice_length_y;
	  latice_length[2] = params.lattice_length_z;
	  
	  for (unsigned int j=0; j<=2; j++){
	    while (fibre.start[j]+fibre.length*fibre.direction[j]>latice_length[j] || fibre.start[j]+fibre.length*fibre.direction[j]<0){
	      fibre.direction[j]=-fibre.direction[j];
	    }
	  }
	  
	  //cout << " fibre " << ll << " placed at: " << fibre.start[0] << "," << fibre.start[1] << "," << fibre.start[2] << endl;
	  //cout << " length: " << fibre.length << " direction " << fibre.direction[0] << "," << fibre.direction[1] << "," << fibre.direction[2]  << endl;
	}
      }
      }*/
    fibres.push_back(fibre);
  }
  }
}

// ***** initialise vessels ******* NEW 5/06/17
void CoupledModel::set_ic_vessels()
{
  Vessel vessel;

  /// @TODO re-write for other dimensions?

  for (int v=0; v < this->params.n_initial_vessels ; v++){
    vessel.vessel_name = v;
    vessel.ves_radius = params.vessel_radius[v];

    /// @brief vessel start position
    vessel.ves_start[0] = params.vessel_startx[v];
    vessel.ves_start[1] = params.vessel_starty[v];
    vessel.ves_start[2] = params.vessel_startz[v];
    
    /// @brief vessel length
    vessel.ves_length = params.vessel_length[v];
   
    /// @brief vessel direction
    vessel.ves_direction[0] = params.vessel_directionx[v];
    vessel.ves_direction[1] = params.vessel_directiony[v];
    vessel.ves_direction[2] = params.vessel_directionz[v];
   
    vessels.push_back(vessel);
  }
}

  // ***********************************************************
void CoupledModel::update_maximum() 
{

  bool stopRun = false;
  if (this->new_maxx < this->new_minx) {
    cout << " ** WARNING: new_maxx " << new_maxx << ", new_minx " << new_minx << endl;
    cout << " ** WARNING: maxx " << maxx << ", minx " << minx << endl;
    stopRun = true;
  }
  if (this->new_maxy < this->new_miny) {
    cout << " ** WARNING: new_maxy " << new_maxy << ", new_miny " << new_miny << endl;
    cout << " ** WARNING: maxy " << maxy << ", miny " << miny << endl;
    stopRun = true;
  }
  if (params.dimension == 3) {
    if (this->new_maxz < this->new_minz) {
      cout << " ** WARNING: new_maxz " << new_maxz << ", new_minz " << new_minz << endl;
      cout << " ** WARNING: maxz " << maxz << ", minz " << minz << endl;
      stopRun = true;
    }
  }
  if (stopRun) exit(1);
  this->maxx = this->new_maxx;
  this->maxy = this->new_maxy;
  this->maxz = this->new_maxz;
  
  this->minx = this->new_minx;
  this->miny = this->new_miny;
  this->minz = this->new_minz;
}

// ***********************************************************
/* Cleans the elements of the box and update the new cells */
void CoupledModel::update_box()
{

  for(int u=this->minx; u<=this->maxx; u++) {
    for(int v=this->miny; v<=this->maxy; v++) {
      for(int w=this->minz; w<=this->maxz; w++) {

	// clear cell contacts
	for(unsigned int j=0; j<this->boxes_new_A[u][v][w].cells.size(); j++)   {
	   this->boxes_new_A[u][v][w].cells[j].clear_contacts();
	}
	
	// boxes.cells = boxes_new.cells
	this->boxes_A[u][v][w].cells = boxes_new_A[u][v][w].cells;

	//for(unsigned int j=0; j<this->boxes_A[u][v][w].cells.size(); j++)   {
	//  // clear cell contacts
	//  this->boxes_A[u][v][w].cells[j].clear_contacts();
	//}


	// clear new box arrays
	this->boxes_new_A[u][v][w].cells.clear();
	// new: free also the used memory
	vector<Cell> swap(this->boxes_new_A[u][v][w].cells);
	//this->boxes_new_A[u][v][w].cells.resize(0);
	  
      }  
    }  
  }
  
}

/***********************************
 For each box, find the elements with barycenter inside the box
 *********************************** */
void CoupledModel::setElementsInBox(const Mesh& _mesh)
{
  if (_mesh.dim==2) {

    // triangular mesh
    for(unsigned int l=0; l<_mesh.nTria; l++) {
      // barycenter of triangle
      double xT = (_mesh.xp[_mesh.tria[3*l]-1] +  _mesh.xp[_mesh.tria[3*l+1]-1] +  
		   _mesh.xp[_mesh.tria[3*l+2]-1])/3.;
      double yT = (_mesh.yp[_mesh.tria[3*l]-1] +  _mesh.yp[_mesh.tria[3*l+1]-1] +  
		   _mesh.yp[_mesh.tria[3*l+2]-1])/3.;
      double zT = params.lattice_length_z/2.;

      /// @todo floor() not needed here
      int u= (int)(floor( xT/this->box_sizex ));
      int v= (int)(floor( yT/this->box_sizey ));
      int w= (int)(floor( zT/this->box_sizez ));

      this->boxes_A[u][v][w].v_triangles.push_back(l);
      
    }   

  } else {

    // tetrahedral mesh
    for(unsigned int l=0; l<_mesh.nTetra; l++) {
      double xT = (_mesh.xp[_mesh.tetra[4*l]-1] +  _mesh.xp[_mesh.tetra[4*l+1]-1] +  
		   _mesh.xp[_mesh.tetra[4*l+2]-1] + _mesh.xp[_mesh.tetra[4*l+3]-1])/4.;
      
      double yT = (_mesh.yp[_mesh.tetra[4*l]-1] +  _mesh.yp[_mesh.tetra[4*l+1]-1] +  
		   _mesh.yp[_mesh.tetra[4*l+2]-1] + _mesh.yp[_mesh.tetra[4*l+3]-1])/4.;
      
      double zT = (_mesh.zp[_mesh.tetra[4*l]-1] +  _mesh.zp[_mesh.tetra[4*l+1]-1] +  
		   _mesh.zp[_mesh.tetra[4*l+2]-1] + _mesh.zp[_mesh.tetra[4*l+3]-1])/4.;

      /// @todo floor() not needed here
      int u= (int)(floor( xT/this->box_sizex ));
      int v= (int)(floor( yT/this->box_sizey ));
      int w= (int)(floor( zT/this->box_sizez ));
      
      this->boxes_A[u][v][w].v_triangles.push_back(l);

    }
  }
  
}

// good for debugging: check how many elements
// have been assigned to a box
void CoupledModel::checkElementsInBoxes()
{

  for(unsigned int k=0; k<params.boxesx; k++) {
    for(unsigned int l=0; l<params.boxesy; l++) {
      for(unsigned int n=0; n<params.boxesz; n++) {
	cout << "box [" << k << " , "  << l << " , "  << n << "].tetra = ";
	for (unsigned int ij=0; ij< this->boxes_A[k][l][n].v_triangles.size(); ij++) {
	  cout << this->boxes_A[k][l][n].v_triangles[ij] << " ";
	}
	cout << " (tot : "<<  this->boxes_A[k][l][n].v_triangles.size() << ")"<< endl;
      }
    }
  }

}

/*
***************************************************************************
change the status of cell according to O2 concentration
***************************************************************************
*/
/// @todo add state change in the Cell variable
void CoupledModel::oxy_in_cell(Cell& cell) {

  if(cell.type != 3) {

    //Process of normoxic-> hypoxic conversion	    
    if(cell.O2 >= params.threshold_death && cell.O2 < cell.phenotype) { 

      cell.type = 2;
      cell.hypoxic_count++;
    }

    //Process of hypoxic->anoxic conversion
    if(cell.O2 < params.threshold_death) {
      cell.type =3;	
    }	
    
    if(cell.hypoxic_count > params.time_death) {
      cell.type =3;	
    }
  }		
}

/************************************************************************************
 revert cell phenotype hypoxic->normoxic if cell has enough oxygen
 ************************************************************************************/
void CoupledModel::reverse_phenotype(Cell& cell) {

  /// @todo add a parameter for regulating reversion
  double time_steps_per_day = 60.*24./(params.time_step+0.0);
  float valor= this->aleatorio();  
  if( valor < 1./(2.*time_steps_per_day)) {

    cell.type = 1;
    cout << " ** reverse_phenotype: reverting cell " << cell.name << " to normoxic ** " << endl;

  }
  
}

/*************************************************************************************
Caculates the new cells in the system for different phenotypes  
*************************************************************************************/
void CoupledModel::cell_birth(Cell& cell)
{
  // coordinate of the box occupied by the cell
  //int u=cell.box[0];
  //int v=cell.box[1];
  //int w=cell.box[2];

  //int position_in_box_array = search_in_box(cell.box,cell.name); 
  
  /// @todo specify max n. of contacts in input file
  /// @todo put the conditions in a separate bool function
  //this->b_energy = params.birth_energy[cell.interaction_phenotype];

  double E = fmax(params.alpha_YoungM[0],params.alpha_YoungM[1]);
  double nu = fmax(params.alpha_PoissonNo[0],params.alpha_PoissonNo[1]);
  double eff_modulus = E/(2.0*(1-nu*nu));
  double eff_radius = this->Ro*this->Ro/(2*this->Ro);
  double be_d = params.be_displacement;
  double force_rep = 4./3. * eff_modulus * sqrt(eff_radius) * pow(be_d,be_d);
  this->birth_energy = params.be_multiplier * force_rep / (2 * PIG * be_d);
  //cout << "birth energy is " << birth_energy << endl;
  
  if(cell.contacts<16 && // birth allowed only if n.contact < threshold
     cell.type==1 && // type = 1: normoxic cells
     cell.radius>0.99*max_radius_cell // birth allowed only if radius is big enough
     && cell.energy<this->birth_energy// birth allowed only if contact energy is below a threshold
     ) {

    double birth_probability = aleatorio();

    if(birth_probability < params.alpha_birthrate[cell.interaction_phenotype]) {

      if (params.verbose>1) {
	cout << " ** BIRTH: new number of cells " <<  this->total_no_of_cells+1 << " ** " << endl;	
      }

      // 3D: 2.rnew^3 = rold^3
      // 2D: 2.rnew^2 = rold^2
      double new_radius=cell.radius/pow(2.,1./params.dimension);
      //We reset the radius and intracellular concentrations of the mother cell  
      //this->boxes_A[cell.box[0]][cell.box[1]][cell.box[2]].cells[position_in_box_array].radius=new_radius;
      cell.radius = new_radius;
      
      // compute the position of the daugther cell
      double newpositionx,newpositiony,newpositionz;

      if(cell.contacts!=0) {
	//We take the preferred position calculated by the neighbours
	//Noise is necessary in order not to repeat births.
	// todo...
	//newpositionx=newpositionx + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
	//newpositiony=newpositiony + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
	//newpositionz=newpositionz + this->birth_step + box_muller(0,0.3,cell.radius,-cell.radius);
      }

      //If no cells around the we choose a position at random
      double phi = box_muller(0,6.28);
      double theta = box_muller(0,3.14);
      if (params.dimension == 2) theta = acos(-1.0)/2.;
      newpositionx = cell.position[0] + this->birth_step*sin(theta)*cos(phi);
      newpositiony = cell.position[1] + this->birth_step*sin(theta)*sin(phi);
      newpositionz = cell.position[2] + this->birth_step*cos(theta);
          
      //check errors:
      if(newpositionx==cell.position[0] && 
	 newpositiony==cell.position[1] 
	 && newpositionz==cell.position[2]) {
	cout << " *** error in CoupledModel::birthday_phenotype : "
	     << " same position for daughter cell, file "
	     << __FILE__ << " line " << __LINE__ << endl;
	exit(1);
      }
    
      if ((newpositionx<0)||(newpositionx>params.lattice_length_x)||
	  (newpositiony<0)||(newpositiony>params.lattice_length_y)||
	  (newpositionz<0)||(newpositionz>params.lattice_length_z)) {
	cout << " ** error: the new cell created from " << cell.name << " is moving out of the domain (not supported) " << endl;
	exit(1);
      }

      
      /*
      // use periodic BC
      if (newpositionx<0) 
	newpositionx=params.lattice_length_x + newpositionx;
      if (newpositionx>params.lattice_length_x) 
	newpositionx=newpositionx - params.lattice_length_x;
      if (newpositiony<0) 
	newpositiony=params.lattice_length_y + newpositiony;
      if (newpositiony>params.lattice_length_y) 
	newpositiony=newpositiony - params.lattice_length_y;
      if (newpositionz<0) 
	newpositionz=params.lattice_length_z + newpositionz;
      if (newpositionz>params.lattice_length_z) 
	newpositionz=newpositionz - params.lattice_length_z;
      */
    
      // coordinate of the box of the daughter cell
      int c1=(int)(floor(newpositionx/this->box_sizex)); 
      int c2=(int)(floor(newpositiony/this->box_sizey)); 
      int c3=(int)(floor(newpositionz/this->box_sizez)); 
      
      if (this->boxes_A[c1][c2][c3].cells.size()>=max_cell_in_box) {
	cout << " *** ERROR (time step " << reloj
	     << "): !! could not create new cell !! *** "<< endl;
	cout << " *** in box " << c1 << "," << c2 << "," << c3 << endl;
	cout << " *** I found " << this->boxes_A[c1][c2][c3].cells.size() << " cells " << endl;
	cout << " *** (max_cell_in_box = " << max_cell_in_box << ")" << endl;
	this->end();
	exit(1);
      }

      //id of the new cell in the box 
      unsigned int new_cell_id = this->boxes_A[c1][c2][c3].cells.size();
      
      // create a new cell
      Cell newcell;
      boxes_A[c1][c2][c3].cells.push_back(newcell);
      
      // ============================
      // set the properties of the new cell
      this->boxes_A[c1][c2][c3].cells[new_cell_id].birthday = this->reloj;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].mother_name = cell.name;
      // fill a new element of the array box[][][].cells
      this->boxes_A[c1][c2][c3].cells[new_cell_id].position[0] = newpositionx;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].position[1] = newpositiony;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].position[2] = newpositionz;

      // old position: position of mother cell
      this->boxes_A[c1][c2][c3].cells[new_cell_id].position_old[0] = cell.position[0];
      this->boxes_A[c1][c2][c3].cells[new_cell_id].position_old[1] = cell.position[1];
      this->boxes_A[c1][c2][c3].cells[new_cell_id].position_old[2] = cell.position[2];

      this->boxes_A[c1][c2][c3].cells[new_cell_id].box[0]=c1;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].box[1]=c2;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].box[2]=c3;
      
      this->boxes_A[c1][c2][c3].cells[new_cell_id].name=this->total_no_of_cells+cells_counter;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].type=cell.type;

      // cell velocity (polarity) CHECK THIS
      this->boxes_A[c1][c2][c3].cells[new_cell_id].polarity[0] = cell.polarity[0];
      this->boxes_A[c1][c2][c3].cells[new_cell_id].polarity[1] = cell.polarity[1];
      this->boxes_A[c1][c2][c3].cells[new_cell_id].polarity[2] = cell.polarity[2];
      
      // oxygen concentration
      this->boxes_A[c1][c2][c3].cells[new_cell_id].O2 = 100.;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].dxO2 = 0.;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].dyO2 = 0.;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].dzO2 = 0.;

      // adhesion constant
      this->boxes_A[c1][c2][c3].cells[new_cell_id].adhesion = cell.adhesion;  

      this->boxes_A[c1][c2][c3].cells[new_cell_id].radius = new_radius;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].clear_contacts();
      this->boxes_A[c1][c2][c3].cells[new_cell_id].type = 1;
      this->boxes_A[c1][c2][c3].cells[new_cell_id].hypoxic_count = 0;
      // initialize counter of hypoxic state
      this->boxes_A[c1][c2][c3].cells[new_cell_id].phenotype_counter = 0;

      // follower-leader behavior
      this->boxes_A[c1][c2][c3].cells[new_cell_id].follow_lead = cell.follow_lead;
      // determine oxygen phenotype of new cell
      double alea = aleatorio();
      float random_threshold=0.5;
      if(alea>random_threshold) {

	double valor=params.threshold_hypo+params.variance_phenotype * (1.-2.*rand())/(RAND_MAX+0.0);

	if (valor>params.threshold_death) {
	  this->boxes_A[c1][c2][c3].cells[new_cell_id].phenotype=valor;
	} else {
	  this->boxes_A[c1][c2][c3].cells[new_cell_id].phenotype=cell.phenotype;
	}
      } else {
	this->boxes_A[c1][c2][c3].cells[new_cell_id].phenotype=cell.phenotype;
      }

      // determine interaction phenotype of the new cell
      this->boxes_A[c1][c2][c3].cells[new_cell_id].interaction_phenotype = cell.interaction_phenotype;

      // determine polarisation of new cell
      this->boxes_A[c1][c2][c3].cells[new_cell_id].polarised = 0;


      // if (total_no_of_cells>100){
      // 	if (cell.interaction_phenotype==1){
      // 	  double prob= aleatorio();
      // 	  if (prob>0.5) {
      // 	    this->boxes_A[c1][c2][c3].cells[new_cell_id].interaction_phenotype = 0;
      // 	  }
      // 	}
      // }
      // ============================


      
      if (this->boxes_A[c1][c2][c3].cells[new_cell_id].interaction_phenotype==0)
      {
	this->phenotype1_count=this->phenotype1_count+1;
      }
      else
      {
	this->phenotype2_count=this->phenotype2_count+1;
       }

      // counts the number of newborn cells
      cells_counter++;
      
      
      // updating the box domain
      if(this->maxx<c1) {
	this->new_maxx=c1;

	if (c1==(int) params.boxesx){
	  cout << " ** error: Cell " << this->boxes_A[c1][c2][c3].cells[new_cell_id].name
	       << " is moving out of the domain (not supported) " << endl;
	  exit(1);
	  //this->new_maxx=0;	
	  //c1=0;
	}

      }
      if(this->maxy<c2) {
	this->new_maxy=c2;	

	if (c2==(int) params.boxesy) {
	  cout << " ** error: Cell " << this->boxes_A[c1][c2][c3].cells[new_cell_id].name
	       << " is moving out of the domain (not supported) " << endl;
	  exit(1);
	  //this->new_maxy=0;
	  //c2=0;
	}
      }
      
      if(this->maxz<c3) {
	this->new_maxz=c3;
	
	if (c3==(int) params.boxesz) { 
	  cout << " ** error: Cell " << this->boxes_A[c1][c2][c3].cells[new_cell_id].name
	       << " is moving out of the domain (not supported) " << endl;
	  exit(1);
	  //this->new_maxz=0;	
	  //c3=0;
	}
      }	
      if(this->minx>c1) {
	this->new_minx=c1;

	if (c1<0) {
	  cout << " ** error: Cell " << this->boxes_A[c1][c2][c3].cells[new_cell_id].name
	       << " is moving out of the domain (not supported) " << endl;
	  exit(1);
	  //this->new_minx=params.boxesx;	
	  //c1=params.boxesx-1;
	}
      }
      
      if(this->miny>c2) {
	this->new_miny=c2;

	if (c2<0) {
	  cout << " ** error: Cell " << this->boxes_A[c1][c2][c3].cells[new_cell_id].name
	       << " is moving out of the domain (not supported) " << endl;
	  exit(1);
	  //this->new_miny=params.boxesy;
	  //c2=params.boxesy-1;
	}
      }
      
      if(this->minz>c3) {
	this->new_minz=c3;	

	if (c3<0) {
	  cout << " ** error: Cell " << this->boxes_A[c1][c2][c3].cells[new_cell_id].name
	       << " is moving out of the domain (not supported) " << endl;
	  exit(1);
	  //this->new_minz=params.boxesz;	
	  //c3=params.boxesz-1;
	}
      }
    }
  } // if birth_allowed

  //grow if still possible
  if(cell.radius<max_radius_cell) {
    cell.radius += params.growth_rate[cell.interaction_phenotype]*params.time_step;
  }	

}

/*****************************************************************************/
/* velocity correction due to fibres                                         */
/*****************************************************************************/
void CoupledModel::cell_fibres_interaction(Cell& cell) {

  vector<double> fibre_adhesion, fibre_contact;
  fibre_adhesion.resize(params.dimension);
  fibre_adhesion.clear();
  fibre_contact.resize(params.dimension);
  fibre_contact.clear();
  unsigned int n_fibres = cell.contact_fibres.size();
  double velocity_dot_direction = 0.;
  double cell_velocity = 0.;
  for(unsigned int k=0; k<n_fibres; k++) {
    
    // compare fibre orientation and previous cell velocity (polarity)
    for (unsigned int j=0; j<params.dimension; j++) {
      velocity_dot_direction += cell.contact_fibres[k]->direction[j]*cell.polarity[j];
      //cout << velocity_dot_direction << endl;
      cell_velocity += cell.polarity[j]*cell.polarity[j];
      //cout << cell_velocity << endl;
    }

    cell_velocity = sqrt(cell_velocity);
    
    // average fibre direction
    for (unsigned int j=0; j<params.dimension; j++) {
      // vel_adhesion=1e-3
      fibre_adhesion[j] += params.vel_adhesion*velocity_dot_direction*
	cell.contact_fibres[k]->direction[j];
      double x = velocity_dot_direction/(cell_velocity+1e-8);
      fibre_contact[j] += -params.vel_contact*cell.polarity[j]*(1-x*x);
    }
  }


  // compute new velocity
  for (unsigned int j=0; j<params.dimension; j++) {
    
     cell.vel[j] += fibre_adhesion[j] + fibre_contact[j] + cell.polarity[j];
   }
  /*cout << "cell polarity: " << cell.polarity[0] << " " << cell.polarity[1] << " " << cell.polarity[2] << endl;
  cout << cell_velocity << " " << cell.vel[0] << " " << cell.vel[1] << " " << cell.vel[2] <<endl;
  cout << " n_fibres = " << n_fibres << " adhesion force: "
       << sqrt(fibre_adhesion[0]*fibre_adhesion[0]+fibre_adhesion[1]*fibre_adhesion[1]+fibre_adhesion[2]*fibre_adhesion[2])
       << " contact force: " 
       << sqrt(fibre_contact[0]*fibre_contact[0]+fibre_contact[1]*fibre_contact[1]+fibre_contact[2]*fibre_contact[2])
       << " velocity: "
       << sqrt(cell.vel[0]*cell.vel[0] + cell.vel[1]*cell.vel[1] + cell.vel[2]*cell.vel[2]) << endl;*/

}
/*****************************************************************************/
/* hertz: gives the force from the hertz model                               */
/*****************************************************************************/
void CoupledModel::hertz(Cell& cell) {

  ///@todo this should be done somewhere else
  oxy_in_cell(cell);
  
  if (params.verbose > 3) {
    if (cell.neighbors.size()) {
      cout << " compute contact forces for cell " << cell.name << endl;
      cout << " contact with: ";
      for (unsigned int j=0; j < cell.neighbors.size(); j++) {
	cout << cell.neighbors[j]->name << " ";
      }
      cout << endl;
    }
  }

  cell.energy = 0;

  // ******* vessel-cell forces ************************************

  // attractive-repulsive forces
  double  c_v_fx=0.;
  double  c_v_fy=0.;
  double  c_v_fz=0.;

  for(unsigned int kk=0; kk<cell.contact_vessels.size(); kk++) {

    double cell_pois = params.alpha_PoissonNo[cell.interaction_phenotype];
    double cell_K =  (1.- cell_pois*cell_pois)/params.alpha_YoungM[cell.interaction_phenotype];
    double ves_pois = params.vessel_PoissonNo;
    double ves_K =  (1.- ves_pois*ves_pois)/params.vessel_YoungM;
    double eff_mod = 1.0/(cell_K = ves_K);
    //cout << eff_mod << endl;

    double c_v_disp = 0.0;
    for (unsigned int j=0; j < vessels.size(); j++) {

      if (vessels[j].vessel_name == cell.contact_vessels[kk]->vessel_name) {

	double cell_vessel_dist = DISTANCE(cell,this->vessels[j]);
	// displacement
	c_v_disp = cell.radius + cell.contact_vessels[kk]->ves_radius - cell_vessel_dist;

	this->VESSELPOINT(cell,this->vessels[j]);	
      }
    }
    //cout << c_v_disp << endl;

    double c_v_adh_coeff = fmin(cell.adhesion,params.vessel_adhesion);
    double c_v_f_ij =
      //c_v_repulsion
      4./3. * eff_mod * sqrt(cell.radius) * pow(c_v_disp,1.5)
      //c_v_adhesion
      -c_v_adh_coeff * (cell.radius*c_v_disp - c_v_disp*c_v_disp/4.);
    // NOT SURE ABOUT THIS BIT 
    double surf1 = 2 * PIG * c_v_disp;
    cell.energy += fabs(c_v_f_ij)/surf1;

    // co-ordinate distances between cell centre and vessel
    double c_v_x = this->vessel_point[0] - cell.position[0];
    double c_v_y = this->vessel_point[1] - cell.position[1];
    double c_v_z = this->vessel_point[2] - cell.position[2];
    double dist_cell_vessel = sqrt(c_v_x*c_v_x + c_v_y*c_v_y + c_v_z*c_v_z);
    //cout << dist_cell_vessel << endl;
    
    c_v_x = c_v_x/dist_cell_vessel;
    c_v_y = c_v_y/dist_cell_vessel;
    c_v_z = c_v_z/dist_cell_vessel;

    c_v_fx = c_v_fx + c_v_f_ij * (-c_v_x);
    c_v_fy = c_v_fy + c_v_f_ij * (-c_v_y);
    c_v_fz = c_v_fz + c_v_f_ij * (-c_v_z);

  }
  
  // attractive-repulsive forces
  //double  fx=0.;
  //double  fy=0.;
  //double  fz=0.;
  double  fx=c_v_fx;
  double  fy=c_v_fy;
  double  fz=c_v_fz;
  //cell.energy=0;

  // forces for cell migration (leader-follower)
  double migration_force_x=0.;
  double migration_force_y=0.;
  double migration_force_z=0.;
  
  
  for(unsigned int k=0; k<cell.contacts; k++) {

    // original version (Caiazzo & Ramis, JTB 2015)
    //double K_tap= 3./4.* 2 * (1.-params.PoissonNo*params.PoissonNo) / params.YoungM;
    // = 3./4. (1.-params.PoissonNo*params.PoissonNo) / (params.YoungM)
    //+ (1.-params.PoissonNo*params.PoissonNo) / (params.YoungM)
    
    // K = (1-nu^2)/E (for considered cell and neighbor)
    double nu = params.alpha_PoissonNo[cell.interaction_phenotype];
    double cell_K =  (1.- nu*nu)/params.alpha_YoungM[cell.interaction_phenotype];
    
    nu = params.alpha_PoissonNo[cell.neighbors[k]->interaction_phenotype];
    double neighbour_K = (1.-nu*nu)/params.alpha_YoungM[cell.neighbors[k]->interaction_phenotype];

    // sum of Ks
    double eff_K = (cell_K + neighbour_K);

    // effective radius
    double eff_radius = cell.radius * cell.neighbors[k]->radius /
      (cell.radius + cell.neighbors[k]->radius);

    double adhesion_coeff = fmin(cell.adhesion,cell.neighbors[k]->adhesion);
     // ******************************
    
    // distance between centers
    double p_x = cell.neighbors[k]->position[0] - cell.position[0];
    double p_y = cell.neighbors[k]->position[1] - cell.position[1];
    double p_z = cell.neighbors[k]->position[2] - cell.position[2];
    double dist_cells = sqrt(p_x*p_x + p_y*p_y + p_z*p_z);
    
    double e_x = p_x/dist_cells;
    double e_y = p_y/dist_cells;
    double e_z = p_z/dist_cells;
   
    // distance between surfaces
    double d_ij = cell.radius + cell.neighbors[k]->radius - dist_cells;

    /**
       Theoretically d_ij>0 when cells are in contact. However, in some cases,
       (e.g. follow-leader) we also track the neighbors at time n-1 that
       are no longer neighbors at time n. In this case, we have to set manually
       d_ij = 0 so that the following terms vanished
    **/
    if (d_ij<1e-8) {
      d_ij = 0;
    }
    
    // compute force
    double f_ij =
      // repulsion
      4./3. * (1./eff_K) * sqrt(eff_radius) * pow(d_ij,1.5) // NEW
      // adhesion
      //-cell.adhesion * (cell.radius*d_ij - d_ij*d_ij/4.);
      -adhesion_coeff * (cell.radius*d_ij - d_ij*d_ij/4.);
      //-adhesion_coeff * surface_area_new; // IGNACIO WONDERED ABOUT CHANGING THE ADHESION FORCE
    double surface = 2 * PIG * d_ij;
    cell.energy += fabs(f_ij)/surface;

    // // =================================================================
    // // write cell energies to data file
    // // =================================================================
    // ofstream cellenergy;
    // string ce_list =  "cell_energies.txt";
    // cellenergy.open(ce_list.c_str(),ios::app);
    // cellenergy << cell.energy << " " ;
    // cellenergy.close();
    // // ********************************

    //Projection of the force over the three dimensions
    fx = fx + f_ij * (-e_x);
    fy = fy + f_ij * (-e_y);
    fz = fz + f_ij * (-e_z);

    
    // migration: based on variation in contact area
    ///@todo need to handle the case when the cells where not in contact before?
    double p_x_old = cell.neighbors[k]->position_old[0] - cell.position_old[0];
    double p_y_old = cell.neighbors[k]->position_old[1] - cell.position_old[1];
    double p_z_old = cell.neighbors[k]->position_old[2] - cell.position_old[2];
    double dist_cells_old = sqrt(p_x_old*p_x_old + p_y_old*p_y_old + p_z_old*p_z_old);
    double d_ij_old = cell.radius + cell.neighbors[k]->radius - dist_cells_old;
    
    //double contact_area_new = 2 * PIG * cell.radius * (d_ij)/2.0;
    double contact_area_new = 2 * PIG * this->Ro * (d_ij)/2.0;
    double contact_area_old = 2 * PIG * this->Ro * (d_ij_old)/2.0;
    
    double variation_area = contact_area_new-contact_area_old;

    /*cout << "     cell " << cell.name << " in contact with cell " << cell.neighbors[k]->name << endl;
    cout << "          cell " << cell.name << " is in position " << cell. position[0] << endl;
    cout << "          cell " << cell.neighbors[k]->name << " is in position "
	 << cell.neighbors[k]->position[0] << endl;
    cout << "          current contact area is " << contact_area_new << endl;
    cout << "          previous contact area is " << contact_area_old << endl;
    cout << "          variation in contact areas is " << variation_area << endl;*/

    /* The additional force is only applied only if:
       - the cells are greater than 7 units apart 
       (otherwise cells envelop other cells)
       - variation_area < 0
    */
    
    if ((variation_area<0)&&(dist_cells>7)) {

      // intensity = f* a^2/(1+a^2)
      double add_denominator = 0.1;
      double intensity_migration_force = params.follower_force*
	variation_area*variation_area/(add_denominator+variation_area*variation_area);
      
      migration_force_x = migration_force_x + intensity_migration_force*e_x;
      migration_force_y = migration_force_y + intensity_migration_force*e_y;
      migration_force_z = migration_force_z + intensity_migration_force*e_z;
    }
    
    // Alfonso: this normalization looks strange here:
    // normalise_force = variation_area
    // I moved it out of the loop
    /*double normalise_force = sqrt(migration_force_x*migration_force_x +
      migration_force_y*migration_force_y +
      migration_force_z*migration_force_z);
      if (normalise_force!=0){
      migration_force_x = migration_force_x/normalise_force;
      migration_force_y = migration_force_y/normalise_force;
      migration_force_z = migration_force_z/normalise_force;
      }*/
    
  }

  /*
  // normalize migration forces
  double normalise_force = sqrt(migration_force_x*migration_force_x +
				migration_force_y*migration_force_y +
				migration_force_z*migration_force_z);
  if (normalise_force!=0){
    migration_force_x = migration_force_x/normalise_force;
    migration_force_y = migration_force_y/normalise_force;
    migration_force_z = migration_force_z/normalise_force;
  }
  */
  
  //cout << " m.force: " << migration_force_x << " " << migration_force_y << endl;
  
  // compute cell velocity
  switch (cell.type) {
  case 3:
    // dead cells
    cell.vel[0] = fx/(params.Gcm);
    cell.vel[1] = fy/(params.Gcm);
    cell.vel[2] = fz/(params.Gcm);
    break;
    
  case 2:
    // hypoxic cells
    // - might respond to oxygem gradien 
    // - higher diffusion coefficient (10x)
    if (cell.O2<cell.phenotype ) {
      double normgrad = sqrt(cell.dxO2*cell.dxO2 + cell.dyO2*cell.dyO2 + cell.dzO2*cell.dzO2)+1;
      // renormalized chemotaxis term: psi* gradc/(1+psi*|gradc|)  
      cell.vel[0] += (fx+ box_muller(0,vf*10))/(params.Gcm) + params.oxygen_response * 
	cell.dxO2/(1 + normgrad*params.oxygen_response)/params.Gcm;
      cell.vel[1] += (fy+ box_muller(0,vf*10))/(params.Gcm) + params.oxygen_response * 
	cell.dxO2/(1 + normgrad*params.oxygen_response)/params.Gcm;
      cell.vel[2] += (fz+ box_muller(0,vf*10))/(params.Gcm) + params.oxygen_response * 
      cell.dxO2/(1 + normgrad*params.oxygen_response)/params.Gcm;	      
    } else {
      // @attention this looks strange?
      // cell behaves as normoxic with enough oxygen   
      cell.vel[0] += (fx+ box_muller(0,vf))/(params.Gcm);
      cell.vel[1] += (fy+ box_muller(0,vf))/(params.Gcm);
      cell.vel[2] += (fz+ box_muller(0,vf))/(params.Gcm);
    }
    break;
    
  case 1: {
    // normoxic    
    double friction = params.alpha_gcm[cell.interaction_phenotype]*
      (params.Gcm);
    //double friction = params.alpha_gcm[cell.interaction_phenotype]*
    //(params.Gcm + 0.05*(1.+tanh(cell.position[0]-300)));

    // FOLLOWER CELLS
    /* for follower cells an additional force is included */
    //double f_multiplier = params.follower_force;
    //double add_denominator = 1.0;
    if(cell.follow_lead==1){
      cell.vel[0] = (fx + box_muller(0,vf) + migration_force_x)/friction;
      cell.vel[1] = (fy + box_muller(0,vf) + migration_force_y)/friction;
      cell.vel[2] = (fz + box_muller(0,vf) + migration_force_z)/friction;
      /*cell.vel[0] = (fx + box_muller(0,vf)+
		     (f_multiplier*migration_force_x*migration_force_x)/
		     (add_denominator+migration_force_x*migration_force_x))/friction;
      cell.vel[1] = (fy + box_muller(0,vf)+
		     (f_multiplier*migration_force_y*migration_force_y)/
		     (add_denominator+migration_force_y*migration_force_y))/friction;
      cell.vel[2] = (fz + box_muller(0,vf)+
		     (f_multiplier*migration_force_z*migration_force_z)/
		     (add_denominator+migration_force_z*migration_force_z))/friction;*/
      /*cout << " migration force: " << f_multiplier*migration_force_x*migration_force_x/
      (add_denominator+migration_force_x*migration_force_x) << " " <<
	f_multiplier*migration_force_y*migration_force_y/
	(add_denominator+migration_force_y*migration_force_y) << endl;*/
    }
    // LEADER CELLS
    else {
      cell.vel[0] = (fx + box_muller(0,vf))/friction;
      cell.vel[1] = (fy + box_muller(0,vf))/friction;
      cell.vel[2] = (fz + box_muller(0,vf))/friction;
      }
    
    break;
  }
  default:
    cout << " *** ERROR in file " << __FILE__ << ", line " << __LINE__
	 << " cell " << cell.name << " has type: " << cell.type  << ", unknown. " << endl;
    exit(1);
    break;

  }

  // handle the 2-dimensional case
  if (params.dimension == 2)  cell.vel[2] = 0.;
				
  //  cout << " dx[ " << cell.name << "] = " << this->dx[cell.name] << ","
  //   << this->dy[cell.name] << " " << this->dz[cell.name] << endl;
  //cout << " force: " << fx << "," << fy << "," << fz << endl;
  //cout << " velocity: " << cell.vel[0] << "," << cell.vel[1] << "," << cell.vel[2] << endl;
  
}



/*************************************************************************
 * moves the cells acording to the velocity previously computed                     
*************************************************************************/
void CoupledModel::movement(const Cell& cell, 
	      const int u, const int v, const int w, 
	      const unsigned int cont_cell)
{
  // copy the cell
  Cell celula_nueva = cell;

  // save the old position
  for (unsigned int j=0; j<3; j++) {
    celula_nueva.position_old[j] = cell.position[j];
  }
  double A = 0.0; 
  //double omega = 1.0/10.0;
  //if(this->total_no_of_cells==100.0){
  // WE CAN SWITCH ON THE TRAJECTORY OF THE LEADER CELLS ONCE WE HAVE GROWN A POPULATION OF CELLS
  A = 5.0;
  //}
      
  // compute new position
  //int idcell = this->boxes_A[u][v][w].cells[cont_cell].name;
  //A*sin(omega*reloj)*params.time_step // AN ALTERNATIVE TRAJECTORY
  //LEADER CELL
  /* the force-based movement of cells now includes an additional term for a basic linear trajectory of leader cells. */
  if (cell.follow_lead==0){
    celula_nueva.position[0]= this->boxes_A[u][v][w].cells[cont_cell].position[0]  + A*params.time_step +
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[0];
    celula_nueva.position[1]= this->boxes_A[u][v][w].cells[cont_cell].position[1] + A*params.time_step + 
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[1];
    celula_nueva.position[2]= this->boxes_A[u][v][w].cells[cont_cell].position[2] + 
      movez * params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[2];
  }
  // FOLLOWER CELLS
  else {
    celula_nueva.position[0]= this->boxes_A[u][v][w].cells[cont_cell].position[0] + 
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[0];
    celula_nueva.position[1]= this->boxes_A[u][v][w].cells[cont_cell].position[1] + 
      params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[1];
    celula_nueva.position[2]= this->boxes_A[u][v][w].cells[cont_cell].position[2] + 
      movez * params.time_step * this->boxes_A[u][v][w].cells[cont_cell].vel[2];
  }
      
  // periodic bc
  if (celula_nueva.position[0]<0) {
    celula_nueva.position[0]=params.lattice_length_x+celula_nueva.position[0];
  } 
  if (celula_nueva.position[0]>params.lattice_length_x) {
    celula_nueva.position[0]=celula_nueva.position[0]-params.lattice_length_x;
  }
  
  if (celula_nueva.position[1]<0) {
    celula_nueva.position[1]=params.lattice_length_y+celula_nueva.position[1];
  } 
  if (celula_nueva.position[1]>params.lattice_length_y) {
    celula_nueva.position[1]=celula_nueva.position[1]-params.lattice_length_y;
  }
  
  
  if (celula_nueva.position[2]<0) {
    celula_nueva.position[2]=params.lattice_length_z+celula_nueva.position[2];
  } 
  if (celula_nueva.position[2]>params.lattice_length_z) {
    celula_nueva.position[2]=celula_nueva.position[2]-params.lattice_length_z;
  }
  
  //cout << " set new box " << endl;
  
  //nueva caja
  int c1=(int)(floor(celula_nueva.position[0]/this->box_sizex)); 
  int c2=(int)(floor(celula_nueva.position[1]/this->box_sizey)); 
  int c3=(int)(floor(celula_nueva.position[2]/this->box_sizez)); 
  
  celula_nueva.box[0]=c1;
  celula_nueva.box[1]=c2;
  celula_nueva.box[2]=c3;

  if( ((unsigned int) c1 >= params.boxesx) || (c1 < 0) ||
      ((unsigned int) c2 >= params.boxesy) || (c2 < 0) ||
      ((unsigned int) c3 >= params.boxesz) || (c3 < 0) ) {
    
    cout << " ** error: the cell " << cell.name << " is moving out of the domain (not supported) "
	 << endl;
    exit(1);
    }


  //cout << " replace the cell in the box " << endl;
  /// @todo this should be done without copying the cell
  this->boxes_new_A[c1][c2][c3].cells.push_back(celula_nueva);
  
  // update maximum and minumum boxes if necessary
  if(this->maxx<c1) 
    this->new_maxx=c1;
  
  if(this->maxy<c2) 
    this->new_maxy=c2;	
  
  if(this->maxz<c3) 
    this->new_maxz=c3;	
  
  if(this->minx>c1) 
    this->new_minx=c1;
  
  if(this->miny>c2) 
    this->new_miny=c2;
  
  if(this->minz>c3) 
    this->new_minz=c3;	
  	

  
}

// FIBRE INDUCED MOVEMENT - CICELY NEW
// void CoupledModel::fibre_induced_movement(Cell& cell)
// {
//   //double pol_mag = sqrt(pow(cell.pol_axis[0],2)+pow(cell.pol_axis[1],2)+pow(cell.pol_axis[2],2));   // TODO include pol_axis[3] as a cell variable in Cell.h
//   double pol_axis[3];
//   double prob[fibres.size()];
//   double rand_num = aleatorio(); /// TODO check this generates the right type of random number
//   double sum = 0;
//   double size = 0;
//   double mtmx, mtmy, mtmz;
  
//   for (unsigned int j=0; j<fibres.size(); j++)
//   {
//     for (unsigned int i=0; i<=2; i++)
//     {
//       pol_axis[i]=1.0;
// 	//cell.pol_axis[i]/pol_mag;   // TODO include pol_axis[3] as a cell variable in Cell.h
//       prob[j] = prob[j] + pol_axis[i]*fibres[j].direction[i];
//       size = size + prob[j]; 
//     }
    
//     for (unsigned int i=0; i<=2; i++)
//     {
//       prob[j] = prob[j]/size;
//     }

//     if (rand_num > sum && rand_num <= sum+prob[j])
//     {
//       mtmx = fibres[j].direction[0];
//       mtmy = fibres[j].direction[1];
//       mtmz = fibres[j].direction[2];
//     }
//     else
//     {
//       sum = sum + prob[j];
//     }
//   }

//   // double num_cells = this->total_no_of_cells;
//   // double fibspeed[numcells];
//   // for (unsigned int k=0; k<num_cells; k++)
//   // {
//   //   fibspeed[k] = sqrt(pow(mtmx,2)+pow(mtmy,2)+pow(mtmz,2));
//   //   mtmx = (mtmx/fibspeed[k])*(max_cell_speed - (max_cell_speed/2500)*pow(integrins[k],2));
//   //   mtmy = (mtmy/fibspeed[k])*(max_cell_speed - (max_cell_speed/2500)*pow(integrins[k],2));
//   //   mtmz = (mtmz/fibspeed[k])*(max_cell_speed - (max_cell_speed/2500)*pow(integrins[k],2));
//   // } 
// }



/**********************************************************************************************/
/* calculate the local distance between cells and fibres - CICELY NEW */ 
/**********************************************************************************************/
void CoupledModel::compute_cell_fibres_contact(const int u,
		 const int v,const int w) 
{
  unsigned int total_cells_in_box=this->boxes_A[u][v][w].cells.size();
  for(unsigned int i=0; i<total_cells_in_box; i++){
      //cout << " compute_cell_fibres_contact process cell: " << boxes_A[u][v][w].cells[i].name << endl;
      // find neighbouring fibres
      /// @brief number of fibres in the neighbouring box
      //unsigned int fibres_in_box=this->boxes_A[u+h][v+l][w+m].b_fibres.size();
      // cout << fibres_in_box << endl;
	    
      for(unsigned int j=0; j<fibres.size(); j++){
	  //cout << " process fibre: " << fibres[j].fibre_name << endl;


	  if (fibres[j].fibre_exists>0.5){
	    double fibre_centre[3];
	    double c_f_centre=0;
	    for (unsigned int jj=0; jj<=2; jj++){
	      fibre_centre[jj] = fibres[j].start[jj]+0.5*fibres[j].length*fibres[j].direction[jj];
	      c_f_centre = c_f_centre+pow(boxes_A[u][v][w].cells[i].position[jj]-fibre_centre[jj],2);
	    }
      
	    double cell_fibre_dist_check = sqrt(c_f_centre);
	    //cout << cell_fibre_dist_check << endl;
	    double dist_ref = boxes_A[u][v][w].cells[i].radius + fibres[j].length;
	    if(cell_fibre_dist_check < dist_ref)
	      {
		// compute distance between cells and fibres
		double cell_fibre_min_dist = DISTANCE(this->boxes_A[u][v][w].cells[i],this->fibres[j]);
		double cell_rad = this->boxes_A[u][v][w].cells[i].radius;

		
		if(cell_fibre_min_dist < cell_rad*2.0)
		  {
		    fibre_degradation(this->boxes_A[u][v][w].cells[i],fibres[j],cell_fibre_min_dist);
		    /*double c_to_f_start = 0.0;
		    double c_to_f_end = 0.0;
		    for (unsigned int l=0; l<3; l++) {
		      c_to_f_start += (this->boxes_A[u][v][w].cells[i].position[l]-fibres[j].start[l])*this->boxes_A[u][v][w].cells[i].polarity[l];
		      c_to_f_end += (this->boxes_A[u][v][w].cells[i].position[l]-(fibres[j].start[l]+fibres[j].length*fibres[j].direction[l]))*this->boxes_A[u][v][w].cells[i].polarity[l];
		    }
		    
		    if (c_to_f_start>0.0 || c_to_f_end>0.0){
		      
		      if(cell_fibre_min_dist < cell_rad)
			{
			  double velocity_dot_direction = 0.0;
			  double cell_velocity = 0.0;
			  for (unsigned int k=0; k<3; k++) {
			    velocity_dot_direction += fibres[j].direction[k]*this->boxes_A[u][v][w].cells[i].polarity[k];
			    cell_velocity += this->boxes_A[u][v][w].cells[i].polarity[k]*this->boxes_A[u][v][w].cells[i].polarity[k];
			  }
			  cell_velocity = sqrt(cell_velocity);
			  velocity_dot_direction = velocity_dot_direction/cell_velocity;
			  if (fabs(velocity_dot_direction)<sqrt(2.0)/2.0){
			    double rand_contact_degrade = aleatorio();
			    if (rand_contact_degrade<0.5){
			      fibres[j].fibre_exists=0.0;
			    }
			  }
			}
		    }
		    double rand_diffuse_degrade = aleatorio();
		    if (rand_diffuse_degrade<0.05){
		      fibres[j].fibre_exists=0.0;
		    }
		    */
		    if (params.verbose>1) {
		      cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " in contact with fibre " << this->fibres[j].fibre_name << endl;
		      cout << " distance between the cell and fibre is " << cell_fibre_min_dist << endl;
		    }
		    this->boxes_A[u][v][w].cells[i].contact_fibres.push_back(&fibres[j]);
		    
		  }
	      }
	  }
      }
  }
}

/**********************************************************************************************/
/* calculate the local distance between cells and vessels - CICELY NEW 06/06/17 */ 
/**********************************************************************************************/
void CoupledModel::compute_cell_vessels_contact(const int u,
		 const int v,const int w) 
{
  unsigned int total_cells_in_box = this->boxes_A[u][v][w].cells.size();
  for(unsigned int i=0; i<total_cells_in_box; i++){

      this->boxes_A[u][v][w].cells[i].interaction_phenotype = 0.0;
	    
      for(unsigned int j=0; j<vessels.size(); j++){

	    double vessel_centre[3];
	    double c_to_v_centre=0;
	    for (unsigned int jj=0; jj<=2; jj++){
	      vessel_centre[jj] = vessels[j].ves_start[jj]+0.5*vessels[j].ves_length*vessels[j].ves_direction[jj];
	      c_to_v_centre = c_to_v_centre+pow(boxes_A[u][v][w].cells[i].position[jj]-vessel_centre[jj],2);
	    }
      
	  double cell_vessel_dist_check = sqrt(c_to_v_centre);
	  double dist_ref = boxes_A[u][v][w].cells[i].radius + vessels[j].ves_length;
	  if(cell_vessel_dist_check < dist_ref){
	      // compute distance between cells and vessels
	      double cell_vessel_min_dist = DISTANCE(this->boxes_A[u][v][w].cells[i],this->vessels[j]);
	      double cell_rad = this->boxes_A[u][v][w].cells[i].radius;
	      double vessel_rad = vessels[j].ves_radius;
	
	      if(cell_vessel_min_dist < (cell_rad+vessel_rad)){
		  if (params.verbose>1) {
		    cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " in contact with vessel " << this->vessels[j].vessel_name << endl;
		    cout << " distance between the cell and vessel is " << cell_vessel_min_dist-(cell_rad+vessel_rad) << endl;
		  }
		  this->boxes_A[u][v][w].cells[i].contact_vessels.push_back(&this->vessels[j]);

		  this->boxes_A[u][v][w].cells[i].interaction_phenotype = 1.0;
	      }
	  }
      }
  }
}

/**********************************************************************************************/
/* calculate the local distance between cells in the lattice and produces the conexion matrices*/
/**********************************************************************************************/
void CoupledModel::contact_forces(const int u,
		 const int v,const int w,
		 const unsigned int cells_number) 
{
  int borderx,bordery,borderz;

  
  for(unsigned int i=0; i<cells_number; i++) {

    // loop in neighbor boxes to find cells in contact with cell i
    // we consider a layer of [-1,1]x[-1,1]x[-1,1] 
    for(int h=-1;h<2;h++) {
      borderx=u+h;
      for(int l=-1;l<2;l++) {
	bordery=v+l;
	for(int m=-1;m<2;m++) {
	  borderz=w+m;		   		   
	  
	  // check that we are looking inside the computational domain
	  if( borderx<=(int) params.boxesx-1 && borderx>=0 && 
	      bordery<=(int) params.boxesy-1 && bordery>=0 && 
	      borderz<=(int) params.boxesz-1 && borderz>=0 ) {

	    for(unsigned int j=0; j<this->boxes_A[u+h][v+l][w+m].cells.size(); j++) {
	      // check that we are not looking at the same cell
	      if(this->boxes_A[u][v][w].cells[i].name != this->boxes_A[u+h][v+l][w+m].cells[j].name) {
		
		// compute distance between cells
		double cell_cell_min = DISTANCE(this->boxes_A[u][v][w].cells[i],
					  this->boxes_A[u+h][v+l][w+m].cells[j]); 
		double cell_cell_center = this->boxes_A[u+h][v+l][w+m].cells[j].radius +
		  this->boxes_A[u][v][w].cells[i].radius;
		
		// cells in contact
		if(cell_cell_min < cell_cell_center){

		  // store the pointer to the neighbor cell in the cell.neighbors vector
		  this->boxes_A[u][v][w].cells[i].neighbors.push_back(&this->boxes_A[u+h][v+l][w+m].cells[j]);   
		  // increase number of contacts
		  this->boxes_A[u][v][w].cells[i].contacts++;

		  if (params.verbose>3) {
		    cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " in  contact with "
			 << this->boxes_A[u+h][v+l][w+m].cells[j].name << endl;
		  }

		} else {

		  if (params.follower_force>0) {

		    if ( (this->boxes_A[u][v][w].cells[i].birthday < reloj)&&
			 (this->boxes_A[u+h][v+l][w+m].cells[j].birthday < reloj)
			 ) {
		      // check if cells were in contact before
		      // compute distance between cells
		      double cell_cell_dist_old = 
			pow(this->boxes_A[u][v][w].cells[i].position_old[0] -
			    this->boxes_A[u+h][v+l][w+m].cells[j].position_old[0],2)+
			pow(this->boxes_A[u][v][w].cells[i].position_old[1] -
			    this->boxes_A[u+h][v+l][w+m].cells[j].position_old[1],2)+
			pow(this->boxes_A[u][v][w].cells[i].position_old[2] -
			    this->boxes_A[u+h][v+l][w+m].cells[j].position_old[2],2);
		      cell_cell_dist_old = sqrt(cell_cell_dist_old);
		      
		      double cell_cell_center = this->boxes_A[u+h][v+l][w+m].cells[j].radius +
			this->boxes_A[u][v][w].cells[i].radius;
		      
		      if( cell_cell_dist_old < cell_cell_center){
			cout << " time step: " << this->reloj
			     << ", cells " << this->boxes_A[u][v][w].cells[i].name
			     << " and " << this->boxes_A[u+h][v+l][w+m].cells[j].name
			     << " were in contact at the previous time step " << endl;

			this->boxes_A[u][v][w].cells[i].neighbors.
			  push_back(&this->boxes_A[u+h][v+l][w+m].cells[j]);   
			// increase number of contacts
			this->boxes_A[u][v][w].cells[i].contacts++;
		      }
		    }
		  }
		  
		}//end if
		
	      }
	    } // for(unsigned int j=0;j<neighbors box cells;j++)
	    
	  }//end if (inside domain)
	  
	}  
      }
    }//end for loops
    
    // for debugging: print all contacts
    if (params.verbose > 3) {
      if (this->boxes_A[u][v][w].cells[i].contacts) {
 	cout << " cell " << this->boxes_A[u][v][w].cells[i].name << " has " 
	     << this->boxes_A[u][v][w].cells[i].contacts << " contacts with ";
      }
      for (unsigned int j=0; j<this->boxes_A[u][v][w].cells[i].neighbors.size(); j++) {
	cout << this->boxes_A[u][v][w].cells[i].neighbors[j]->name << " ";
      }
      cout << endl;
    }

    // polarise cells - CICELY NEW
    //polarise(this->boxes_A[u][v][w].cells[i]);
    
    // compute cell-cell forces and new cell velocity

    // Equation (for the velocity)
    // x_dot = 1/tissue_friction * sum_{cells in contact} f_contact
    //         + 1/n_fibres_in_contact * sum_{fibres in contact} direction
    //         + polarity_velocity

    // contact forces
    hertz(this->boxes_A[u][v][w].cells[i]);
    // velocity correction (fibres and polarization)
    //if (this->fibres.size()){
      //cout << "number of fibres is: " << fibres.size() << endl;
    // cell_fibres_interaction(this->boxes_A[u][v][w].cells[i]);
    //for (int k=0; k<3; k++) {
    //	this->boxes_A[u][v][w].cells[i].polarity[k] = this->boxes_A[u][v][w].cells[i].vel[k];
    //}
    // }
    //cout << " --after fibres " << endl;
    //this->boxes_A[u][v][w].cells[0].printInfo(); 


    ///@todo this should be done somewhere else
    // check if hypoxic cells can become normoxic again
    if (this->boxes_A[u][v][w].cells[i].type==2) {
      reverse_phenotype(this->boxes_A[u][v][w].cells[i]);
    }


  }//End cells loop
}



/*****************************************************************************/
void CoupledModel::compare_elements(int u, int v, int w,
		      unsigned int cells_number,Mesh& _mesh) 
{
  
  int borderx,bordery,borderz;
  
  //Now we run over the number of cells in box
  for(unsigned int i=0;i<cells_number;i++) {
    float minimaDistanciaElem=200;
    int minimoElem=-1;
    
    for (unsigned int ii = 0; ii < this->boxes_A[u][v][w].v_triangles.size(); ii++) {
      double x,y,z,dist_tri;
      unsigned int it = this->boxes_A[u][v][w].v_triangles[ii];

      x = (_mesh.xp[_mesh.tetra[4*it]-1] +  _mesh.xp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.xp[_mesh.tetra[4*it+2]-1] + _mesh.xp[_mesh.tetra[4*it+3]-1])/4.;
      
      y = (_mesh.yp[_mesh.tetra[4*it]-1] +  _mesh.yp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.yp[_mesh.tetra[4*it+2]-1] + _mesh.yp[_mesh.tetra[4*it+3]-1])/4.;
      
      z = (_mesh.zp[_mesh.tetra[4*it]-1] +  _mesh.zp[_mesh.tetra[4*it+1]-1] 
	   +  _mesh.zp[_mesh.tetra[4*it+2]-1] + _mesh.zp[_mesh.tetra[4*it+3]-1])/4.;

      dist_tri=sqrt(pow((this->boxes_A[u][v][w].cells[i].position[0]-x),2)+
		    pow((this->boxes_A[u][v][w].cells[i].position[1]-y),2)
		    +pow((this->boxes_A[u][v][w].cells[i].position[2]-z),2));
      
      if(dist_tri<minimaDistanciaElem){
	minimaDistanciaElem=dist_tri;
	minimoElem = it;
      }
      
    }
    if (minimoElem>=0) {
      _mesh.cellsInTria[minimoElem] = _mesh.cellsInTria[minimoElem]+1;
      // ------------------------
      // fill array per type
      int cellType =  this->boxes_A[u][v][w].cells[i].type;
      if (cellType==1) {
	_mesh.cellsInTriaNorm[minimoElem] = _mesh.cellsInTriaNorm[minimoElem]+1; 
      } else if (cellType==2) {
	_mesh.cellsInTriaHypo[minimoElem] = _mesh.cellsInTriaHypo[minimoElem]+1; 
      } else {
	_mesh.cellsInTriaDead[minimoElem] = _mesh.cellsInTriaDead[minimoElem]+1; 
      }
      // ------------------------

    }
    //If there is no triangle we look in the nearby boxes
    if(minimaDistanciaElem==200){
      for(int h=-1;h<2;h++) {
	borderx=u+h;
	for(int l=-1;l<2;l++) {
	  bordery=v+l;							   
	  for(int m=-1;m<2;m++) {
	    borderz=w+m;		   		   
	    if(borderx<=(int) params.boxesx-1 && borderx>=0 && 
	       bordery<=(int) params.boxesy-1 && bordery>=0 && 
	       borderz<=(int) params.boxesz-1 && borderz>=0) //We are inside of the boxes domain
	      {//begin if  
		for (unsigned int ii = 0; ii < this->boxes_A[u+h][v+l][w+m].v_triangles.size(); ii++) {
		  double x,y,z,dist_tri;
		  unsigned int it2 = this->boxes_A[u+h][v+l][w+m].v_triangles[ii];
		  x = (_mesh.xp[_mesh.tetra[4*it2]-1] +  _mesh.xp[_mesh.tetra[4*it2+1]-1] 
		       +  _mesh.xp[_mesh.tetra[4*it2+2]-1] + _mesh.xp[_mesh.tetra[4*it2+3]-1])/4.;
		  
		  y = (_mesh.yp[_mesh.tetra[4*it2]-1] +  _mesh.yp[_mesh.tetra[4*it2+1]-1] 
		       +  _mesh.yp[_mesh.tetra[4*it2+2]-1] + _mesh.yp[_mesh.tetra[4*it2+3]-1])/4.;
		  
		  z = (_mesh.zp[_mesh.tetra[4*it2]-1] +  _mesh.zp[_mesh.tetra[4*it2+1]-1] 
		       +  _mesh.zp[_mesh.tetra[4*it2+2]-1] + _mesh.zp[_mesh.tetra[4*it2+3]-1])/4.;
		  /// @attention there was a bug here: boxes_A[u+h][v+l][w+m] instead of boxes_A[u][v][w]
		  dist_tri=sqrt(pow((this->boxes_A[u][v][w].cells[i].position[0]-x),2)+
				pow((this->boxes_A[u][v][w].cells[i].position[1]-y),2)
				+pow((this->boxes_A[u][v][w].cells[i].position[2]-z),2));
	  
		  if(dist_tri<minimaDistanciaElem){
		    minimaDistanciaElem=dist_tri;
		    minimoElem = it2;
		  }
		}
	      }
	    
	  }
	}
      }
      if (minimoElem>=0) {
	_mesh.cellsInTria[minimoElem] = _mesh.cellsInTria[minimoElem]+1; 

	// ------------------------
	// fill array per type
	int cellType =  this->boxes_A[u][v][w].cells[i].type;
	if ( cellType==1) {
	  _mesh.cellsInTriaNorm[minimoElem] = _mesh.cellsInTriaNorm[minimoElem]+1; 
	} else if (cellType==2) {
	  _mesh.cellsInTriaHypo[minimoElem] = _mesh.cellsInTriaHypo[minimoElem]+1; 
	} else {
	  _mesh.cellsInTriaDead[minimoElem] = _mesh.cellsInTriaDead[minimoElem]+1; 
	}
	// ------------------------
      }
    }
  } //for i=1..nCells 



}


/*****************************************************************************/
/// @todo write a version search_in_box(Box b, Cell c)
int CoupledModel::search_in_box(const vector<int>& box_number,
				const unsigned int cell_name) {

  int u=box_number[0];
  int v=box_number[1];
  int w=box_number[2];
  
  for(unsigned int i=0;i<this->boxes_A[u][v][w].cells.size();i++) {
    if(cell_name == this->boxes_A[u][v][w].cells[i].name)
      return i;
  }
  cout << "WARNING search_in_box(): cell " << cell_name << " not found in box "
       << u << "," << v << "," << w << endl;
  return -1;

}


/*****************************************************************************/
// MAIN LOOP
void CoupledModel::loop()
{

  /**
     @todo move this to the init function
     for this, we need to replace fe_mesh with oxy_diff.mesh everywhere
  **/
  Mesh fe_mesh;    
  if (this->params.femSolverType>0) {
    /// @todo move this part in init() function
    // ======================
    // read FE-mesh from file
    fe_mesh.read(this->params.meshdir + this->params.meshname);
    fe_mesh.info();
    // set which elements are contained in which box
    setElementsInBox(fe_mesh);
    // initialize arrays containing cell densities
    fe_mesh.initCellsArrays();
    
    // check boxes-tetra initialization
    if (this->params.verbose>2) {
      checkElementsInBoxes();
    }
    
    // ======================
    
    // set the dimension of the problem
    // according to the mesh
    if (fe_mesh.dim != (int) this->params.dimension) {
      cout << " !WARNING! input dimension = " << this->params.dimension
	   << " and mesh dimension = " << fe_mesh.dim << endl;
      cout << " I am changing the input dimension to " <<  fe_mesh.dim 
	   << " (file " << __FILE__ << ", line " 
	   << __LINE__ << ")" << endl;
    }
    
    // initialize PDE solver
    this->oxy_diff.init(this->params, fe_mesh);
    // ======================
  } else {
    // initialize empty PDE solver
    this->oxy_diff.init();
  }

  // initialize cell type counters
  this->totNorm.resize(0);
  this->totHypo.resize(0);
  this->totDead.resize(0);



  cout << " ========= " << endl;
  system(" date\n");
  cout << " starting main loop, end time: " << params.n_steps << endl;
  cout << " ========= " << endl;
  this->reloj=0; // initial time step
  

  // ====================
  // MAIN LOOP
  // ====================
  while(this->total_no_of_cells < this->max_cell && 
	this->reloj < params.n_steps) {

    //cout << "time is " << reloj << endl;

    int n_cells_old=this->total_no_of_cells;

    if (params.verbose>0) {
      if ((this->reloj%10)==0) {
	cout << " *** time step: " << this->reloj << " (of " << params.n_steps << ") " 
	     << " *** n. of cells: " << this->total_no_of_cells << endl;
      }
      // count how many cells per type (and save on file)
      if ((this->reloj%50)==0) {
	cout << " count per type " << endl;
	this->count_cells_per_type();
      }
    }
    
    // increase coupled model iteration number
    this->reloj++;
    // increase PDE iteration number
    this->oxy_diff.time = this->oxy_diff.time + 1.;

    // ********************************
    // Loop of birth
    // ********************************
    cells_counter=0;  // total number of new cells      
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++) {
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	    this->cell_birth(this->boxes_A[k][l][n].cells[i]);
	  }				
	}       
      }
    }

    // update maximum number of occupied boxes after births
    this->update_maximum();

    // update number of cels
    this->total_no_of_cells=this->total_no_of_cells+cells_counter;
    if (cells_counter) {
      if (params.verbose>0) {
	cout << " *** time step " << this->reloj << ", newborn: " << cells_counter
	     << " total number of cells: " << this->total_no_of_cells << endl;
	cout << " phenotype 1 - " << this->phenotype1_count <<  " cells "
	     << " phenotype 2 - " << this->phenotype2_count << " cells " << endl;
      }
    }

    double dist_centre_x = 0.0;
    double dist_centre_y = 0.0;
    double dist_centre_z = 0.0;
    double dist_centre = 0.0;
    double mean_dist = 0.0;
    double variance_dist = 0.0;

    for(int k=this->minx; k<=this->maxx; k++) {
	for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	      dist_centre_x = this->boxes_A[k][l][n].cells[i].position[0]-params.lattice_length_x/2.0;
	      dist_centre_y = this->boxes_A[k][l][n].cells[i].position[1]-params.lattice_length_y/2.0;
	      dist_centre_z = this->boxes_A[k][l][n].cells[i].position[2]-params.lattice_length_z/2.0;
	      
	      dist_centre = sqrt(dist_centre_x*dist_centre_x+dist_centre_y*dist_centre_y+dist_centre_z*dist_centre_z);

	      mean_dist += dist_centre;
	    }}}}
    mean_dist = mean_dist/this->total_no_of_cells;

     for(int k=this->minx; k<=this->maxx; k++) {
	for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	      dist_centre_x = this->boxes_A[k][l][n].cells[i].position[0]-params.lattice_length_x/2.0;
	      dist_centre_y = this->boxes_A[k][l][n].cells[i].position[1]-params.lattice_length_y/2.0;
	      dist_centre_z = this->boxes_A[k][l][n].cells[i].position[2]-params.lattice_length_z/2.0;
	      
	      dist_centre = sqrt(dist_centre_x*dist_centre_x+dist_centre_y*dist_centre_y+dist_centre_z*dist_centre_z);

	      variance_dist += (dist_centre-mean_dist)*(dist_centre-mean_dist);
	    }}}}
     variance_dist = variance_dist/this->total_no_of_cells;


     if (this->params.verbose>1){
    // =================================================================
    // write mean position (relative to centre of domain) to data file
    // =================================================================
    ofstream meandist;
    string meandist_list =  this->params.casedirectory + this->params.casename + "_mean_distance.txt";
    meandist.open(meandist_list.c_str(),ios::app);
    meandist << mean_dist << " " ;
    meandist.close();
    // ********************************

    // ============================================
    // write variance of above to data file
    // ============================================
    ofstream vardist;
    string vardist_list =  this->params.casedirectory + this->params.casename + "_variance_distance.txt";
    vardist.open(vardist_list.c_str(),ios::app);
    vardist << variance_dist << " " ;
    vardist.close();
    // ********************************
  }
    
    // ============================================
    // write intermediary cell numbers to data file
    // ============================================
    ofstream pheno1File;
    string full_list =  this->params.casedirectory + this->params.casename + "_full.txt";
    pheno1File.open(full_list.c_str(),ios::app);
    pheno1File << this->phenotype1_count << " " ;
    pheno1File.close();
    // if (params.ic_phenotype.size()==2){
    //   ofstream pheno2File;
    //   pheno2File.open("phenotype2_cellnumbers.txt",ios::app);
    //   pheno2File << this->phenotype2_count << " " ;
    //   pheno2File.close();
    // }
    // ********************************

    // ********************************
    // Loop of movement
    // ********************************
    // count cells for each box (before movement)
    if (this->params.verbose>0) {
      this->count_cells_per_box();
    }
    
    //clock_t t = clock();
    /*
    //clock_t t = clock();
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++) {
	for(int n=this->minz; n<=this->maxz; n++) {
	  // find fibres in contact with cells in the box (k,l,n)
	  this->compute_cell_fibres_contact(k,l,n);
	}
      }
    }
    //t = clock() - t;
    //cout << " total Wall Clock time:" << ((float)t)/CLOCKS_PER_SEC << endl;
    */
    
    // Compute contacts and  forces
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++) {
	for(int n=this->minz; n<=this->maxz; n++) {

	  this->compute_cell_fibres_contact(k,l,n);
	  this->compute_cell_vessels_contact(k,l,n);
	  
	  // the number of cells in the box must be stored since it
	  // might change during the contact force computation ?
	  initial_cells_in_box=this->boxes_A[k][l][n].cells.size();
	  this->contact_forces(k,l,n,initial_cells_in_box);

	}
      }
    }

    //t = clock() - t;
    //cout << " total Wall Clock time (contacts):" << ((float)t)/CLOCKS_PER_SEC << endl;
    
    // move cells
    // set to zero the counter of cells in each tria
    fe_mesh.initCellsArrays();
    unsigned int m_count=0;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++) {
	for(int n=this->minz; n<=this->maxz; n++) {

	  if(this->params.femSolverType>0) {
	    // -------------------------------
	    // fill the array: fe_mesh.cellsInTria
	    this->compare_elements(k,l,n,this->boxes_A[k][l][n].cells.size(),fe_mesh);
	  }
	  // -------------------------------
	  // move cells
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    //cout << " move cell " << boxes_A[k][l][n].cells[i].name << endl;
	    m_count++;
	    this->movement(this->boxes_A[k][l][n].cells[i],k,l,n,i);
	  }           
	}
      }
    } 
   

    // // ================================
    // // COMMUNICATION: write cell output
    // // ================================    
    // // write cell file
    // ofstream outCellFile;
    // outCellFile.open(this->params.fileCells2FEM.c_str(),ios::out);
    // if (!outCellFile) {
    //   cerr << " *** ERROR, file " << __FILE__ << ", line " << __LINE__ 
    // 	   << " *** could not open file " << this->params.fileCells2FEM << endl;
    //   exit(1);
    // }

    // // header: length, iteration (pde), iteration (total)
    // outCellFile << "# first line: n_cells pde_iteration cells_iteration x_size y_size z_size phenotype_count_1 phenotype_count_2" << endl;
    // outCellFile << "# then: cell database" << endl;
    // outCellFile << this->total_no_of_cells << " "
    // 		<< this->oxy_diff.iter << " "
    // 		<< this->reloj << " "
    // 		<< this->maxx-this->minx << " " << this->maxy-this->miny << " " << this->maxz-this->minz << " "
    // 		<< this->phenotype1_count << " " << this->phenotype2_count << endl;
    // // body: x,y,z,type,r,phenotype,adhesion coeff,id,energy
    // for(int k=this->minx; k<=this->maxx; k++) {
    //   for(int l=this->miny; l<=this->maxy; l++){
    // 	for(int n=this->minz; n<=this->maxz; n++) {
    // 	  for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
    // 	    outCellFile << this->boxes_A[k][l][n].cells[i].position[0] << " "
    // 			<< this->boxes_A[k][l][n].cells[i].position[1] << " "
    // 			<< this->boxes_A[k][l][n].cells[i].position[2] << " "
    // 			<< this->boxes_A[k][l][n].cells[i].type << " "
    // 			<< this->boxes_A[k][l][n].cells[i].radius << " "
    // 			<< this->boxes_A[k][l][n].cells[i].phenotype << " "
    // 			<< this->boxes_A[k][l][n].cells[i].adhesion << " "
    // 			<< this->boxes_A[k][l][n].cells[i].name << " "
    // 			<< this->boxes_A[k][l][n].cells[i].energy << endl;
    // 	  }
    // 	}
    //   }
    // }
    // outCellFile.close();

    
    if(this->params.femSolverType==0) {

      // write output list of cells
      std::stringstream outputFileName;
      std::string s0 = this->params.outputDirectory + this->params.testcase + "_cells.";
      outputFileName << s0  << reloj << ".vtk";
      this->writeVtk(outputFileName.str());      

      if (params.writeVtkFiles) {
	s0 = this->params.outputDirectory + this->params.testcase + "_fibres.";
	std::stringstream fibres_outputFileName;
	fibres_outputFileName << s0  << reloj << ".vtk";
	this->writeFibresVtk(fibres_outputFileName.str());

	s0 = this->params.outputDirectory + this->params.testcase + "_vessels.";
	std::stringstream vessels_outputFileName;
	vessels_outputFileName << s0  << reloj << ".vtk";
	this->writeVesselsVtk(vessels_outputFileName.str());

	
      }
    } else {

      // =================================
      // write piecewise constant cell density (-> FreeFem)
      ofstream o_file;
      o_file.open(this->params.fileCellsDensity2FEM.c_str(),ios::out);
      for (unsigned int i=0; i<fe_mesh.nElem; i++) {
	o_file << fe_mesh.cellsInTria[i] << " "
	       << fe_mesh.cellsInTriaNorm[i] << " " // normoxic
	       << fe_mesh.cellsInTriaHypo[i] << " " // hypoxic
	       << fe_mesh.cellsInTriaDead[i] << " " // dead
	       << endl;
      }
      o_file.close();
      
      // =================================
      // PDE solver
      // =================================
      
      /// @todo generalize the flag for solving PDE
      this->oxy_diff.launch = ((this->reloj%100)==0) || (this->total_no_of_cells-n_cells_old>30);
      
      if (this->oxy_diff.launch) {
	cout << " solve PDE " << endl;
	// solve PDE (and write new concentration)
	this->oxy_diff.solve();
	
	// ------------------------------
	// read the new concentration file
	ifstream o2_conc_file;
	string ifname = this->params.fileFEM2Cells;
	o2_conc_file.open(ifname.c_str(),ios::in);
	
	unsigned int ic = 0;
	for(int k=this->minx; k<=this->maxx; k++) {
	  for(int l=this->miny; l<=this->maxy; l++){
	    for(int n=this->minz; n<=this->maxz; n++) {
	      for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
		o2_conc_file >> this->boxes_A[k][l][n].cells[i].O2;
		o2_conc_file >> this->boxes_A[k][l][n].cells[i].dxO2;
		o2_conc_file >> this->boxes_A[k][l][n].cells[i].dyO2;
		o2_conc_file >> this->boxes_A[k][l][n].cells[i].dzO2;
		ic++;
	      }
	    }
	  }
	}
	o2_conc_file.close();
	// ------------------------------
      } // if (oxy_diff.launch)   
      
      
    }
    
    // print all cells infos
    if (this->params.verbose>4) {
      for(int k=this->minx; k<=this->maxx; k++) {
	for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	      this->boxes_A[k][l][n].cells[i].printInfo();
	    }
	  }
	}
      }
    }

    this->update_maximum();
    this->update_box();
    // ********************************
   
    
  }   //main loop while

  // ============================================
    // parameter list
    // ============================================
    ofstream parameters;
    string paramlist = this->params.casedirectory + this->params.casename + "_parameters.txt";
    parameters.open(paramlist.c_str(),ios::out);
    parameters << "The parameters for the case " << this->params.casename.c_str() << " are:" << endl;
    parameters << "growth rate = " << params.growth_rate[0] << endl <<
      "birth rate = " << params.alpha_birthrate[0] << endl <<
      "Youngs Modulus = " << params.alpha_YoungM[0] << endl <<
      "Poisson number = " << params.alpha_PoissonNo[0] << endl <<
      "GCM = " << params.alpha_gcm[0]*params.Gcm << endl <<
      "adhesion value = " << params.adhesion_value[0] << endl <<
      "random motion = " << params.variance_motion << endl <<
      "birth energy = " << this->birth_energy << endl;
    parameters.close();
    // ********************************

if (this->params.verbose>0){
    // ============================================
    // write mean position to data file
    // ============================================
    ofstream meandist;
    string meandist_list =  this->params.casedirectory + this->params.casename + "_mean_distance.txt";
    meandist.open(meandist_list.c_str(),ios::app);
    meandist << endl;
    meandist.close();
    // ********************************

    // ============================================
    // write variance of above to data file
    // ============================================
    ofstream vardist;
    string vardist_list =  this->params.casedirectory + this->params.casename + "_variance_distance.txt";
    vardist.open(vardist_list.c_str(),ios::app);
    vardist << endl;
    vardist.close();
    // ********************************
 }


  // ============================================
    // intermediary cell numbers end line
    // ============================================
    ofstream pheno1File;
    string full_list = this->params.casedirectory + this->params.casename + "_full.txt";
    pheno1File.open(full_list.c_str(),ios::app);
    pheno1File << endl;
    pheno1File.close();
    // if (params.ic_phenotype.size()==2){
    //   ofstream pheno2File;
    //   pheno2File.open("phenotype2_cellnumbers.txt",ios::app);
    //   pheno2File << endl;
    //   pheno2File.close();
    // }
    // ********************************

  // ======================================
    // write final cell numbers to data file
    // ====================================== 
    ofstream cellnumbersFile;
    string final_list = this->params.casedirectory + this->params.casename + "_final_timestep.txt";
    cellnumbersFile.open(final_list.c_str(),ios::app);
    cellnumbersFile << this->total_no_of_cells << " "
		<< this->reloj << " "
		<< this->maxx-this->minx << " " << this->maxy-this->miny << " " << this->maxz-this->minz << " "
		    << this->phenotype1_count << " " << this->phenotype2_count << " "
		    << (this->maxx-this->minx)*(this->maxy - this->miny)*(this->maxz - this->minz) << endl;
    cellnumbersFile.close();
    
    if (params.alpha_birthrate[0]==0){
      ofstream cellpositionFile;
      cellpositionFile.open("cell_positions.txt",ios::app);
      for(int k=this->minx; k<=this->maxx; k++) {
	for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	    for(unsigned int i=0; i<this->boxes_A[k][l][n].cells.size(); i++) {
	      cellpositionFile << this->boxes_A[k][l][n].cells[i].position[0]-params.ic_cell_x[i] << " "
			       << this->boxes_A[k][l][n].cells[i].position[1]-params.ic_cell_y[i] << " "
			       << this->boxes_A[k][l][n].cells[i].position[2]-params.ic_cell_z[i] << " " << endl;
	    }}}}
      cellpositionFile.close();
      }
}


void CoupledModel::count_cells_per_box()
{
      for(int k=this->minx; k<=this->maxx; k++) {
	for(int l=this->miny; l<=this->maxy; l++){
	  for(int n=this->minz; n<=this->maxz; n++) {
	if (this->boxes_A[k][l][n].cells.size()) {
	  cout << " box [" << k << "," << l << "," << n << "] has " 
	       << this->boxes_A[k][l][n].cells.size() << " cells " << endl;
	}
      }
    }
  }
  
}

void CoupledModel::count_cells_per_type()
{
  unsigned int countNorm = 0;
  unsigned int countHypo = 0;
  unsigned int countDead = 0;
  unsigned int countPhenotypeOxygen_0_10 = 0;
  unsigned int countPhenotypeOxygen_0_12 = 0;
  unsigned int countPhenotypeOxygen_14_24 = 0;
  unsigned int countPhenotypeOxygen_12_24 = 0;
  
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	  
	  unsigned int cellType = this->boxes_A[k][l][n].cells[i].type;
	  countNorm = countNorm + (cellType==1);
	  countHypo = countHypo + (cellType==2);
	  countDead = countDead + (cellType==3);
	  if (this->boxes_A[k][l][n].cells[i].phenotype < 10) {
	    countPhenotypeOxygen_0_10 ++;
	  } else if ( (this->boxes_A[k][l][n].cells[i].phenotype > 14) &&
		      (this->boxes_A[k][l][n].cells[i].phenotype < 24) ){
	    countPhenotypeOxygen_14_24 ++;
	  }
	  
	  if (this->boxes_A[k][l][n].cells[i].phenotype < 12) {
	    countPhenotypeOxygen_0_12 ++;
	  } else if ( (this->boxes_A[k][l][n].cells[i].phenotype > 12) &&
		      (this->boxes_A[k][l][n].cells[i].phenotype < 24) ){
	    countPhenotypeOxygen_12_24 ++;}
	  
	}
	
      }
    }
  }
  this->totNorm.push_back(countNorm);
  this->totHypo.push_back(countHypo);
  this->totDead.push_back(countDead);
  
  ofstream nCellFile;
  nCellFile.open("cell_counter.txt",ios::app);
  if (!nCellFile) {
    cerr << " *** ERROR *** could not open file *** " << endl;
    exit(1);
  }
  
  unsigned int nt = this->totNorm.size()-1;
  nCellFile << this->totNorm[nt] << " " 
	    << this->totHypo[nt] << " " << this->totDead[nt] <<  endl;
  cout << " cells (norm,hypo,dead): " 
       << this->totNorm[nt] << " " << this->totHypo[nt] << " " << this->totDead[nt] <<  endl;
  if ( (this->totNorm[nt] + this->totHypo[nt] + this->totDead[nt]) < this->total_no_of_cells ) {
    cout << " ERROR in CoupledModel::count_cells_per_type(): not all cells found ! " << endl;
    this->end();
    exit(1);
  }

}
  

void CoupledModel::end() {

  std::string s0;
  if (params.writeVtkFiles == 0) {
    s0 = this->params.outputDirectory + this->params.testcase + "_fibres.";
    std::stringstream fibres_outputFileName;
    fibres_outputFileName << s0  << reloj << ".vtk";
    this->writeFibresVtk(fibres_outputFileName.str());
  }
  s0 = this->params.outputDirectory + this->params.testcase + "_boxes.vtk";
  this->writeBoxesVtk(s0);
}

void CoupledModel::writeVtk(string filename,unsigned int onlyCoord)
{

  ofstream outfile(filename.c_str());

  // rename (locally) the total number of cells
  unsigned int nCells = this->total_no_of_cells;
  
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Unstructured grid legacy vtk file with point scalar data" << endl;
  outfile << "ASCII\n\n";
  
  outfile << "DATASET UNSTRUCTURED_GRID\n";
  outfile << "POINTS " << nCells << " double\n";
  
  // write positions: we loop on all boxes
  //@warning the 2-dimensional output is not supported
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	  outfile << this->boxes_A[k][l][n].cells[i].position[0] << " "
		  << this->boxes_A[k][l][n].cells[i].position[1] << " "
		  << this->boxes_A[k][l][n].cells[i].position[2] << endl;
	}
      }
    }
  }
  
  outfile << endl;

  outfile << "POINT_DATA " << nCells << endl;
  // write radii
  outfile << "SCALARS radius double" << endl;
  outfile << "LOOKUP_TABLE default" << endl;
  for(int k=this->minx; k<=this->maxx; k++) {
    for(int l=this->miny; l<=this->maxy; l++){
      for(int n=this->minz; n<=this->maxz; n++) {
	for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	  outfile << this->boxes_A[k][l][n].cells[i].radius << endl;
	}
      }
    }
  }
  outfile << endl;

  if (onlyCoord==0) {
    // write cell type
    outfile << "SCALARS status double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].type << endl;
	  }
	}
      }
    }
    outfile << endl;
    
    // write oxygen concentration
    outfile << "SCALARS concentration double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].O2 << endl;
	  }
	}
      }
    }
    outfile << endl;

    //write fibre contacts
    outfile << "SCALARS fibres double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].contact_fibres.size() << endl;
	  }
	}
      }
    }
    outfile << endl;
    
    //write phenoype
    outfile << "SCALARS phenot double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].phenotype << endl;
	  }
	}
      }
    }
    outfile << endl;
    
    
    
    //write adhesion
    outfile << "SCALARS name double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].name << endl;
	  }
	}
      }
    }
    outfile << endl;

     //write interaction phenotype
    outfile << "SCALARS interaction_phenotype double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].interaction_phenotype << endl;
	  }
	}
      }
    }
    outfile << endl;

    //write interaction phenotype
    outfile << "SCALARS follow_lead double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].follow_lead << endl;
	  }
	}
      }
    }
    outfile << endl;


     //write polarised - CICELY NEW
    outfile << "SCALARS polarised double" << endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for(int k=this->minx; k<=this->maxx; k++) {
      for(int l=this->miny; l<=this->maxy; l++){
	for(int n=this->minz; n<=this->maxz; n++) {
	  for(unsigned int i=0;i<this->boxes_A[k][l][n].cells.size();i++) {
	    outfile << this->boxes_A[k][l][n].cells[i].polarised << endl;
	  }
	}
      }
    }
    outfile << endl;
    
  }


}

void CoupledModel::writeBoxesVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  // header
  ofile << "# vtk DataFile Version 3.0" << endl;
  ofile << "Boxes output - vtk " << endl;
  ofile << "ASCII" << endl;
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  // nodes
  int nNodesX = params.boxesx+1;
  int nNodesY = params.boxesy+1;
  int nNodesZ = params.boxesz+1;
  int nNodes = nNodesX*nNodesY*nNodesZ;
  ofile << "POINTS " << nNodes << " float" << endl;
  for (int k=0; k<nNodesZ; k++){
    for (int j=0; j<nNodesY; j++){
      for (int i=0; i<nNodesX; i++){
	ofile << i*this->box_sizex << " " << j*this->box_sizey << " "
	      <<  k*this->box_sizez << endl;
      }
    }  
  }
  // cells
  int nCells = (params.boxesx)*(params.boxesy)*(params.boxesz);
  ofile << "CELLS " << nCells << " " << 9*nCells << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	// coordinates of box_(i,j,k)
	// the points must be ordered according to vtk data structures
	// see e.g. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
	ofile << 8 << " " <<  i + j*nNodesX + k*nNodesX*nNodesY << " "
	      << (i+1) + j*nNodesX + k*nNodesX*nNodesY << " "
	      << i + (j+1)*nNodesX + k*nNodesX*nNodesY << " "
	      << (i+1) + (j+1)*nNodesX + k*nNodesX*nNodesY << " "
	      <<  i + j*nNodesX + (k+1)*nNodesX*nNodesY << " "
	      << (i+1) + j*nNodesX + (k+1)*nNodesX*nNodesY << " "
	      << i + (j+1)*nNodesX + (k+1)*nNodesX*nNodesY << " "
	      << (i+1) + (j+1)*nNodesX + (k+1)*nNodesX*nNodesY
	      << endl;
      }
    }  
  }
  // nodes
  ofile << "CELL_TYPES " << nCells  << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	ofile << 11 << endl;
      }
    }  
  }
  ofile << "CELL_DATA  " << nCells << endl;
  ofile << "SCALARS nCells float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int k=0; k<params.boxesz; k++){
    for (unsigned int j=0; j<params.boxesy; j++){
      for (unsigned int i=0; i<params.boxesx; i++){
	ofile << this->boxes_A[i][j][k].cells.size()<< endl;
      }
    }
  }
  ofile.close();
}

// CICELY NEW
void CoupledModel::writeFibresVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  // file version and identifier in the form "# vtk DataFile Version x.x"
  ofile << "# vtk DataFile Version 3.0" << endl;
  // header used to describe the data
  ofile << "Fibres output - vtk " << endl;
  // file format either "ASCII" or "BINARY"
  ofile << "ASCII" << endl;
  // dataset structure in the form "DATASET_type"
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  // dataset attributes using "POINTS " or "CELLS " followed by an integer number specifying the number of points or cells.
  ofile << "POINTS " << 2*fibres.size() << " double" << endl; // Each fibre has two points (start and end)
  for(unsigned int ff=0; ff < fibres.size(); ff++)
  {
    ofile << this->fibres[ff].start[0] << " " << this->fibres[ff].start[1] << " " << this->fibres[ff].start[2] << endl;
    ofile << this->fibres[ff].start[0]+fibres[ff].length*fibres[ff].direction[0] << " " << this->fibres[ff].start[1]+fibres[ff].length*fibres[ff].direction[1] << " " << this->fibres[ff].start[2]+fibres[ff].length*fibres[ff].direction[2] << endl;
  }
  ofile << "CELLS " << fibres.size() << " " << 3*fibres.size() << endl; // Each Fibre element has to detail its type and its points
  // 	// coordinates of each fibre
  // 	// the points must be ordered according to vtk data structures
  // 	// see e.g. http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
  for (unsigned int ff=0; ff<fibres.size(); ff++)
  {
    ofile << 2 << " " << 2*ff << " " << 2*ff + 1 << endl;
  }
  ofile << "CELL_TYPES " << fibres.size() << endl;
  for(unsigned int ff=0; ff < fibres.size(); ff++)
  {
    ofile << 3 << endl;
  }
  ofile << "CELL_DATA  " << fibres.size() << endl;
  ofile << "SCALARS length float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++)
  {
    ofile << fibres[ff].length << endl;
  }
  ofile << "SCALARS fibre_name double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++)
  {
   ofile << fibres[ff].fibre_name << endl;
  }
  ofile << "SCALARS fibre_exists double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++)
  {
   ofile << fibres[ff].fibre_exists << endl;
  }
  ofile << "VECTORS orientation double" << endl;
  for (unsigned int ff=0; ff < fibres.size(); ff++)
  {
   ofile << fibres[ff].direction[0] << " " << fibres[ff].direction[1] << " " << fibres[ff].direction[2] << endl;
  }
  ofile.close();
}

// *********** Write Vessels Vtk ********************
void CoupledModel::writeVesselsVtk(string filename)
{
  ofstream ofile;
  ofile.open(filename.c_str());
  ofile << "# vtk DataFile Version 3.0" << endl;
  ofile << "Vessels output - vtk " << endl;
  ofile << "ASCII" << endl;
  ofile << "DATASET UNSTRUCTURED_GRID" << endl;
  ofile << "POINTS " << 2*vessels.size() << " double" << endl; // Each vessel has two points (start and end)
  for(unsigned int ff=0; ff < vessels.size(); ff++)
  {
    ofile << this->vessels[ff].ves_start[0] << " " << this->vessels[ff].ves_start[1] << " " << this->vessels[ff].ves_start[2] << endl;
    ofile << this->vessels[ff].ves_start[0]+vessels[ff].ves_length*vessels[ff].ves_direction[0] << " " << this->vessels[ff].ves_start[1]+vessels[ff].ves_length*vessels[ff].ves_direction[1] << " " << this->vessels[ff].ves_start[2]+vessels[ff].ves_length*vessels[ff].ves_direction[2] << endl;
  }
  ofile << "CELLS " << vessels.size() << " " << 3*vessels.size() << endl; // Each Vessel element has to detail its type and its points
  for (unsigned int ff=0; ff<vessels.size(); ff++)
  {
    ofile << 2 << " " << 2*ff << " " << 2*ff + 1 << endl;
  }
  ofile << "CELL_TYPES " << vessels.size() << endl;
  for(unsigned int ff=0; ff < vessels.size(); ff++)
  {
    ofile << 3 << endl;
  }
  ofile << "CELL_DATA  " << vessels.size() << endl;
  ofile << "SCALARS length float" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++)
  {
    ofile << vessels[ff].ves_length << endl;
  }
  ofile << "SCALARS vessel_name double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
  for (unsigned int ff=0; ff < vessels.size(); ff++)
  {
   ofile << vessels[ff].vessel_name << endl;
  }
  ofile << "SCALARS vessel_radius double" << endl;
  ofile << "LOOKUP_TABLE default" << endl;
 for (unsigned int ff=0; ff < vessels.size(); ff++)
  {
    ofile << vessels[ff].ves_radius << endl;
  }
  ofile.close();
}

void CoupledModel::fibre_degradation(Cell& cell, Fibre& fibre, double cell_fibre_min_dist) {
  
  double c_to_f_start = 0.0;
  double c_to_f_end = 0.0;
  for (unsigned int l=0; l<3; l++) {
    c_to_f_start += (cell.position[l]-fibre.start[l])*cell.polarity[l];
    c_to_f_end += (cell.position[l]-(fibre.start[l]+fibre.length*fibre.direction[l]))*cell.polarity[l];
  }
  
  if (c_to_f_start>0.0 || c_to_f_end>0.0) {
    
    if(cell_fibre_min_dist < cell.radius)
      {

	double velocity_dot_direction = 0.0;
	double cell_velocity = 0.0;
	for (unsigned int k=0; k<3; k++) {
	  velocity_dot_direction += fibre.direction[k]*cell.polarity[k];
	  cell_velocity += cell.polarity[k]*cell.polarity[k];
	}
	cell_velocity = sqrt(cell_velocity);
	velocity_dot_direction = velocity_dot_direction/cell_velocity;

	//if (fabs(velocity_dot_direction)<sqrt(2.0)/2.0){
	  double rand_contact_degrade = aleatorio();
	  if (rand_contact_degrade<0.5){
	    fibre.fibre_exists=0.0;
	    cout << "fibre " << fibre.fibre_name << " eliminated by contact" << endl;
	  }
	  //}
	
      }
  }
  
  double rand_diffuse_degrade = aleatorio();
  if (rand_diffuse_degrade<0.5){
    fibre.fibre_exists=0.5;
    cout << "fibre " << fibre.fibre_name << " eliminated by diffusion" << endl;
  }
  
}		    


int CoupledModel::get_sub_domain_id(double x,double y,double z)
{
  int sub_domain_model_type = 0;
  
  switch (sub_domain_model_type) {
  case 0:
    {  
      if (x <= params.lattice_length_x/2.0) 
	return 0;
      else
	return 1;
      
      break;
    }
  default:
    {
      cout << " ** ERROR in CoupledModel::get_sub_domain_id: sub_domain_model_type = "
	   << sub_domain_model_type << " not implemented yet" << endl;
      exit(1);
      return -1;
      break;
    }
    return -1;
  }
 
}
