#include "Cell.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

Cell::Cell()
{
  // allocate double vectors when creating new cell
  box.resize(3);
  new_box.resize(3);
  position.resize(3);
  position_old.resize(3);
  vel.resize(3);
  polarity.resize(3);

  contacts = 0;
  neighbors.clear();
  contact_fibres.clear();
  contact_vessels.clear();
}

void Cell::printInfo() const 
{
  
  cout << " *** print cell " << name << " info *** " << endl;
  cout << " - position: " << position[0] << "," << position[1] << "," << position[2] << endl;
  cout << " - velocity: " << vel[0] << "," << vel[1] << "," << vel[2] << endl;
  cout << " - contacts: " << contacts;
  if (contacts) {
    cout << " *** with cell: ";
    for (unsigned int i = 0; i < contacts ; i++) {
      cout <<  neighbors[i]->name << ", ";
    }
    // cout << " *** with fibre: ";
    // for (unsigned int i = 0; i < contacts ; i++) {
    //   cout <<  contact_fibres[i]->name << ", ";
    // }
  }
  cout << endl;
  cout << " - box: " << box[0] << " " << box[1] << " " << box[2] << endl;
  cout << " - type: " << type << endl;
  cout << " - radius: " << radius << endl;
  //cout << " - polarised: " << polarised << endl;
  cout << " - energy: " << energy << endl;
  cout << " - phenotype: " << phenotype << endl;
  cout << " - adhesion: " << adhesion << endl;
  //cout << " O2: " << O2 << ", grad = " << dxO2 <<" " << dyO2 < " " << dzO2 << endl;
  cout << " ***********************\n " << endl;

}

void Cell::clear_contacts(){
    contacts = 0;  
    // note: swap frees the memory and removes the elements
    vector<Cell*>().swap(neighbors);
    vector<Fibre*>().swap(contact_fibres);
    vector<int>().swap(contact_vessels);
}
