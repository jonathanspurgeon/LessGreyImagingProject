#include "Fibre.h"
#include <stdlib.h>
#include <iostream>

using namespace std;

Fibre::Fibre()
{
  //fibre_clear_contacts();
}

void Fibre::fibre_printInfo() const 
{
  
  cout << " *** print fibre " << fibre_name << " info *** " << endl;
  cout << " - start position: " << start[0] << "," 
        << start[1] << "," << start[2] << endl;
  cout << " - length: " << length << endl;
  cout << " - direction: " << direction[0] << "," 
        << direction[1] << "," << direction[2] << endl;
  // cout << " - contacts: " << contacts;
  // if (contacts) {
  //   cout << " *** with cell: ";
  //   for (int i = 0; i < contacts ; i++) {
  //     cout <<  neighbors[i]->name << ", ";
  //   }
  // }
  // cout << endl;
  // cout << " - box: " << box[0] << " " << box[1] << " " 
  //      << box[2] << endl;
  // cout << " - type: " << type << endl;
  // cout << " - radius: " << radius << endl;
  // cout << " - energy: " << energy << endl;
  // cout << " - phenotype: " << phenotype << endl;
  // cout << " - adhesion: " << adhesion << endl;
  // cout << " O2: " << O2 << ", grad = " << dxO2 <<" " << dyO2 
  //      << " " << dzO2 << endl;
  cout << " ***********************\n " << endl;
  
}


void Fibre::fibre_clear_contacts()
{
  // fibre_contacts = 0;
  // fibre_neighbors.resize(0);
}
