#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

#define BLOCKDIMENTION 4

/**
 * Number of cell we have per axis
 */
int numberOfCellsPerAxisX;
int numberOfCellsPerAxisY;
int numberOfCellsPerAxisZ;

/**
 * Velocity of fluid along x,y and z direction
 */
double* ux;
double* uy;
double* uz;

/**
 * Helper along x,y and z direction
 */
double* Fx;
double* Fy;
double* Fz;

/*
 * Array to store if the block contains a cell that is in the object.
 */
bool* objectInBlock;

/**
 * Pressure in the cell
 */
double* p;
double* rhs;

/**
 * A marker that is required for the Scientific Computing module.
 */
double* ink;


/**
 * Is cell inside domain
 */
bool* cellIsInside;

double timeStepSize;

double ReynoldsNumber = 0;

const int MaxComputePIterations                  = 20000;
const int MinAverageComputePIterations           = 200;
const int IterationsBeforeTimeStepSizeIsAltered  = 64;
const double ChangeOfTimeStepSize                = 0.1;
const double PPESolverThreshold                  = 1e-6;


/**
 * Switch on to have a couple of security checks
 */
#define CheckVariableValues


/**
 * This is a small macro that I use to ensure that something is valid. There's
 * a more sophisticated version of this now in the C++ standard (and there are
 * many libs that have way more mature assertion functions), but it does the
 * job. Please comment it out if you do runtime measurements - you can do so by
 * commenting out the define statement above.
 */
void assertion( bool condition, int line ) {
  #ifdef CheckVariableValues
  if (!condition) {
    std::cerr << "assertion failed in line " << line << std::endl;
    exit(-1);
  }
  #endif
}

/**
 * Maps the three coordinates onto one cell index.
 * This is converting from the raw version to normal one
 */
int getCellIndex(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<numberOfCellsPerAxisX+2,__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<numberOfCellsPerAxisY+2,__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<numberOfCellsPerAxisZ+2,__LINE__);
  return ix+iy*(numberOfCellsPerAxisX+2)+iz*(numberOfCellsPerAxisX+2)*(numberOfCellsPerAxisY+2);
}

/**
 * Maps the three coordinates onto one vertex index.
 * Please note that we hold only the inner and boundary vertices.
 */
int getVertexIndex(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<numberOfCellsPerAxisX+1,__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<numberOfCellsPerAxisY+1,__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<numberOfCellsPerAxisZ+1,__LINE__);
  return ix+iy*(numberOfCellsPerAxisX+1)+iz*(numberOfCellsPerAxisX+1)*(numberOfCellsPerAxisY+1);
}

/**
 * Maps the three coordinates onto the corresponding block.
 * ADDED BY ME
 * This will take in a raw number that will need to be converted to the corresponding block
 */
int getBlockIndexFromCellCoordinates(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<(numberOfCellsPerAxisX+2),__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<(numberOfCellsPerAxisY+2),__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<(numberOfCellsPerAxisZ+2),__LINE__);
  return (ix-1)/BLOCKDIMENTION+(iy-1)/BLOCKDIMENTION*(numberOfCellsPerAxisX+2)+(iz-1)/BLOCKDIMENTION*(numberOfCellsPerAxisX+2)*(numberOfCellsPerAxisY+2);
}

/**
 * Maps the three coordinates onto the corresponding block.
 * ADDED BY ME
 * This will take in a raw number that will need to be converted to the corresponding block
 */
int getBlockIndexFromBlockCoordinates(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<(numberOfCellsPerAxisX+2)/BLOCKDIMENTION,__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<(numberOfCellsPerAxisY+2)/BLOCKDIMENTION,__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<(numberOfCellsPerAxisZ+2)/BLOCKDIMENTION,__LINE__);
  return (ix/BLOCKDIMENTION)+(iy/BLOCKDIMENTION)*(numberOfCellsPerAxisX+2)+(iz/BLOCKDIMENTION)*(numberOfCellsPerAxisX+2)*(numberOfCellsPerAxisY+2);
}

/**
 * Gives you the face with the number ix,iy,iz.
 *
 * Takes into account that there's one more face in X direction than numberOfCellsPerAxisX.
 */
int getFacex_cell_WithHalo(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<numberOfCellsPerAxisX+3,__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<numberOfCellsPerAxisY+2,__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<numberOfCellsPerAxisZ+2,__LINE__);
  return ix+iy*(numberOfCellsPerAxisX+3)+iz*(numberOfCellsPerAxisX+3)*(numberOfCellsPerAxisY+2);
}


int getFacey_cell_WithHalo(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<numberOfCellsPerAxisX+2,__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<numberOfCellsPerAxisY+3,__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<numberOfCellsPerAxisZ+2,__LINE__);
  return ix+iy*(numberOfCellsPerAxisX+2)+iz*(numberOfCellsPerAxisX+2)*(numberOfCellsPerAxisY+2);
}


int getFacez_cell_WithHalo(int ix, int iy, int iz) {
  assertion(ix>=0,__LINE__);
  assertion(ix<numberOfCellsPerAxisX+2,__LINE__);
  assertion(iy>=0,__LINE__);
  assertion(iy<numberOfCellsPerAxisY+2,__LINE__);
  assertion(iz>=0,__LINE__);
  assertion(iz<numberOfCellsPerAxisZ+3,__LINE__);
  return ix+iy*(numberOfCellsPerAxisX+2)+iz*(numberOfCellsPerAxisX+2)*(numberOfCellsPerAxisY+2);
}


/**
 * We use always numberOfCellsPerAxisX=numberOfCellsPerAxisY and we make Z 5
 * times bigger, i.e. we always work with cubes. We also assume that the whole
 * setup always has exactly the height 1.
 */
double getH() {
  return 1.0/static_cast<double>(numberOfCellsPerAxisY);
}


/**
 * There are two types of errors/bugs that really hunt us in these codes: We
 * either have programmed something wrong (use wrong indices) or we make a
 * tiny little error in one of the computations. The first type of errors is
 * covered by assertions. The latter type is realised in this routine, where
 * we do some consistency checks.
 */
void validateThatEntriesAreBounded(const std::string&  callingRoutine) {
  #ifdef CheckVariableValues
  for (int ix=0; ix<numberOfCellsPerAxisX+2; ix++)
  for (int iy=0; iy<numberOfCellsPerAxisY+2; iy++)
  for (int iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
    if ( std::abs(p[ getCellIndex(ix,iy,iz)])>1e10 ) {
      std::cerr << "error in routine " + callingRoutine + " in p[" << ix << "," << iy << "," << iz << "]" << std::endl;
      exit(-1);
    }
  }

  for (int ix=0; ix<numberOfCellsPerAxisX+3; ix++)
  for (int iy=0; iy<numberOfCellsPerAxisY+2; iy++)
  for (int iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
    if ( std::abs(ux[ getFacex_cell_WithHalo(ix,iy,iz)])>1e10 ) {
      std::cerr << "error in routine " + callingRoutine + " in ux[" << ix << "," << iy << "," << iz << "]" << std::endl;
      exit(-1);
    }
    if ( std::abs(Fx[ getFacex_cell_WithHalo(ix,iy,iz)])>1e10 ) {
      std::cerr << "error in routine " + callingRoutine + " in Fx[" << ix << "," << iy << "," << iz << "]" << std::endl;
      exit(-1);
    }
  }

  for (int ix=0; ix<numberOfCellsPerAxisX+2; ix++)
  for (int iy=0; iy<numberOfCellsPerAxisY+3; iy++)
  for (int iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
    if ( std::abs(uy[ getFacey_cell_WithHalo(ix,iy,iz)])>1e10 ) {
        std::cerr << "error in routine " + callingRoutine + " in uy[" << ix << "," << iy << "," << iz << "]" << std::endl;
        exit(-1);
    }
    if ( std::abs(Fy[ getFacey_cell_WithHalo(ix,iy,iz)])>1e10 ) {
        std::cerr << "error in routine " + callingRoutine + " in Fy[" << ix << "," << iy << "," << iz << "]" << std::endl;
        exit(-1);
    }
  }

  for (int ix=0; ix<numberOfCellsPerAxisX+2; ix++)
  for (int iy=0; iy<numberOfCellsPerAxisY+2; iy++)
  for (int iz=0; iz<numberOfCellsPerAxisZ+3; iz++) {
    if ( std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz)])>1e10 ) {
      std::cerr << "error in in routine " + callingRoutine + " uz[" << ix << "," << iy << "," << iz << "]: " << uz[ getFacez_cell_WithHalo(ix,iy,iz)] << std::endl;
      exit(-1);
    }
    if ( std::abs(Fz[ getFacez_cell_WithHalo(ix,iy,iz)])>1e10 ) {
      std::cerr << "error in in routine " + callingRoutine + " Fz[" << ix << "," << iy << "," << iz << "]" << std::endl;
      exit(-1);
    }
  }
  #endif
}


/**
 * Plot a vtk file. This function probably never has to be changed when you do
 * your assessment.
 */
void plotVTKFile() {
  static int vtkFileCounter = 0;

  std::ostringstream  outputFileName;
  outputFileName << "res/output-" << vtkFileCounter << ".vtk";
  std::ofstream out;
  out.open( outputFileName.str().c_str() );

  std::cout << "\t write " << outputFileName.str();

  out << "# vtk DataFile Version 2.0" << std::endl
      << "Tobias Weinzierl: CPU, Manycore and Cluster Computing" << std::endl
      << "ASCII" << std::endl << std::endl;

  out << "DATASET STRUCTURED_POINTS" << std::endl
      << "DIMENSIONS  "
        << numberOfCellsPerAxisX+1 << " "
        << numberOfCellsPerAxisY+1 << " "
        << numberOfCellsPerAxisZ+1
        << std::endl << std::endl;
  out << "ORIGIN 0 0 0 " << std::endl << std::endl;
  out << "SPACING "
        << getH() << " "
        << getH() << " "
        << getH() << " "
        << std::endl << std::endl;

  const int numberOfVertices = (numberOfCellsPerAxisX+1) * (numberOfCellsPerAxisY+1) * (numberOfCellsPerAxisZ+1);
  out << "POINT_DATA " << numberOfVertices << std::endl << std::endl;


  out << "VECTORS velocity float" << std::endl;
  for (int iz=1; iz<numberOfCellsPerAxisZ+2; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2; ix++) {
        out << 0.25 * ( ux[ getFacex_cell_WithHalo(ix,iy-1,iz-1) ] + ux[ getFacex_cell_WithHalo(ix,iy-0,iz-1) ] + ux[ getFacex_cell_WithHalo(ix,iy-1,iz-0) ] + ux[ getFacex_cell_WithHalo(ix,iy,iz) ] ) << " ";
        out << 0.25 * ( uy[ getFacey_cell_WithHalo(ix-1,iy,iz-1) ] + uy[ getFacey_cell_WithHalo(ix-0,iy,iz-1) ] + uy[ getFacey_cell_WithHalo(ix-1,iy,iz-0) ] + uy[ getFacey_cell_WithHalo(ix,iy,iz) ] ) << " ";
        out << 0.25 * ( uz[ getFacez_cell_WithHalo(ix-1,iy-1,iz) ] + uz[ getFacez_cell_WithHalo(ix-0,iy-1,iz) ] + uz[ getFacez_cell_WithHalo(ix-1,iy-0,iz) ] + uz[ getFacez_cell_WithHalo(ix,iy,iz) ] ) << std::endl;
      }
    }
  }
  out << std::endl << std::endl;

  //
  // For debugging, it sometimes pays off to see F. Usually not required - notably not for this year's
  // assignments
  //
  #ifdef CheckVariableValues
  out << "VECTORS F float" << std::endl;
  for (int iz=1; iz<numberOfCellsPerAxisZ+2; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2; ix++) {
        out << 0.25 * ( Fx[ getFacex_cell_WithHalo(ix,iy-1,iz-1) ] + Fx[ getFacex_cell_WithHalo(ix,iy-0,iz-1) ] + Fx[ getFacex_cell_WithHalo(ix,iy-1,iz-0) ] + Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] ) << " ";
        out << 0.25 * ( Fy[ getFacey_cell_WithHalo(ix-1,iy,iz-1) ] + Fy[ getFacey_cell_WithHalo(ix-0,iy,iz-1) ] + Fy[ getFacey_cell_WithHalo(ix-1,iy,iz-0) ] + Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] ) << " ";
        out << 0.25 * ( Fz[ getFacez_cell_WithHalo(ix-1,iy-1,iz) ] + Fz[ getFacez_cell_WithHalo(ix-0,iy-1,iz) ] + Fz[ getFacez_cell_WithHalo(ix-1,iy-0,iz) ] + Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] ) << std::endl;
      }
    }
  }
  out << std::endl << std::endl;
  #endif

  out << "SCALARS ink float 1" << std::endl;
  out << "LOOKUP_TABLE default" << std::endl;
  for (int iz=0; iz<numberOfCellsPerAxisZ+1; iz++) {
    for (int iy=0; iy<numberOfCellsPerAxisY+1; iy++) {
      for (int ix=0; ix<numberOfCellsPerAxisX+1; ix++) {
        out << ink[ getVertexIndex(ix,iy,iz) ] << std::endl;
      }
    }
  }
  out << std::endl << std::endl;

  const int numberOfCells = numberOfCellsPerAxisX * numberOfCellsPerAxisY * numberOfCellsPerAxisZ;
  out << "CELL_DATA " << numberOfCells << std::endl << std::endl;

  out << "SCALARS pressure float 1" << std::endl;
  out << "LOOKUP_TABLE default" << std::endl;

  for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+1; ix++) {
        out << p[ getCellIndex(ix,iy,iz) ] << std::endl;
      }
    }
  }

  out << "SCALARS obstacle float 1" << std::endl;
  out << "LOOKUP_TABLE default" << std::endl;

  for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+1; ix++) {
        out << cellIsInside[ getCellIndex(ix,iy,iz) ] << std::endl;
      }
    }
  }

  //
  // For debugging, it sometimes pays off to see the rhs of the pressure
  // Poisson equation. Usually not required - notably not for this year's
  // assignments
  //
  #ifdef CheckVariableValues
  out << "SCALARS rhs float 1" << std::endl;
  out << "LOOKUP_TABLE default" << std::endl;

  for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+1; ix++) {
        out << rhs[ getCellIndex(ix,iy,iz) ] << std::endl;
      }
    }
  }
  #endif

  out.close();

  vtkFileCounter++;
}


/**
 * Computes a helper velocity. See book of Griebel for details.
 */
void computeF() {
  const double alpha = timeStepSize / getH();

  for (int iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
      for (int ix=2; ix<numberOfCellsPerAxisX+3-2; ix++) {
        if (
          cellIsInside[getCellIndex(ix-1,iy,iz)]
          &&
          cellIsInside[getCellIndex(ix,iy,iz)]
        ) {
          const double diffusiveTerm =
            + (-1.0 * ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ] + 2.0 * ux[ getFacex_cell_WithHalo(ix,iy,iz) ] - 1.0 * ux[ getFacex_cell_WithHalo(ix+1,iy,iz) ] )
            + (-1.0 * ux[ getFacex_cell_WithHalo(ix,iy-1,iz) ] + 2.0 * ux[ getFacex_cell_WithHalo(ix,iy,iz) ] - 1.0 * ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ] )
            + (-1.0 * ux[ getFacex_cell_WithHalo(ix,iy,iz-1) ] + 2.0 * ux[ getFacex_cell_WithHalo(ix,iy,iz) ] - 1.0 * ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ] );

          const double convectiveTerm =
            + ( (ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix+1,iy,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix+1,iy,iz) ]) - (ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz) ])    *(ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz) ]) )
            + ( (uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ]) - (uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix+1,iy-1,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy-1,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz) ]) )
            + ( (uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ]) - (uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix+1,iy,iz-1) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz-1) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz) ]) )
            + alpha * ( std::abs(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix+1,iy,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]-ux[ getFacex_cell_WithHalo(ix+1,iy,iz) ]) - std::abs(ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz) ])    *(ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]-ux[ getFacex_cell_WithHalo(ix,iy,iz) ]) )
            + alpha * ( std::abs(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]-ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ]) - std::abs(uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix+1,iy-1,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy-1,iz) ]-ux[ getFacex_cell_WithHalo(ix,iy,iz) ]) )
            + alpha * ( std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]-ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ]) - std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix+1,iy,iz-1) ])*(ux[ getFacex_cell_WithHalo(ix,iy,iz-1) ]-ux[ getFacex_cell_WithHalo(ix,iy,iz) ]) )
            ;

          Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] =
           ux[ getFacex_cell_WithHalo(ix,iy,iz) ]
           - timeStepSize/ReynoldsNumber * 1.0/getH()/getH() * diffusiveTerm
           - timeStepSize * 1.0/getH()/4.0 * convectiveTerm;
        }
      }
    }
  }

  for (int iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
    for (int iy=2; iy<numberOfCellsPerAxisY+3-2; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2-1; ix++) {
        if (
          cellIsInside[getCellIndex(ix,iy-1,iz)]
          &&
          cellIsInside[getCellIndex(ix,iy,iz)]
        ) {
          const double diffusiveTerm =
           + (-1.0 * uy[ getFacey_cell_WithHalo(ix-1,iy,iz) ] + 2.0 * uy[ getFacey_cell_WithHalo(ix,iy,iz) ] - 1.0 * uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ] )
           + (-1.0 * uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ] + 2.0 * uy[ getFacey_cell_WithHalo(ix,iy,iz) ] - 1.0 * uy[ getFacey_cell_WithHalo(ix,iy+1,iz) ] )
           + (-1.0 * uy[ getFacey_cell_WithHalo(ix,iy,iz-1) ] + 2.0 * uy[ getFacey_cell_WithHalo(ix,iy,iz) ] - 1.0 * uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ] )
           ;

          const double convectiveTerm =
           + ( (ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ])*(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ]) - (ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix-1,iy+1,iz) ]) *(uy[ getFacey_cell_WithHalo(ix-1,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz) ]) )
           + ( (uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy+1,iz) ])*(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy+1,iz) ]) - (uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz) ])     *(uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz) ]) )
           + ( (uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ])*(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ]) - (uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix,iy+1,iz-1) ]) *(uy[ getFacey_cell_WithHalo(ix,iy,iz-1) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz) ]) )
           + alpha * ( std::abs(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ])*(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]-uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ]) - std::abs(ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix-1,iy+1,iz) ]) *(uy[ getFacey_cell_WithHalo(ix-1,iy,iz) ]-uy[ getFacey_cell_WithHalo(ix,iy,iz) ]) )
           + alpha * ( std::abs(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy+1,iz) ])*(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]-uy[ getFacey_cell_WithHalo(ix,iy+1,iz) ]) - std::abs(uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz) ])     *(uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]-uy[ getFacey_cell_WithHalo(ix,iy,iz) ]) )
           + alpha * ( std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ])*(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]-uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ]) - std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix,iy+1,iz-1) ]) *(uy[ getFacey_cell_WithHalo(ix,iy,iz-1) ]-uy[ getFacey_cell_WithHalo(ix,iy,iz) ]) )
           ;

          Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] =
           uy[ getFacey_cell_WithHalo(ix,iy,iz) ]
           - timeStepSize/ReynoldsNumber * 1.0/getH()/getH() * diffusiveTerm
           - timeStepSize * 1.0/getH()/4.0 * convectiveTerm;
        }
      }
    }
  }

  for (int iz=2; iz<numberOfCellsPerAxisZ+3-2; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2-1; ix++) {
        if (
          cellIsInside[getCellIndex(ix,iy,iz-1)]
          &&
          cellIsInside[getCellIndex(ix,iy,iz)]
        ) {
          const double diffusiveTerm =
           + (-1.0 * uz[ getFacez_cell_WithHalo(ix-1,iy,iz) ] + 2.0 * uz[ getFacez_cell_WithHalo(ix,iy,iz) ] - 1.0 * uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ] )
           + (-1.0 * uz[ getFacez_cell_WithHalo(ix,iy-1,iz) ] + 2.0 * uz[ getFacez_cell_WithHalo(ix,iy,iz) ] - 1.0 * uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ] )
           + (-1.0 * uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ] + 2.0 * uz[ getFacez_cell_WithHalo(ix,iy,iz) ] - 1.0 * uz[ getFacez_cell_WithHalo(ix,iy,iz+1) ] )
           ;

          const double convectiveTerm =
           + ( (ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ])*(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ]) - (ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix-1,iy,iz+1) ]) *(uz[ getFacez_cell_WithHalo(ix-1,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz) ]) )
           + ( (uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ])*(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ]) - (uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy-1,iz+1) ]) *(uz[ getFacez_cell_WithHalo(ix,iy-1,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz) ]) )
           + ( (uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz+1) ])*(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz+1) ]) - (uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz) ])     *(uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz) ]) )
           + alpha * ( std::abs(ux[ getFacex_cell_WithHalo(ix,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ])*(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]-uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ]) - std::abs(ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ]+ux[ getFacex_cell_WithHalo(ix-1,iy,iz+1) ]) *(uz[ getFacez_cell_WithHalo(ix-1,iy,iz) ]-uz[ getFacez_cell_WithHalo(ix,iy,iz) ]) )
           + alpha * ( std::abs(uy[ getFacey_cell_WithHalo(ix,iy,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ])*(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]-uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ]) - std::abs(uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ]+uy[ getFacey_cell_WithHalo(ix,iy-1,iz+1) ]) *(uz[ getFacez_cell_WithHalo(ix,iy-1,iz) ]-uz[ getFacez_cell_WithHalo(ix,iy,iz) ]) )
           + alpha * ( std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz+1) ])*(uz[ getFacez_cell_WithHalo(ix,iy,iz) ]-uz[ getFacez_cell_WithHalo(ix,iy,iz+1) ]) - std::abs(uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]+uz[ getFacez_cell_WithHalo(ix,iy,iz) ])     *(uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ]-uz[ getFacez_cell_WithHalo(ix,iy,iz) ]) )
           ;

          Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] =
           uz[ getFacez_cell_WithHalo(ix,iy,iz) ]
           - timeStepSize/ReynoldsNumber * 1.0/getH()/getH() * diffusiveTerm
           - timeStepSize * 1.0/getH()/4.0 * convectiveTerm;
        }
      }
    }
  }

  validateThatEntriesAreBounded( "computeF" );
}


/**
 * Compute the right-hand side. This basically means how much a flow would
 * violate the incompressibility if there were no pressure.
 */
void computeRhs() {
  for (int iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2-1; ix++) {
        if ( cellIsInside[getCellIndex(ix,iy,iz)] ) {
          rhs[ getCellIndex(ix,iy,iz) ] = 1.0/timeStepSize/getH()*
            (
              Fx[getFacex_cell_WithHalo(ix+1,iy,iz)] - Fx[getFacex_cell_WithHalo(ix,iy,iz)] +
              Fy[getFacey_cell_WithHalo(ix,iy+1,iz)] - Fy[getFacey_cell_WithHalo(ix,iy,iz)] +
              Fz[getFacez_cell_WithHalo(ix,iy,iz+1)] - Fz[getFacez_cell_WithHalo(ix,iy,iz)]
            );
        }
      }
    }
  }
}


/**
 * Set boundary conditions for pressure. The values of the pressure at the
 * domain boundary might depend on the pressure itself. So if we update it
 * in the algorithm, we afterwards have to reset the boundary conditions
 * again.
 */
void setPressureBoundaryConditions() {
  int ix, iy, iz;

  // Neumann Boundary conditions for p
  for (iy=0; iy<numberOfCellsPerAxisY+2; iy++) {
    for (iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
      ix=0;
      p[ getCellIndex(ix,iy,iz) ]   = p[ getCellIndex(ix+1,iy,iz) ];
      ix=numberOfCellsPerAxisX+1;
      p[ getCellIndex(ix,iy,iz) ]   = p[ getCellIndex(ix-1,iy,iz) ];
    }
  }
  for (ix=0; ix<numberOfCellsPerAxisX+2; ix++) {
    for (iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
      iy=0;
      p[ getCellIndex(ix,iy,iz) ]   = p[ getCellIndex(ix,iy+1,iz) ];
      iy=numberOfCellsPerAxisY+1;
      p[ getCellIndex(ix,iy,iz) ]   = p[ getCellIndex(ix,iy-1,iz) ];
    }
  }
  for (ix=0; ix<numberOfCellsPerAxisX+2; ix++) {
    for (iy=0; iy<numberOfCellsPerAxisY+2; iy++) {
      iz=0;
      p[ getCellIndex(ix,iy,iz) ]   = p[ getCellIndex(ix,iy,iz+1) ];
      iz=numberOfCellsPerAxisZ+1;
      p[ getCellIndex(ix,iy,iz) ]   = p[ getCellIndex(ix,iy,iz-1) ];
    }
  }

  // Normalise pressure at rhs to zero
  for (iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
    for (iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
      p[ getCellIndex(numberOfCellsPerAxisX+1,iy,iz) ]   = 0.0;
    }
  }

  // Pressure conditions around obstacle
  for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
      for (int ix=2; ix<numberOfCellsPerAxisX+1; ix++) {
        if (cellIsInside[getCellIndex(ix,iy,iz)]) {
          if ( !cellIsInside[getCellIndex(ix-1,iy,iz)] ) { // left neighbour
            p[getCellIndex(ix-1,iy,iz)]     = p[getCellIndex(ix,iy,iz)];
          }
          if ( !cellIsInside[getCellIndex(ix+1,iy,iz)] ) { // right neighbour
            p[getCellIndex(ix+1,iy,iz)]     = p[getCellIndex(ix,iy,iz)];
          }
          if ( !cellIsInside[getCellIndex(ix,iy-1,iz)] ) { // bottom neighbour
            p[getCellIndex(ix,iy-1,iz)]     = p[getCellIndex(ix,iy,iz)];
          }
          if ( !cellIsInside[getCellIndex(ix,iy+1,iz)] ) { // right neighbour
            p[getCellIndex(ix,iy+1,iz)]     = p[getCellIndex(ix,iy,iz)];
          }
          if ( !cellIsInside[getCellIndex(ix,iy,iz-1)] ) { // front neighbour
            p[getCellIndex(ix,iy,iz-1)]     = p[getCellIndex(ix,iy,iz)];
          }
          if ( !cellIsInside[getCellIndex(ix,iy,iz+1)] ) { // right neighbour
            p[getCellIndex(ix,iy,iz+1)]     = p[getCellIndex(ix,iy,iz)];
          }
        }
      }
    }
  }
}


/**
 * Determine the new pressure. The operation depends on properly set pressure
 * boundary conditions. See setPressureBoundaryConditions().
 *
 * @return Number of iterations required or max number plus one if we had to
 *         stop iterating as the solver diverged.
 * EDITED MY ME
 */
int computeP() {
  double       globalResidual         = 1.0;
  double       firstResidual          = 1.0;
  double       previousGlobalResidual = 2.0;
  int          iterations             = 0;

  while(
   (
    std::abs(globalResidual-previousGlobalResidual)>PPESolverThreshold
    &&
    iterations<MaxComputePIterations
    &&
    std::abs(globalResidual)>PPESolverThreshold
    &&
    (globalResidual/firstResidual>PPESolverThreshold)
   )
   ||
   (iterations%2==1) // we have alternating omega, so we allow only even iteration counts
  ) {
    const double omega = iterations%2==0 ? 1.2 : 0.8;
    setPressureBoundaryConditions();

    previousGlobalResidual = globalResidual;
    globalResidual         = 0.0;

    // EDIT STARTED HERE

    /*
    for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
      for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
        for (int ix=1; ix<numberOfCellsPerAxisX+1; ix++) {
          if ( cellIsInside[getCellIndex(ix,iy,iz)] ) { // THIS IS THE STEP 2 BIT
            double residual = rhs[ getCellIndex(ix,iy,iz) ] +
              1.0/getH()/getH()*
              (
                - 1.0 * p[ getCellIndex(ix-1,iy,iz) ]
                - 1.0 * p[ getCellIndex(ix+1,iy,iz) ]
                - 1.0 * p[ getCellIndex(ix,iy-1,iz) ]
                - 1.0 * p[ getCellIndex(ix,iy+1,iz) ]
                - 1.0 * p[ getCellIndex(ix,iy,iz-1) ]
                - 1.0 * p[ getCellIndex(ix,iy,iz+1) ]
                + 6.0 * p[ getCellIndex(ix,iy,iz) ]
              );
            globalResidual              += residual * residual;
            p[ getCellIndex(ix,iy,iz) ] += -omega * residual / 6.0 * getH() * getH();
          }
        }
      }
    }

    for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
      for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
        for (int ix=1; ix<numberOfCellsPerAxisX+1; ix++) {
          if (iz <= numberOfCellsPerAxisZ and iy <= numberOfCellsPerAxisY and ix <= numberOfCellsPerAxisX){
            printf("----\n");
            printf("%i\n", iz);
            printf("%i\n", iy);
            printf("%i\n", ix);
          }
        }
      }
    }
    */
    
    const int numberOfBlocksX = (numberOfCellsPerAxisX/BLOCKDIMENTION);
    const int numberOfBlocksY = (numberOfCellsPerAxisY/BLOCKDIMENTION);
    const int numberOfBlocksZ = (numberOfCellsPerAxisZ/BLOCKDIMENTION);

    /*
    printf("-------------\n");
    printf("%i\n", numberOfBlocksX);
    printf("%i\n", numberOfBlocksY);
    printf("%i\n", numberOfBlocksZ);
    */

    for (int x_Block_NoHalo = 0; x_Block_NoHalo < numberOfBlocksX; ++x_Block_NoHalo){
      for (int y_Block_NoHalo = 0; y_Block_NoHalo < numberOfBlocksY; ++y_Block_NoHalo){
        for (int z_Block_NoHalo = 0; z_Block_NoHalo < numberOfBlocksZ; ++z_Block_NoHalo){
          // iterate through every block

          // !objectInBlock[getBlockIndexFromBlockCoordinates(x_Block_NoHalo,y_Block_NoHalo,z_Block_NoHalo)]

          if (false){
            // this means that it is just fluid meaning that it can be SIMD

            for (int iz = 0; iz < BLOCKDIMENTION; ++iz){
              for (int iy = 0; iy < BLOCKDIMENTION; ++iy){
                for (int ix = 0; ix < BLOCKDIMENTION; ++ix){
                  // iterate through all indexes in the block
                  int x_cell_WithHalo = x_Block_NoHalo * BLOCKDIMENTION + ix;
                  int y_cell_WithHalo = y_Block_NoHalo * BLOCKDIMENTION + iy;
                  int z_cell_WithHalo = z_Block_NoHalo * BLOCKDIMENTION + iz;
                  // the index of current

                  if ( cellIsInside[getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo)] ) {
                    double residual = rhs[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo) ] +
                      1.0/getH()/getH()*
                      (
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo-1,y_cell_WithHalo,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo+1,y_cell_WithHalo,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo-1,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo+1,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo-1) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo+1) ]
                        + 6.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo) ]
                    );
                    globalResidual              += residual * residual;
                    p[ getCellIndex(ix,iy,iz) ] += -omega * residual / 6.0 * getH() * getH();
                  }
                }
              }
            }
          }
          else{
            // this means that it is a mixture of fluid and object this cannot be SIMD

            for (int ix = 0; ix < BLOCKDIMENTION; ++ix){
              for (int iy = 0; iy < BLOCKDIMENTION; ++iy){
                for (int iz = 0; iz < BLOCKDIMENTION; ++iz){
                  // iterate through all indexes in the block
                  int x_cell_WithHalo = (x_Block_NoHalo) * BLOCKDIMENTION + ix + 1;
                  int y_cell_WithHalo = (y_Block_NoHalo) * BLOCKDIMENTION + iy + 1;
                  int z_cell_WithHalo = (z_Block_NoHalo) * BLOCKDIMENTION + iz + 1;
                  // the index of current
                  /*
                  if (z_cell_WithHalo <= numberOfCellsPerAxisZ and y_cell_WithHalo <= numberOfCellsPerAxisY and x_cell_WithHalo <= numberOfCellsPerAxisX)
                  {
                    printf("------\n");
                    printf("%i\n", z_cell_WithHalo);
                    printf("%i\n", y_cell_WithHalo);
                    printf("%i\n", x_cell_WithHalo);
                    printf("------\n");
                  }
                  */

                  if ( cellIsInside[getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo)] ) {
                    double residual = rhs[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo) ] +
                      1.0/getH()/getH()*
                      (
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo-1,y_cell_WithHalo,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo+1,y_cell_WithHalo,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo-1,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo+1,z_cell_WithHalo) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo-1) ]
                        - 1.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo+1) ]
                        + 6.0 * p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo) ]
                    );
                    globalResidual              += residual * residual;
                    p[ getCellIndex(x_cell_WithHalo,y_cell_WithHalo,z_cell_WithHalo) ] += -omega * residual / 6.0 * getH() * getH();
                  }
                }
              }
            }
          }
          // here is where I would add the skip any blocks inside the object
        }
      }
    }

    // EDIT ENDED HERE

    globalResidual        = std::sqrt(globalResidual);
    firstResidual         = firstResidual==0 ? globalResidual : firstResidual;
    iterations++;

  }

  std::cout << "iterations n=" << iterations
            << ", |res(n)|_2=" << globalResidual
            << ", |res(n-1)|_2=" << previousGlobalResidual
            << ", |res(n-1)|_2-|res(n)|_2=" << (previousGlobalResidual-globalResidual);

  return iterations;
}


/**
 * @todo Your job if you attend the Scientific Computing submodule. Otherwise empty.
 */
void updateInk() {
}



/**
 * Once we have F and a valid pressure p, we may update the velocities.
 */
void setNewVelocities() {
  for (int iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
      for (int ix=2; ix<numberOfCellsPerAxisX+3-2; ix++) {
        ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] - timeStepSize/getH() * ( p[getCellIndex(ix,iy,iz)] - p[getCellIndex(ix-1,iy,iz)]);
      }
    }
  }

  for (int iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
    for (int iy=2; iy<numberOfCellsPerAxisY+3-2; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2-1; ix++) {
        uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] - timeStepSize/getH() * ( p[getCellIndex(ix,iy,iz)] - p[getCellIndex(ix,iy-1,iz)]);
      }
    }
  }

  for (int iz=2; iz<numberOfCellsPerAxisZ+3-2; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+2-1; ix++) {
        uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] - timeStepSize/getH() * ( p[getCellIndex(ix,iy,iz)] - p[getCellIndex(ix,iy,iz-1)]);
      }
    }
  }
}


/**
 * Setup our scenario, i.e. initialise all the big arrays and set the
 * right boundary conditions. This is something you might want to change in
 * part three of the assessment.
 * EDITED BY ME
 */
void setupScenario() {
  const int numberOfCells = (numberOfCellsPerAxisX+2) * (numberOfCellsPerAxisY+2) * (numberOfCellsPerAxisZ+2);
  const int numberOfBlocks = (numberOfCellsPerAxisX/BLOCKDIMENTION) * (numberOfCellsPerAxisY/BLOCKDIMENTION) * (numberOfCellsPerAxisZ/BLOCKDIMENTION);
  const int numberOfBlocksX = (numberOfCellsPerAxisX/BLOCKDIMENTION);
  const int numberOfBlocksY = (numberOfCellsPerAxisY/BLOCKDIMENTION);
  const int numberOfBlocksZ = (numberOfCellsPerAxisZ/BLOCKDIMENTION);

  const int numberOfFacesX = (numberOfCellsPerAxisX+3) * (numberOfCellsPerAxisY+2) * (numberOfCellsPerAxisZ+2);
  const int numberOfFacesY = (numberOfCellsPerAxisX+2) * (numberOfCellsPerAxisY+3) * (numberOfCellsPerAxisZ+2);
  const int numberOfFacesZ = (numberOfCellsPerAxisX+2) * (numberOfCellsPerAxisY+2) * (numberOfCellsPerAxisZ+3);

  ux  = 0;
  uy  = 0;
  uz  = 0;
  Fx  = 0;
  Fy  = 0;
  Fz  = 0;

  p   = 0;
  rhs = 0;
  ink = 0;

  ux  = new (std::nothrow) double[numberOfFacesX];
  uy  = new (std::nothrow) double[numberOfFacesY];
  uz  = new (std::nothrow) double[numberOfFacesZ];
  Fx  = new (std::nothrow) double[numberOfFacesX];
  Fy  = new (std::nothrow) double[numberOfFacesY];
  Fz  = new (std::nothrow) double[numberOfFacesZ];

  objectInBlock = new (std::nothrow) bool[numberOfBlocksX];

  p   = new (std::nothrow) double[numberOfCells];
  rhs = new (std::nothrow) double[numberOfCells];

  ink = new (std::nothrow) double[(numberOfCellsPerAxisX+1) * (numberOfCellsPerAxisY+1) * (numberOfCellsPerAxisZ+1)];

  cellIsInside = new (std::nothrow) bool[numberOfCells];

  if (
    ux  == 0 ||
    uy  == 0 ||
    uz  == 0 ||
    Fx  == 0 ||
    Fy  == 0 ||
    Fz  == 0 ||
    p   == 0 ||
    rhs == 0
  ) {
    std::cerr << "could not allocate memory. Perhaps not enough memory free?" << std::endl;
    exit(-1);
  }


  for (int i=0; i<(numberOfCellsPerAxisX+1) * (numberOfCellsPerAxisY+1) * (numberOfCellsPerAxisZ+1); i++) {
    ink[i]          = 0.0;
  }

  for (int i=0; i<numberOfCells; i++) {
    p[i]            = 0.0;
    cellIsInside[i] = true;
  }

  for (int i=0; i<numberOfFacesX; i++) {
    ux[i]=0;
    Fx[i]=0;
  }
  for (int i=0; i<numberOfFacesY; i++) {
    uy[i]=0;
    Fy[i]=0;
  }
  for (int i=0; i<numberOfFacesZ; i++) {
    uz[i]=0;
    Fz[i]=0;
  }

  for (int i=0; i<numberOfBlocks; i++) {
    objectInBlock[i]           = false;
  }

  //
  // Insert the obstacle that forces the fluid to do something interesting.
  //
  int sizeOfObstacle    = numberOfCellsPerAxisY/3;
  int xOffsetOfObstacle = sizeOfObstacle*2;
  if (sizeOfObstacle<2) sizeOfObstacle = 2;
  int zDelta = numberOfCellsPerAxisZ<=8 ? 0 : sizeOfObstacle/3;
  for (int iz=1 + zDelta; iz<numberOfCellsPerAxisZ+2-zDelta; iz++) {
    cellIsInside[ getCellIndex(xOffsetOfObstacle,    sizeOfObstacle+1,iz) ] = false;
    cellIsInside[ getCellIndex(xOffsetOfObstacle+1,  sizeOfObstacle+1,iz) ] = false;
    for (int ii=0; ii<sizeOfObstacle; ii++) {
      cellIsInside[ getCellIndex(xOffsetOfObstacle+ii,  sizeOfObstacle+ii+2,iz) ] = false;
      cellIsInside[ getCellIndex(xOffsetOfObstacle+ii+1,sizeOfObstacle+ii+2,iz) ] = false;
      cellIsInside[ getCellIndex(xOffsetOfObstacle+ii+2,sizeOfObstacle+ii+2,iz) ] = false;
    }
    cellIsInside[ getCellIndex(xOffsetOfObstacle+sizeOfObstacle+0,  2*sizeOfObstacle+2,iz) ] = false;
    cellIsInside[ getCellIndex(xOffsetOfObstacle+sizeOfObstacle+1,  2*sizeOfObstacle+2,iz) ] = false;
  }

  int blockX;
  int blockY;
  int blockZ;

  for (int iz=1; iz<numberOfCellsPerAxisZ+1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+1; iy++) {
      for (int ix=1; ix<numberOfCellsPerAxisX+1; ix++) {
        // iterate through every cell (they are padded)
        if (cellIsInside[ getCellIndex(ix, iy, iz) ]){
          // if the cell is in the object
          objectInBlock[getBlockIndexFromCellCoordinates(ix,iy,iz)] = true;
          // get the corresponding block
          // change the blocks index to true in the objectInBlock array
        }
      }
    }
  }
/*
  for (int i = 0; i < numberOfBlocks; ++i){
    printf("%d\n", objectInBlock[i]);
  }
*/
  validateThatEntriesAreBounded("setupScenario()");
}


/**
 * Clean up the system
 * EDITED BY ME
 */
void freeDataStructures() {
  delete[] p;
  delete[] ink;
  delete[] rhs;
  delete[] ux;
  delete[] uy;
  delete[] uz;
  delete[] Fx;
  delete[] Fy;
  delete[] Fz;
  delete[] cellIsInside;

  ux  = 0;
  uy  = 0;
  uz  = 0;
  Fx  = 0;
  Fy  = 0;
  Fz  = 0;
  p   = 0;
  rhs = 0;
  ink = 0;
}


/**
 * - Handle all the velocities at the domain boundaries. We either
 *   realise no-slip or free-slip.
 *
 * - Set the inflow and outflow profile.
 *
 * - Fix all the F values. Along the boundary, the F values equal the
 *   velocity values.
 *
 */
void setVelocityBoundaryConditions(double time) {
  int ix, iy, iz;

  validateThatEntriesAreBounded("setVelocityBoundaryConditions(double)[in]");

  const bool UseNoSlip = true;

  // We ensure that no fluid leaves the domain. For this, we set the velocities
  // along the boundary to zero. Furthermore,  we ensure that the tangential
  // components of all velocities are zero. For this, we set the ghost/virtual
  // velocities to minus the original one. If we interpolate linearly in-between,
  // we obtain zero tangential speed.
  for (iy=0; iy<numberOfCellsPerAxisY+2; iy++) {
    for (iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
      ix=0;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = 0.0;
      ix=1;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = 0.0;

      ix=0;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ];
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ];

      ix=numberOfCellsPerAxisX+2;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = 0.0;
      ix=numberOfCellsPerAxisX+1;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = 0.0;

      ix=numberOfCellsPerAxisX+1;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uy[ getFacey_cell_WithHalo(ix-1,iy,iz) ];
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uz[ getFacez_cell_WithHalo(ix-1,iy,iz) ];
    }
  }

  for (ix=0; ix<numberOfCellsPerAxisX+2; ix++) {
    for (iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
      iy=0;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = 0.0;
      iy=1;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = 0.0;

      iy=0;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ];
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ];

      iy=numberOfCellsPerAxisY+2;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = 0.0;
      iy=numberOfCellsPerAxisY+1;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = 0.0;

      iy=numberOfCellsPerAxisY+1;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * ux[ getFacex_cell_WithHalo(ix,iy-1,iz) ];
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uz[ getFacez_cell_WithHalo(ix,iy-1,iz) ];
    }
  }

  for (ix=0; ix<numberOfCellsPerAxisX+2; ix++) {
    for (iy=0; iy<numberOfCellsPerAxisY+2; iy++) {
      iz=0;
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = 0.0;
      iz=1;
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = 0.0;

      iz=0;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ];
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ];

      iz=numberOfCellsPerAxisZ+2;
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = 0.0;
      iz=numberOfCellsPerAxisZ+1;
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = 0.0;

      iz=numberOfCellsPerAxisZ+1;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * ux[ getFacex_cell_WithHalo(ix,iy,iz-1) ];
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = (UseNoSlip ? -1.0 : 1.0) * uy[ getFacey_cell_WithHalo(ix,iy,iz-1) ];
    }
  }

  validateThatEntriesAreBounded("setVelocityBoundaryConditions(double)[no slip set]");


  // Don't switch on in-flow immediately but slowly induce it into the system
  //
  const double inputProfileScaling = std::min(time*10.0,1.0);
  // const double inputProfileScaling = 1.0;

  for (iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
    for (iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
      const double yDistance = (iy-1) * (1.0/numberOfCellsPerAxisY);
      const double zDistance = (iz-1) * (1.0/numberOfCellsPerAxisZ);
      const double inflow    = UseNoSlip ? inputProfileScaling * yDistance * (1.0-yDistance) * zDistance * (1.0-zDistance) : inputProfileScaling;

      ix=0;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = inflow;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = 0.0;
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = 0.0;
      ix=1;
      ux[ getFacex_cell_WithHalo(ix,iy,iz) ] = inflow;
      uy[ getFacey_cell_WithHalo(ix,iy,iz) ] = 0.0;
      uz[ getFacez_cell_WithHalo(ix,iy,iz) ] = 0.0;

      // outflow
      ux[ getFacex_cell_WithHalo(numberOfCellsPerAxisX+2,iy,iz) ] = ux[ getFacex_cell_WithHalo(numberOfCellsPerAxisX,iy,iz) ];
      ux[ getFacex_cell_WithHalo(numberOfCellsPerAxisX+1,iy,iz) ] = ux[ getFacex_cell_WithHalo(numberOfCellsPerAxisX,iy,iz) ];
      uy[ getFacey_cell_WithHalo(numberOfCellsPerAxisX+1,iy,iz) ] = uy[ getFacey_cell_WithHalo(numberOfCellsPerAxisX+0,iy,iz) ];
      uz[ getFacez_cell_WithHalo(numberOfCellsPerAxisX+1,iy,iz) ] = uz[ getFacez_cell_WithHalo(numberOfCellsPerAxisX+0,iy,iz) ];
    }
  }

  validateThatEntriesAreBounded("setVelocityBoundaryConditions(double)[inflow set]");

  //
  //  Once all velocity boundary conditions are set, me can fix the F values
  //
  for (iy=0; iy<numberOfCellsPerAxisY+2; iy++) {
    for (iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
      ix=0;
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];

      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fy[ getFaceIndex^(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];

      ix=1;
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];

      ix=numberOfCellsPerAxisX+1;
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uz[ getFacey_cell_WithHalo(ix,iy,iz) ];
      //Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uy[ getFacez_cell_WithHalo(ix,iy,iz) ];

      ix=numberOfCellsPerAxisX+2;
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
    }
  }

  for (ix=0; ix<numberOfCellsPerAxisX+2; ix++) {
    for (iz=0; iz<numberOfCellsPerAxisZ+2; iz++) {
      iy=0;
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];

      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      
      iy=1;
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ]
                                          ;
      iy=numberOfCellsPerAxisY+1;
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      //Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];

      iy=numberOfCellsPerAxisY+2;
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
    }
  }


  for (ix=0; ix<numberOfCellsPerAxisX+2; ix++) {
    for (iy=0; iy<numberOfCellsPerAxisY+2; iy++) {
      iz=0;
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      //Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      //Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      
      iz=1;
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];

      iz=numberOfCellsPerAxisZ+1;
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
      Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
      Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      //Fx[ getFacex_cell_WithHalo(ix,iy,iz) ] = uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
      //Fy[ getFacey_cell_WithHalo(ix,iy,iz) ] = ux[ getFacex_cell_WithHalo(ix,iy,iz) ];

      iz=numberOfCellsPerAxisZ+2;
      Fz[ getFacez_cell_WithHalo(ix,iy,iz) ] = uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
    }
  }

  //
  // Handle the obstacle
  //
  for (int iz=1; iz<numberOfCellsPerAxisZ+2-1; iz++) {
    for (int iy=1; iy<numberOfCellsPerAxisY+2-1; iy++) {
      for (int ix=2; ix<numberOfCellsPerAxisX+3-2; ix++) {
        if (cellIsInside[getCellIndex(ix,iy,iz)]) {
          if ( !cellIsInside[getCellIndex(ix-1,iy,iz)] ) { // left neighbour
            ux[ getFacex_cell_WithHalo(ix,iy,iz) ]   = 0.0;
            ux[ getFacex_cell_WithHalo(ix-1,iy,iz) ] = 0.0;
            Fx[ getFacex_cell_WithHalo(ix,iy,iz) ]   = 0.0;
            uy[ getFacey_cell_WithHalo(ix-1,iy,iz) ] = -uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
            uz[ getFacez_cell_WithHalo(ix-1,iy,iz) ] = -uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
          }
          if ( !cellIsInside[getCellIndex(ix+1,iy,iz)] ) { // right neighbour
            ux[ getFacex_cell_WithHalo(ix+1,iy,iz) ] = 0.0;
            ux[ getFacex_cell_WithHalo(ix+2,iy,iz) ] = 0.0;
            Fx[ getFacex_cell_WithHalo(ix+1,iy,iz) ] = 0.0;
            uy[ getFacey_cell_WithHalo(ix+1,iy,iz) ] = -uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
            uz[ getFacez_cell_WithHalo(ix+1,iy,iz) ] = -uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
          }
          if ( !cellIsInside[getCellIndex(ix,iy-1,iz)] ) { // bottom neighbour
            uy[ getFacey_cell_WithHalo(ix,iy,iz) ]   = 0.0;
            uy[ getFacey_cell_WithHalo(ix,iy-1,iz) ] = 0.0;
            Fy[ getFacey_cell_WithHalo(ix,iy,iz) ]   = 0.0;
            ux[ getFacex_cell_WithHalo(ix,iy-1,iz) ] = -ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
            uz[ getFacez_cell_WithHalo(ix,iy-1,iz) ] = -uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
          }
          if ( !cellIsInside[getCellIndex(ix,iy+1,iz)] ) { // top neighbour
            uy[ getFacey_cell_WithHalo(ix,iy+1,iz) ] = 0.0;
            uy[ getFacey_cell_WithHalo(ix,iy+2,iz) ] = 0.0;
            Fy[ getFacey_cell_WithHalo(ix,iy+1,iz) ] = 0.0;
            ux[ getFacex_cell_WithHalo(ix,iy+1,iz) ] = -ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
            uz[ getFacez_cell_WithHalo(ix,iy+1,iz) ] = -uz[ getFacez_cell_WithHalo(ix,iy,iz) ];
          }
          if ( !cellIsInside[getCellIndex(ix,iy,iz-1)] ) { // front neighbour
            uz[ getFacez_cell_WithHalo(ix,iy,iz) ]   = 0.0;
            uz[ getFacez_cell_WithHalo(ix,iy,iz-1) ] = 0.0;
            Fz[ getFacez_cell_WithHalo(ix,iy,iz) ]   = 0.0;
            ux[ getFacex_cell_WithHalo(ix,iy,iz-1) ] = -ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
            uy[ getFacey_cell_WithHalo(ix,iy,iz-1) ] = -uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
          }
          if ( !cellIsInside[getCellIndex(ix,iy,iz+1)] ) { // right neighbour
            uz[ getFacez_cell_WithHalo(ix,iy,iz+1) ] = 0.0;
            uz[ getFacez_cell_WithHalo(ix,iy,iz+2) ] = 0.0;
            Fz[ getFacez_cell_WithHalo(ix,iy,iz+1) ] = 0.0;
            ux[ getFacex_cell_WithHalo(ix,iy,iz+1) ] = -ux[ getFacex_cell_WithHalo(ix,iy,iz) ];
            uy[ getFacey_cell_WithHalo(ix,iy,iz+1) ] = -uy[ getFacey_cell_WithHalo(ix,iy,iz) ];
          }
        }
      }
    }
  }

  validateThatEntriesAreBounded("setVelocityBoundaryConditions(double)[out]");
}



int main (int argc, char *argv[]) {
  if (argc!=4) {
      std::cout << "Usage: executable number-of-elements-per-axis time-steps-between-plots reynolds-number" << std::endl;
      std::cout << "    number-of-elements-per-axis  Resolution. Start with 20, but try to increase as much as possible later." << std::endl;
      std::cout << "    time-between-plots           Determines how many files are written. Set to 0 to switch off plotting (for performance studies)." << std::endl;
      std::cout << "    reynolds-number              Use something in-between 1 and 1000. Determines viscosity of fluid." << std::endl;
      return 1;
  }

  numberOfCellsPerAxisY    = atoi(argv[1]);
  numberOfCellsPerAxisZ    = numberOfCellsPerAxisY / 2;
  numberOfCellsPerAxisX    = numberOfCellsPerAxisY * 5;
  double timeBetweenPlots  = atof(argv[2]);
  ReynoldsNumber           = atof(argv[3]);

  std::cout << "Re=" << ReynoldsNumber << std::endl;

  std::cout << "create " << numberOfCellsPerAxisX << "x" << numberOfCellsPerAxisY << "x" << numberOfCellsPerAxisZ << " grid" << std::endl;
  setupScenario();

  //   dt <= C Re dx^2
  // whereas the CFD lab at TUM uses
  //   const double MaximumTimeStepSize                 = 0.8 * std::min( ReynoldsNumber/2.0/(3.0/numberOfCellsPerAxisY/numberOfCellsPerAxisY), 1.0/numberOfCellsPerAxisY );
  const double TimeStepSizeConstant = 1e-4;
  const double MaximumTimeStepSize  = TimeStepSizeConstant * ReynoldsNumber / numberOfCellsPerAxisY / numberOfCellsPerAxisY;
  const double MinimalTimeStepSize  = MaximumTimeStepSize / 800;

  timeStepSize = MaximumTimeStepSize;
  std::cout << "start with time step size " << timeStepSize << std::endl;

  setVelocityBoundaryConditions(0.0);
  std::cout << "velocity start conditions are set";
  if (timeBetweenPlots>0.0) {
    plotVTKFile();
  }
  std::cout << std::endl;

  double t = 0.0;
  double tOfLastSnapshot                       = 0.0;
  int    timeStepCounter                       = 0;
  int    numberOfTimeStepsWithOnlyOneIteration = 0;
  while (t<0.01) {
    std::cout << "time step " << timeStepCounter << ": t=" << t << "\t dt=" << timeStepSize << "\t";

    setVelocityBoundaryConditions(t);
    computeF();
    computeRhs();
    // 1/100th of the startup phase is tackled with overrelaxation; afterwards,
    // we use underrelaxation as we are interested in the pressure gradient, i.e.
    // we want to have a smooth solution
    int innerIterationsRequired = computeP();
    setNewVelocities();
    updateInk();

    if (timeBetweenPlots>0.0 && (t-tOfLastSnapshot>timeBetweenPlots)) {
      plotVTKFile();
      tOfLastSnapshot = t;
    }

    if (innerIterationsRequired>=MinAverageComputePIterations) {
      numberOfTimeStepsWithOnlyOneIteration--;
    }
    else if (innerIterationsRequired<=std::max(MinAverageComputePIterations/10,2) ) {
      numberOfTimeStepsWithOnlyOneIteration++;
    }
    else {
      numberOfTimeStepsWithOnlyOneIteration /= 2;
    }

    if (numberOfTimeStepsWithOnlyOneIteration>IterationsBeforeTimeStepSizeIsAltered && timeStepSize < MaximumTimeStepSize) {
      timeStepSize *= (1.0+ChangeOfTimeStepSize);
      numberOfTimeStepsWithOnlyOneIteration = 0;
      std::cout << "\t time step size seems to be too small. Increased to " << timeStepSize << " to speed up simulation";
    }
    else if (numberOfTimeStepsWithOnlyOneIteration<-IterationsBeforeTimeStepSizeIsAltered && timeStepSize>MinimalTimeStepSize) {
      timeStepSize /= 2.0;
      numberOfTimeStepsWithOnlyOneIteration = 0;
      std::cout << "\t time step size seems to be too big. Reduced to " << timeStepSize << " to keep simulation stable";
    }

    t += timeStepSize;
    timeStepCounter++;

    std::cout << std::endl;
  }

  std::cout << "free data structures" << std::endl;
  freeDataStructures();

  return 0;
}