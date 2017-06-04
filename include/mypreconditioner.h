/*
 *    An MIT forward solver code based on the deal.II (www.dealii.org) library.
 *    Copyright (C) 2013-2015 Ross Kynch & Paul Ledger, Swansea Unversity.
 * 
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.    
*/

#ifndef MYPRECONDITIONER_H
#define MYPRECONDITIONER_H

// deal.II includes:
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/lac/sparse_decomposition.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

// std includes:

// My includes:
#include <all_data.h>

using namespace dealii;

//******************************************************Precondotioner*************************************************
namespace Preconditioner
{
  //_______________________EddyCurrentPreconditionerBase________________________________________________________-
  // base class for all eddy current preconditioners:
  class EddyCurrentPreconditionerBase : public Subscriptor //Clase de Deal.ii
  {
  public:
  //-----------------------------EddycurrentPreconditionerBase()------------------------------------------------------
    EddyCurrentPreconditionerBase (const BlockSparseMatrix<double> &A); //BlockSparseMatrix<> clase de Deal.ii
  //-----------------------------vmult()------------------------------------------------------------------------------
    virtual void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<> clase de Deal.ii
  //-----------------------------Tvmult()------------------------------------------------------------------------------------
    virtual void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<> clase de Deal.ii
  protected:
  //-----------------------------precon_matrix-------------------------------------------------------------------------------
    const SmartPointer<const BlockSparseMatrix<double> > precon_matrix; //BlockSparseMatrix<> clase de Deal.ii
  };
  //_______________________EddyCurrentPreconditionerBase,end________________________________________________________-
  
  //_______________________EddyCurrentPreconditioner_1x1_lowOrder________________________________________________________-
  // Preconditioner for iterative solver, made up of a single block.
  // This essentially assumes the use of p=0.
  class EddyCurrentPreconditioner_1x1_lowOrder : public virtual EddyCurrentPreconditionerBase
  {
  public:
  //-----------------------------EddycurrentPreconditioner_1x1_lowOrder()------------------------------------------------------
    EddyCurrentPreconditioner_1x1_lowOrder (const BlockSparseMatrix<double> &A);
  //-----------------------------vmult()------------------------------------------------------------------------------
    virtual void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  //-----------------------------Tvmult()------------------------------------------------------------------------------------
    virtual void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<> Clase de Deal.ii
  private: 
  //-----------------------------block0----------------------------------------------------------------------------------
        SparseDirectUMFPACK block0; //SparseDirectUMFPACK Clase de Deal.ii
  };
  //_______________________EddyCurrentPreconditioner_1x1_lowOrder,end________________________________________________________-


  //_______________________EddyCurrentPreconditioner_2x2_lowHighOrder________________________________________________________-
  // Preconditioning for use with an FE_Nedelec-based 2x2 block matrix
  // made up of:
  // block 0: lowest
  // block 1: higher order edges/faces/cells.
  class EddyCurrentPreconditioner_2x2_lowHighOrder : public EddyCurrentPreconditionerBase
  {
  public:
  //-----------------------------EddycurrentPreconditioner_2x2_lowHighOrder()------------------------------------------------------
    EddyCurrentPreconditioner_2x2_lowHighOrder (const BlockSparseMatrix<double> &A); //BlockSparseMatrix<>, Clase de Deal.ii
  //-----------------------------vmult()------------------------------------------------------------------------------
    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  //-----------------------------Tvmult()------------------------------------------------------------------------------------
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  private:
  //-----------------------------block0----------------------------------------------------------------------------------
  // Block 0: lowest order edges.    
    SparseDirectUMFPACK block0; //SparseDirectUMFPACK Clase de Deal.ii
  //-----------------------------block1----------------------------------------------------------------------------------
  // Block 1: higher order edges/faces/cells (or just edges if Block2/3 enabled:
    SparseILU<double> block1; //SparseILU Clase de Deal.ii
  //-----------------------------block1_data----------------------------------------------------------------------------------
    SparseILU<double>::AdditionalData block1_data; //SparseILU<double>::AdditionalData Clase de Deal.ii
  };
  //_______________________EddyCurrentPreconditioner_2x2_lowHighOrder,end________________________________________________________-


  //_______________________EddyCurrentPreconditioner_3x3_lowHighOrderGradients________________________________________________________
  // Preconditioning for use with an FE_Nedelec-based 3x3 block matrix
  // made up of:
  // block 0: lowest
  // block 1: gradient-based higher order edges/faces/cells
  // block 2: non-gradient-based higher order edges/faces/cells
  class EddyCurrentPreconditioner_3x3_lowHighOrderGradients : public EddyCurrentPreconditionerBase
  {
  public:
  //-----------------------------EddycurrentPreconditioner_3x3_lowHighOrderGradients()------------------------------------------------------
    EddyCurrentPreconditioner_3x3_lowHighOrderGradients (const BlockSparseMatrix<double> &A); //BlockSparseMatrix<>, Clase de Deal.ii
  //-----------------------------vmult()------------------------------------------------------------------------------
    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  //-----------------------------Tvmult()------------------------------------------------------------------------------------
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  private:
  //-----------------------------block0----------------------------------------------------------------------------------
  // Block 0: lowest order:
    SparseDirectUMFPACK block0; //SparseDirectUMFPACK Clase de Deal.ii
  //-----------------------------block1----------------------------------------------------------------------------------
  // Block 1: higher order gradients:
    SparseDirectUMFPACK block1; //SparseDirectUMFPACK Clase de Deal.ii
  //-----------------------------block2----------------------------------------------------------------------------------
  // Block 2: higher order non-gradients:
    SparseDirectUMFPACK block2; //SparseDirectUMFPACK Clase de Deal.ii
  };
  //_______________________EddyCurrentPreconditioner_3x3_lowHighOrderGradients,end________________________________________________________


  //_______________________EddyCurrentPreconditioner_5x5________________________________________________________-
  // Preconditioning for use with an FE_Nedelec-based 3x3 block matrix
  // made up of:
  // block 0: lowest
  // block 1: gradient-based higher order edges/faces/cells
  // block 2: non-gradient-based higher order edges/faces/cells
  class EddyCurrentPreconditioner_5x5 : public EddyCurrentPreconditionerBase
  {
  public:
  //-----------------------------EddycurrentPreconditioner_3x3_lowHighOrderGradients()------------------------------------------------------
    EddyCurrentPreconditioner_5x5 (const BlockSparseMatrix<double> &A); //BlockSparseMatrix<>, Clase de Deal.ii
  //-----------------------------vmult()------------------------------------------------------------------------------
    void vmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  //-----------------------------Tvmult()------------------------------------------------------------------------------------
    void Tvmult (BlockVector<double> &dst, const BlockVector<double> &src) const; //BlockVector<>, Clase de Deal.ii
  private:
  //-----------------------------block0----------------------------------------------------------------------------------
  // Block 0: lowest order:
    SparseDirectUMFPACK block0; //SparseDirectUMFPACK Clase de Deal.ii
  //-----------------------------block1----------------------------------------------------------------------------------
  // Block 1: real gradients:
    SparseDirectUMFPACK block1; //SparseDirectUMFPACK Clase de Deal.ii
  //-----------------------------block2----------------------------------------------------------------------------------
  // Block 2: real non-gradients:
    SparseDirectUMFPACK block2; //SparseDirectUMFPACK Clase de Deal.ii
    
  // NOTE: Don't actually need block3 or block4
  // They will match block1 and block2, so can just reuse them.
  // Block 3: imag gradients:
  // SparseDirectUMFPACK block3;
  // Block 4: imag non-gradients:
  // SparseDirectUMFPACK block4;
  };
  //_______________________EddyCurrentPreconditioner_5x5,end________________________________________________________-

}
#endif
