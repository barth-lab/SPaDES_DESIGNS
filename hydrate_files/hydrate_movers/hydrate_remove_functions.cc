// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/protocols/hydrate/Hydrate.cc
/// @brief The Hydrate Protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai, Lucas S. P. Rudden

// Protocols
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>

// Core
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/ResfileReader.hh>

// Basic
#include <basic/Tracer.hh>

// Utility
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// Numeric
#include <numeric/NumericTraits.hh>

// C++
#include <string>
#include <iterator>

// wym
#include <protocols/hydrate/hydrate_functions.hh>

#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask

static basic::Tracer TR( "protocols.hydrate.hydrate_remove_functions" );

namespace protocols {
namespace hydrate {

using namespace core;

//// Function to remove all de novo waters from pose at end
void
remove_de_novo_waters(
 core::pose::Pose & pose
) {
 for ( core::Size resid = pose.total_residue(); resid > 1; --resid ) {
  if ( pose.residue( resid ).name() == "TP3" ) {
   protocols::grafting::simple_movers::DeleteRegionMover deleter( resid, resid );
   deleter.apply( pose );
   TR << "Removed de novo water " << resid << std::endl;
  }
 }
}
//
//// Function to remove far away waters from pose at end
void
remove_far_away_waters(
 core::pose::Pose & pose
) {
 for ( core::Size resid = pose.total_residue(); resid > 1; --resid ) {
  if ( pose.residue( resid ).name() == "TP3" ) {
   core::Real xcoord = pose.residue( resid ).xyz( 1 ).x();
   core::Real ycoord = pose.residue( resid ).xyz( 1 ).y();
   core::Real zcoord = pose.residue( resid ).xyz( 1 ).z();
   if ( xcoord > 9999 && ycoord > 9999 && zcoord > 9999 ) {
    protocols::grafting::simple_movers::DeleteRegionMover deleter( resid, resid );
    deleter.apply( pose );
    TR << "Removed far away water " << resid << std::endl;
   }
  }
 }
}


} // close hydrate
} // close hydrate
