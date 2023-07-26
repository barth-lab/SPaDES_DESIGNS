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
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_protocols_hydrate_hydrate_remove_functions_HH
#define INCLUDED_protocols_hydrate_hydrate_remove_functions_HH

// Protocols

// Core
#include <core/pose/Pose.fwd.hh>

// Basic

// Utility headers

//wym

namespace protocols {
namespace hydrate {

//// Function to remove all de novo waters from pose at end
void
remove_de_novo_waters(
core::pose::Pose & pose
);
//
//// Function to remove far away waters from pose at end
void
remove_far_away_waters(
core::pose::Pose & pose
);
	
} // close hydrate
} // close hydrate

#endif
