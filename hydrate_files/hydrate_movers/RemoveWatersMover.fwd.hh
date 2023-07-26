// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/protocols/hydrate/RemoveWatersMover.cc
/// @brief The RemoveWaters mover protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_protocols_hydrate_RemoveWatersMover_fwd_hh
#define INCLUDED_protocols_hydrate_RemoveWatersMover_fwd_hh

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace hydrate {

/// @brief
class RemoveWatersMover;

typedef utility::pointer::shared_ptr< RemoveWatersMover > RemoveWatersMoverOP;
typedef utility::pointer::shared_ptr< RemoveWatersMover const > RemoveWatersMoverCOP;

} // namespace hydrate
} // namespace protocols

#endif  // INCLUDED_protocols_hydrate_RemoveWatersMover_FWD_HH

