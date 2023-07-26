// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/hydrate/Hydrate.cc
/// @brief A protocol to add explicit water molecules as part of Hydrate/SPaDES.
/// @detailed
/// @author Joaquin Ambia Garrido, Jason K. Lai, Lucas S. P. Rudden
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitute.org) to permit multi-threaded packing.

// Protocols
#include <protocols/hydrate/Hydrate.hh>
#include <protocols/hydrate/hydrate_functions.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Core
#include <core/pack/pack_rotamers.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#ifdef MULTI_THREADED
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>
#endif

// Utility headers

// Numeric headers

// C++ headers
#include <iostream>
#include <string>

#include <core/kinematics/MoveMap.hh> // AUTO IWYU For MoveMap
#include <core/pack/task/TaskFactory.hh> // AUTO IWYU For TaskFactory
#include <core/pack/task/operation/TaskOperations.hh> // AUTO IWYU For InitializeFromCommandline, ReadResfile, Rest...
#include <core/scoring/ScoreFunctionFactory.hh> // AUTO IWYU For get_score_function
#include <core/scoring/constraints/util.hh> // AUTO IWYU For add_constraints_from_cmdline_to_scorefxn
#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask

//wym

// Construct tracer.
static basic::Tracer TR( "protocols.hydrate.Hydrate" );

namespace protocols {
namespace hydrate {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////

// Default constructor
Hydrate::Hydrate(): Mover(),
	score_fxn_(new core::scoring::ScoreFunction()),
	main_task_factory_(new core::pack::task::TaskFactory())
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	TR << "Initializing hydrate protocol" << std::endl;

	// Setting up task
	main_task_factory_->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( utility::pointer::make_shared< operation::ReadResfile >() );
	} else {
		operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
		main_task_factory_->push_back( rtrop );
	}

	// Setting up scorefxn and constraints
	score_fxn_ = scoring::get_score_function();
	if ( option[ in::file::centroid_input ].value() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	}

	// Setting up options
	if ( option[ OptionKeys::packing::resfile ].user() ) {
		protein_flexibility_ = "resfile";
	} else {
		protein_flexibility_ = option[ protein_flexibility ]();
	}
	hydrate_all_ = option[ hydrate_all ]();
	partial_hydrate_dew_ = option[ partial_hydrate_dew ]();
	short_residence_time_mode_ = option[short_residence_time_mode ]();
	near_water_threshold_ = option[ near_water_threshold ]();  //yumeng
	minimize_bb_where_packing_ = option[ minimize_bb_where_packing ]();  // yumeng

#ifdef MULTI_THREADED
	set_interaction_graph_threads( option[ multithreading::interaction_graph_threads ]() );
#endif

	type("Hydrate");
}


// Constructor that allows to specify a score function and the protein flexibility (region to pack & minimize )
Hydrate::Hydrate(
	scoring::ScoreFunctionOP scorefxn,
	std::string protein_flexibility
): Mover(),
	main_task_factory_(new core::pack::task::TaskFactory())
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	TR << "Initializing hydrate protocol" << std::endl;

	// Setting up task
	main_task_factory_->push_back( utility::pointer::make_shared< operation::InitializeFromCommandline >() );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory_->push_back( utility::pointer::make_shared< operation::ReadResfile >() );
	} else {
		operation::RestrictToRepackingOP rtrop( new operation::RestrictToRepacking );
		main_task_factory_->push_back( rtrop );
	}

	// Setting up scorefxn and constraints
	score_fxn_ = scorefxn;

	if ( option[ in::file::centroid_input ].value() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
	}

	// Setting up options
	protein_flexibility_ = protein_flexibility;
	hydrate_all_ = option[ hydrate_all ]();
	partial_hydrate_dew_ = option[ partial_hydrate_dew ]();
	short_residence_time_mode_ = option[short_residence_time_mode ]();
	near_water_threshold_ = option[ near_water_threshold ]();  //yumeng
	minimize_bb_where_packing_ = option[ minimize_bb_where_packing ]();  // yumeng

	type("Hydrate");
}

// Copy constructor
Hydrate::Hydrate(Hydrate const & hyd): Mover(hyd)
{
	copy_data(*this, hyd);
}

// Assignment operator
Hydrate &
Hydrate::operator=(Hydrate const & hyd)
{
	// Abort self-assignment.
	if ( this == &hyd ) {
		return *this;
	}

	Mover::operator=(hyd);
	copy_data(*this, hyd);
	return *this;
}

// Destructor
Hydrate::~Hydrate() = default;


// Mover methods
std::string
Hydrate::get_name() const
{
	return type();
}

protocols::moves::MoverOP
Hydrate::clone() const
{
	return utility::pointer::make_shared< Hydrate >( *this );
}

protocols::moves::MoverOP
Hydrate::fresh_instance() const
{
	return utility::pointer::make_shared< Hydrate >();
}


/// @details
/// @param    <pose>: the structure to be moved
/// @remarks
void
Hydrate::apply(Pose & pose)
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::hydrate;

	TR << "Applying " << get_name() << std::endl;

	if ( option[ only_remove_non_buried_waters ]() ) {
		remove_non_buried_wat(pose);
		return;
	}

	if ( option[ just_score ]() ) {
		water_specific_hbond_energy( pose, *score_fxn_ ); //wym output water specific hbond energy to pdb
		add_water_overcoordinated_hb_score( pose, *score_fxn_ );
		show_water_hb_network( pose );
		return;
	}

	set_water_info_and_add_de_novo_water( pose, *score_fxn_ );
	if ( hydrate_all_ ) {
		(*score_fxn_)(pose);    // needs to be scored so that the pose has an energy object
		hydrate_cavities( pose );
	}

	if ( option[ force_enforce_all_waters ]() ) {
		enforce_all_waters( pose );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Packing (Double edge water rotamers)
	//
	// Packing, during this first round, de novo waters generate rotamers with two optimized
	// hydrogen bonds (double edge (dew)).
	// Some de novo water molecules will not be included in this dew packing round, they will
	// be included in the sew round.

	if ( short_residence_time_mode_  ) partial_hydrate_dew_ = 0; // If running on short residence time dont add water here
	if ( partial_hydrate_dew_ < 1 ) set_dew_waters_not_to_be_included(pose, partial_hydrate_dew_);
#ifdef MULTI_THREADED
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose, interaction_graph_threads_ ));
#else
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
#endif
	task->initialize_from_command_line(); // -ex1 -ex2  etc.
	kinematics::MoveMap mm;
	set_task_and_movemap( pose, protein_flexibility_, task, mm, minimize_bb_where_packing_ );

	if ( !short_residence_time_mode_ ) { // If running on short residence time dont add water here
		core::Size nloop( option[ pack_nloop ].value() );
		pack::pack_rotamers_loop( pose, *score_fxn_, task, nloop ); // previous default 25
                remove_high_energy_water_molecules( pose, *score_fxn_ );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Single Edge (anchor) Water (sew) rotamers packing
	//
	// Any de novo water molecule that could not find an energetically favorable position near the protein
	// in the previous packing, now attempts to find it by constructing rotamers with a single
	// optimized hydrogen bond (single edge (sew)) to its anchor.
	// This does not imply that the rotamers will have only that one hb.
	if ( short_residence_time_mode_ ) enforce_all_waters( pose );
	// If running in short residence time mode, all waters packed in this step (Single Edge) while enforced to
	// stay active, near the protein, and then they are packed again without being enforced,
	core::pack::palette::CustomBaseTypePackerPaletteOP palette( utility::pointer::make_shared< core::pack::palette::CustomBaseTypePackerPalette >() );
	palette->add_type("TP3");
#ifdef MULTI_THREADED
	pack::task::PackerTaskOP sew_task( pack::task::TaskFactory::create_packer_task( pose, palette, interaction_graph_threads_ ));
#else
	pack::task::PackerTaskOP sew_task( pack::task::TaskFactory::create_packer_task( pose, palette ));
#endif
	sew_task->initialize_from_command_line(); // -ex1 -ex2  etc.
	(*score_fxn_)(pose);
	get_ready_for_sew_packing( pose, sew_task);
	pack::pack_rotamers_loop( pose, *score_fxn_, sew_task, 25);

	if ( !short_residence_time_mode_ ) remove_high_energy_water_molecules( pose, *score_fxn_);

	if ( short_residence_time_mode_ ) {
		// If running in short residence time mode, at this point, water molecules are packed concurrently with
		// the rest of the protein considered flexible, using input water rotamers (oxygen position dependent) and
		// not being enforced to stay near the protein.
		remove_all_anchors_and_ENF( pose );
#ifdef MULTI_THREADED
		pack::task::PackerTaskOP task_strm( pack::task::TaskFactory::create_packer_task( pose, interaction_graph_threads_ ));
#else
		pack::task::PackerTaskOP task_strm( pack::task::TaskFactory::create_packer_task( pose ));
#endif
		task_strm->initialize_from_command_line(); // -ex1 -ex2  etc.
		kinematics::MoveMap mm_strm;
		set_task_and_movemap (pose, protein_flexibility_, task_strm, mm_strm, minimize_bb_where_packing_);
		pack::pack_rotamers_loop( pose, *score_fxn_, task_strm, 25);
		remove_high_energy_water_molecules( pose, *score_fxn_);
	}

	if ( option[ show_rotamer_count ]() ) {
		return;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Minimizing
	//
	if ( option[ min_backbone_file ].user() ) {
		set_bb_movemap( pose, option[ min_backbone_file ](), mm);
	}
	if ( option[ show_derivatives_check ]() ) {
		TR << "Minimizing and showing derivatives" << std::endl;
		optimization::AtomTreeMinimizer().run( pose, mm, *score_fxn_, optimization::MinimizerOptions("dfpmin",0.001,true,
			true, true));
	} else {
		TR << "Minimizing" << std::endl;
		optimization::AtomTreeMinimizer().run( pose, mm, *score_fxn_, optimization::MinimizerOptions("dfpmin",0.001,true));
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Removing non buried water molecules
	//
	remove_non_buried_wat(pose);

	if ( option[ show_residues_near_water ]() ) {
		print_residues_near_water( pose );
	}

	/////////////////////////////////////////////////////////////////////////
	//
	//  Adding over coordinated hydrogen bonds scores
	//  And if the hybrid score function is used, it is considered here
	//  [[This must be the last action of the apply function]]
	//
	water_specific_hbond_energy( pose, *score_fxn_ ); //wym output water specific hbond energy to pdb
	add_water_overcoordinated_hb_score( pose, *score_fxn_ );
	show_water_hb_network( pose );

}

#ifdef MULTI_THREADED
/// @brief Set the number of threads to use for interaction graph precomputation during packing.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
void
Hydrate::set_interaction_graph_threads(
	core::Size const setting
) {
	interaction_graph_threads_ = setting;
}

/// @brief Get the number of threads to use for interaction graph precomputation during packing.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
core::Size
Hydrate::interaction_graph_threads() const {
	return interaction_graph_threads_;
}
#endif


void
Hydrate::copy_data(
	Hydrate hyd_to,
	Hydrate hyd_from)
{
	hyd_to.score_fxn_ = hyd_from.score_fxn_;
	hyd_to.main_task_factory_ = hyd_from.main_task_factory_;
#ifdef MULTI_THREADED
	hyd_to.interaction_graph_threads_ = hyd_from.interaction_graph_threads_;
#endif
}

}  // namespace hydrate
}  // namespace protocols

