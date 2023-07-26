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
/// @author Joaquin Ambia Garrido, Jason K. Lai

// unit headers
#include <protocols/hydrate/HydratePackerMover.hh>
#include <protocols/hydrate/HydratePackerMoverCreator.hh>

// Protocols
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/hydrate/hydrate_functions.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Core
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/task/PackerTask.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Numeric headers

// C++ headers
#include <iostream>
#include <string>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/io/raw_data/ScoreMap.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// LSPR
#include <protocols/filters/filter_schemas.hh>
//for resfile reading
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/kinematics/MoveMap.hh>

// Construct tracer.
static basic::Tracer TR( "protocols.hydrate.HydratePackerMover" );

namespace protocols {
namespace hydrate {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////

// Default constructor
HydratePackerMover::HydratePackerMover(): Mover( "HydratePackerMover" ),
	scorefxn_( ),
	water_hybrid_sf_( true ),
	resfile_( ),
	hyfile_( ),
	hydrate_all_( false ),
	nofasol_( ),
	weights_( ),
	pretalaris2013_( false ),
	hbond_threshold_( -0.1 ),
	water_rotamers_cap_( 500 ),
	extra_res_fa_( ),
	extra_res_cen_( ),
	spanfile_( ),
	membed_init_( false ),
	mhbond_depth_( false ),
	cst_weight_( 1.0 ),
	cst_fa_weight_( 1.0 )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	TR << "Initializing hydrate protocol" << std::endl;

	// Setting up scorefxn and constraints
	if ( scorefxn_ == nullptr ) {
		scorefxn_ = scoring::get_score_function();
	}
	if ( option[ in::file::centroid_input ].user() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	}
	
	// Setting up options
	if ( option[ OptionKeys::packing::resfile ].user() || resfile_ != "" ) {
		protein_flexibility_ = "resfile"; }
	hydrate_all_ = option[ hydrate_all ]();
	partial_hydrate_dew_ = option[ partial_hydrate_dew ]();
	short_residence_time_mode_ = option[short_residence_time_mode ]();
	near_water_threshold_ = option[ near_water_threshold ]();  //yumeng
	minimize_bb_where_packing_ = option[ minimize_bb_where_packing ]();  // yumeng

	type("Hydrate");
}


HydratePackerMover::HydratePackerMover(
	scoring::ScoreFunctionOP scorefxn,
	std::string protein_flexibility
): Mover(),
	water_hybrid_sf_( true ),
	resfile_( ),
	hyfile_( ),
	hydrate_all_( false ),
	nofasol_( ),
	weights_( ),
	pretalaris2013_( false ),
	hbond_threshold_( -0.1 ),
	water_rotamers_cap_( 500 ),
	extra_res_fa_( ),
	extra_res_cen_( ),
	spanfile_( ),
	membed_init_( false ),
	mhbond_depth_( false ),
	cst_weight_( 1.0 ),
	cst_fa_weight_( 1.0 ),
	protein_flexibility_( )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	TR << "Initializing hydrate protocol" << std::endl;

	// Setting up scorefxn and constraints
	if ( scorefxn_ == nullptr ) {
		scorefxn_ = scorefxn;
	}
	if ( option[ in::file::centroid_input ].user() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*scorefxn_);
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
HydratePackerMover::HydratePackerMover(HydratePackerMover const & hyd): protocols::moves::Mover( hyd )
{
	scorefxn_ = hyd.scorefxn_;
	water_hybrid_sf_ = hyd.water_hybrid_sf_;
	resfile_ = hyd.resfile_;
	hyfile_ = hyd.hyfile_;
	hydrate_all_ = hyd.hydrate_all_;
	nofasol_ = hyd.nofasol_;
	weights_ = hyd.weights_;
	pretalaris2013_ = hyd.pretalaris2013_;
	hbond_threshold_ = hyd.hbond_threshold_;
	water_rotamers_cap_ = hyd.water_rotamers_cap_;
	extra_res_fa_ = hyd.extra_res_fa_;
	extra_res_cen_ = hyd.extra_res_cen_;
	spanfile_ = hyd.spanfile_;
	membed_init_ = hyd.membed_init_;
	mhbond_depth_ = hyd.mhbond_depth_;
	cst_weight_ = hyd.cst_weight_;
	cst_fa_weight_ = hyd.cst_fa_weight_;
	//copy_data(*this, hyd);
}

// Assignment operator
HydratePackerMover &
HydratePackerMover::operator=(HydratePackerMover const & hyd)
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
HydratePackerMover::~HydratePackerMover() {}

// Mover methods
std::string
HydratePackerMover::get_name() const
{
	return type();
}

protocols::moves::MoverOP
HydratePackerMover::clone() const
{
	return HydratePackerMoverOP( new HydratePackerMover( *this ) );
}

protocols::moves::MoverOP
HydratePackerMover::fresh_instance() const
{
	return HydratePackerMoverOP( new HydratePackerMover() );
}


/// @details
/// @param    <pose>: the structure to be moved
/// @remarks
void
HydratePackerMover::apply(Pose & pose)
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace basic::options::OptionKeys::score;
	using namespace pack::task;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Preset options
	if ( water_hybrid_sf_ != true) basic::options::option[ basic::options::OptionKeys::score::water_hybrid_sf ].value( water_hybrid_sf_ );
	if ( resfile_ != "" ) basic::options::option[ basic::options::OptionKeys::packing::resfile ].value( resfile_ );
	if ( hyfile_ != "" ) basic::options::option[ basic::options::OptionKeys::hydrate::hyfile ].value( hyfile_ );
	if ( hydrate_all_ != false ) basic::options::option[ basic::options::OptionKeys::hydrate::hydrate_all ].value( hydrate_all_ );
	if ( nofasol_.size() != 0 ) basic::options::option[ basic::options::OptionKeys::hydrate::ignore_fa_sol_at_positions ].value( nofasol_ );
	if ( weights_ != "" ) basic::options::option[ basic::options::OptionKeys::score::weights ].value( weights_ );
	if ( pretalaris2013_ != false ) basic::options::option[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ].value( pretalaris2013_ );
	if ( hbond_threshold_ != -0.1 ) basic::options::option[ basic::options::OptionKeys::hydrate::hbond_threshold ].value( hbond_threshold_ );
	if ( water_rotamers_cap_ != 500 ) basic::options::option[ basic::options::OptionKeys::hydrate::water_rotamers_cap ].value( water_rotamers_cap_ );
	if ( extra_res_fa_ != "" ) basic::options::option[ basic::options::OptionKeys::in::file::extra_res_fa ].value( extra_res_fa_ );
	if ( extra_res_cen_ != "" ) basic::options::option[ basic::options::OptionKeys::in::file::extra_res_cen ].value( extra_res_cen_ );
	if ( spanfile_ != "" ) basic::options::option[ basic::options::OptionKeys::in::file::spanfile ].value( spanfile_ );
	if ( membed_init_ ) basic::options::option[ basic::options::OptionKeys::membrane::Membed_init ].value( membed_init_ );
	if ( mhbond_depth_ ) basic::options::option[ basic::options::OptionKeys::membrane::Mhbond_depth ].value( mhbond_depth_ );
	if ( cst_weight_ != 1.0 ) basic::options::option[ basic::options::OptionKeys::constraints::cst_weight ].value( cst_weight_ );
	if ( cst_fa_weight_ != 1.0 ) basic::options::option[ basic::options::OptionKeys::constraints::cst_fa_weight ].value( cst_fa_weight_ );

	// Setting up scorefxn and constraints
	if ( scorefxn_ == nullptr ) {
		scorefxn_ = scoring::get_score_function();
	}
	if ( option[ OptionKeys::in::file::centroid_input ].user() ) {
		scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	} else {
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*scorefxn_);
	}

	// Setting up options
	if ( option[ OptionKeys::packing::resfile ].user() || resfile_ != "" ) {
		protein_flexibility_ = "resfile";
	}
	hydrate_all_ = option[ hydrate_all ]();
	partial_hydrate_dew_ = option[ partial_hydrate_dew ]();
	short_residence_time_mode_ = option[ short_residence_time_mode ]();
	near_water_threshold_ = option[ near_water_threshold ]();  //yumeng
	minimize_bb_where_packing_ = option[ minimize_bb_where_packing ]();  // yumeng

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	TR << "Applying " << get_name() << std::endl;

	if ( option[ only_remove_non_buried_waters ]() ) {
		remove_non_buried_wat(pose);
		return;
	}

	if ( option[ just_score ]() ) {
		water_specific_hbond_energy( pose, *scorefxn_ ); //wym output water specific hbond energy to pdb
		add_water_overcoordinated_hb_score( pose, *scorefxn_ );
		show_water_hb_network( pose );
		return;
	}

	set_water_info_and_add_de_novo_water( pose, *scorefxn_ );
	if ( hydrate_all_ ) {
		(*scorefxn_)(pose);    // needs to be scored so that the pose has an energy object
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
	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line( ); // -ex1 -ex2  etc.
	core::kinematics::MoveMap mm;
	if ( protein_flexibility_ == "" ) protein_flexibility_ == "all";
	set_task_and_movemap( pose, protein_flexibility_, task, mm, minimize_bb_where_packing_ );

	if ( !short_residence_time_mode_ ) { // If running on short residence time dont add water here
		core::Size nloop( option[ pack_nloop ].value() );
		pack::pack_rotamers_loop( pose, *scorefxn_, task, nloop ); // previous default 25
		remove_high_energy_water_molecules( pose, *scorefxn_ );
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
	pack::task::PackerTaskOP sew_task( pack::task::TaskFactory::create_packer_task( pose ));
	sew_task->initialize_from_command_line(); // -ex1 -ex2  etc.
	(*scorefxn_)(pose);
	get_ready_for_sew_packing( pose, sew_task);
	pack::pack_rotamers_loop( pose, *scorefxn_, sew_task, 25);

	if ( !short_residence_time_mode_ ) remove_high_energy_water_molecules( pose, *scorefxn_);

	if ( short_residence_time_mode_ ) {
		// If running in short residence time mode, at this point, water molecules are packed concurrently with
		// the rest of the protein considered flexible, using input water rotamers (oxygen position dependent) and
		// not being enforced to stay near the protein.
		remove_all_anchors_and_ENF( pose );
		pack::task::PackerTaskOP task_strm( pack::task::TaskFactory::create_packer_task( pose ));
		task_strm->initialize_from_command_line(); // -ex1 -ex2  etc.
                core::kinematics::MoveMap mm_strm;
		set_task_and_movemap (pose, protein_flexibility_, task_strm, mm_strm, minimize_bb_where_packing_);
		pack::pack_rotamers_loop( pose, *scorefxn_, task_strm, 25);
		remove_high_energy_water_molecules( pose, *scorefxn_);
	}

	if ( option[ show_rotamer_count ]() ) {
		return;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Minimizing
	//
	//  note: rosettascripts version of hydrate/spades does not have any minimization

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
	water_specific_hbond_energy( pose, *scorefxn_ ); //wym output water specific hbond energy to pdb
	add_water_overcoordinated_hb_score( pose, *scorefxn_ );
	show_water_hb_network( pose );

}

void
HydratePackerMover::copy_data(
	HydratePackerMover hyd_to,
	HydratePackerMover hyd_from)
{
	hyd_to.scorefxn_ = hyd_from.scorefxn_;
	hyd_to.water_hybrid_sf_ = hyd_from.water_hybrid_sf_;
	hyd_to.resfile_ = hyd_from.resfile_;
	hyd_to.hyfile_ = hyd_from.hyfile_;
	hyd_to.hydrate_all_ = hyd_from.hydrate_all_;
	hyd_to.nofasol_ = hyd_from.nofasol_;
	hyd_to.weights_ = hyd_from.weights_;
	hyd_to.pretalaris2013_ = hyd_from.pretalaris2013_;
	hyd_to.hbond_threshold_ = hyd_from.hbond_threshold_;
	hyd_to.water_rotamers_cap_ = hyd_from.water_rotamers_cap_;
	hyd_to.extra_res_fa_ = hyd_from.extra_res_fa_;
	hyd_to.extra_res_cen_ = hyd_from.extra_res_cen_;
	hyd_to.spanfile_ = hyd_from.spanfile_;
	hyd_to.membed_init_ = hyd_from.membed_init_;
	hyd_to.mhbond_depth_ = hyd_from.mhbond_depth_;
	hyd_to.cst_weight_ = hyd_from.cst_weight_;
	hyd_to.cst_fa_weight_ = hyd_from.cst_fa_weight_;
}

// setters
void HydratePackerMover::score_function( scoring::ScoreFunctionOP sf )
{
	runtime_assert( sf != nullptr );
	scorefxn_ = sf;
}

void HydratePackerMover::task_factory( TaskFactoryOP tf )
{
	runtime_assert( tf != nullptr );
	//task_factory_ = tf;
	//main_task_factory_ = tf;
}

//
void HydratePackerMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap
) 
{
	if ( tag->hasOption("water_hybrid_sf") ) {
		water_hybrid_sf_ = tag->getOption<bool>("water_hybrid_sf");
	}
	if ( tag->hasOption("resfile") ) {
		resfile_ = tag->getOption<std::string>("resfile");
	}
	if ( tag->hasOption("hyfile") ) {
		hyfile_ = tag->getOption<std::string>("hyfile");
	}
	if ( tag->hasOption("hydrate_all") ) {
		hydrate_all_ = tag->getOption<bool>("hydrate_all");
	}
	if ( tag->hasOption("nofasol") ) {
		nofasol_ = utility::string_split( tag->getOption<std::string>("nofasol"), ',', core::Size() );
	}
	if ( tag->hasOption("weights") ) {
		weights_ = tag->getOption<std::string>("weights");
	}
	if ( tag->hasOption("pretalaris2013") ) {
		pretalaris2013_ = tag->getOption<bool>("pretalaris2013");
	}
	if ( tag->hasOption("hbond_threshold") ) {
		hbond_threshold_ = tag->getOption<Real>("hbond_threshold");
	}
	if ( tag->hasOption("water_rotamers_cap") ) {
		water_rotamers_cap_ = tag->getOption<Real>("water_rotamers_cap");
	}
	if ( tag->hasOption("extra_res_fa") ) {
		extra_res_fa_ = tag->getOption<std::string>("extra_res_fa");
	}
	if ( tag->hasOption("extra_res_cen") ) {
		extra_res_cen_ = tag->getOption<std::string>("extra_res_cen");
	}
	if ( tag->hasOption("spanfile") ) {
		spanfile_ = tag->getOption<std::string>("spanfile");
	}
	if ( tag->hasOption("membed_init") ) {
		membed_init_ = tag->getOption<bool>("membed_init");
	}
	if ( tag->hasOption("mhbond_depth") ) {
		mhbond_depth_ = tag->getOption<bool>("mhbond_depth");
	}
	if ( tag->hasOption("cst_weight") ) {
		cst_weight_ = tag->getOption<Real>("cst_weight");
	}
	if ( tag->hasOption("cst_fa_weight") ) {
		cst_fa_weight_ = tag->getOption<Real>("cst_fa_weight");
	}

	HydratePackerMover::parse_score_function( tag, datamap );
	HydratePackerMover::parse_task_operations( tag, datamap );
}

/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
HydratePackerMover::parse_score_function(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap
)
{
	scoring::ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == nullptr ) return;
	score_function( new_score_function );
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
HydratePackerMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap
)
{
	core::pack::task::TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	task_factory( new_task_factory );
}

void HydratePackerMover::register_options() {
	using namespace basic::options::OptionKeys;
}

std::string HydratePackerMover::mover_name() {
	return "HydratePackerMover";
}


utility::tag::XMLSchemaComplexTypeGeneratorOP
HydratePackerMover::complex_type_generator_for_hydrate_packer_mover( utility::tag::XMLSchemaDefinition & )
{

	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "water_hybrid_sf", xsct_rosetta_bool, "score::water_hybrid_sf" )
		+ XMLSchemaAttribute( "resfile", xs_string, "in::file::resfile" )
		+ XMLSchemaAttribute( "hyfile", xs_string, "hydrate::hyfile" )
		+ XMLSchemaAttribute( "hydrate_all", xsct_rosetta_bool, "hydrate::hydrate_all" )
		+ XMLSchemaAttribute( "nofasol", xsct_int_cslist, "hydrate::ignore_fa_sol_at_positions" )
		+ XMLSchemaAttribute( "weights", xs_string, "score::weights" )
		+ XMLSchemaAttribute( "pretalaris2013", xsct_rosetta_bool, "mistakes::restore_pre_talaris_2013_behavior" )
		+ XMLSchemaAttribute( "hbond_threshold", xsct_real, "hydrate::hbond_threshold" )
		+ XMLSchemaAttribute( "water_rotamers_cap", xsct_real, "hydrate::water_rotamers_cap" )
		+ XMLSchemaAttribute( "extra_res_fa", xs_string, "in::file::extra_res_fa" )
		+ XMLSchemaAttribute( "extra_res_cen", xs_string, "in::file::extra_res_cen" )
		+ XMLSchemaAttribute( "spanfile", xs_string, "in::file::spanfile" )
		+ XMLSchemaAttribute( "membed_init", xsct_rosetta_bool, "membrane::Membed_init" )
		+ XMLSchemaAttribute( "mhbond_depth", xsct_rosetta_bool, "membrane::Mhbond_depth" )
		+ XMLSchemaAttribute( "cst_weight", xsct_real, "constraints::cst_weight" )
		+ XMLSchemaAttribute( "cst_fa_weight", xsct_real, "constraints::cst_fa_weight" );

	rosetta_scripts::attributes_for_parse_score_function( attributes );
	rosetta_scripts::attributes_for_parse_task_operations( attributes );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attributes )
		.add_optional_name_attribute();
	return ct_gen;
}


void HydratePackerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_hydrate_packer_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Hydrates and repacks system as a mover." )
		.write_complex_type_to_schema( xsd );
}

std::string HydratePackerMoverCreator::keyname() const {
	return HydratePackerMover::mover_name();
}

protocols::moves::MoverOP
HydratePackerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new HydratePackerMover );
}

void HydratePackerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HydratePackerMover::provide_xml_schema( xsd );
}



}  // namespace hydrate
}  // namespace protocols

