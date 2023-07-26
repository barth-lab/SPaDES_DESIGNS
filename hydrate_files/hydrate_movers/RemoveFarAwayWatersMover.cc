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
#include <protocols/hydrate/RemoveFarAwayWatersMover.hh>
#include <protocols/hydrate/RemoveFarAwayWatersMoverCreator.hh>

// Protocols
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/hydrate/hydrate_functions.hh>
#include <protocols/hydrate/hydrate_remove_functions.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/AtomTreeDiffJobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Core
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

// Basic
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
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

// Construct tracer.
static basic::Tracer TR( "protocols.hydrate.RemoveFarAwayWatersMover" );

namespace protocols {
namespace hydrate {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////

// Default constructor
RemoveFarAwayWatersMover::RemoveFarAwayWatersMover(): protocols::moves::Mover( "RemoveFarAwayWatersMover" )
{
	type("Hydrate");
}

// Copy constructor
RemoveFarAwayWatersMover::RemoveFarAwayWatersMover( RemoveFarAwayWatersMover const & hyd ): protocols::moves::Mover( hyd )
{
	//copy_data(*this, hyd);
}

// Assignment operator
RemoveFarAwayWatersMover &
RemoveFarAwayWatersMover::operator=( RemoveFarAwayWatersMover const & hyd )
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
RemoveFarAwayWatersMover::~RemoveFarAwayWatersMover() {}

// Mover methods
std::string
RemoveFarAwayWatersMover::get_name() const
{
	return type();
}

protocols::moves::MoverOP
RemoveFarAwayWatersMover::clone() const
{
	return RemoveFarAwayWatersMoverOP( new RemoveFarAwayWatersMover( *this ) );
}

protocols::moves::MoverOP
RemoveFarAwayWatersMover::fresh_instance() const
{
	return RemoveFarAwayWatersMoverOP( new RemoveFarAwayWatersMover() );
}


/// @details
/// @param    <pose>: the structure to be moved
/// @remarks
void
RemoveFarAwayWatersMover::apply(Pose & pose)
{
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace pack::task;

	remove_far_away_waters( pose );

}

void
RemoveFarAwayWatersMover::copy_data(
	RemoveFarAwayWatersMover hyd_to,
	RemoveFarAwayWatersMover hyd_from)
{
}

void RemoveFarAwayWatersMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap
) {
}

void RemoveFarAwayWatersMover::register_options() {
	using namespace basic::options::OptionKeys;
}

std::string RemoveFarAwayWatersMover::mover_name() {
	return "RemoveFarAwayWatersMover";
}


utility::tag::XMLSchemaComplexTypeGeneratorOP
RemoveFarAwayWatersMover::complex_type_generator_for_remove_far_away_water_mover( utility::tag::XMLSchemaDefinition & )
{

	using namespace utility::tag;
	AttributeList attributes;

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attributes )
		.add_optional_name_attribute();
	return ct_gen;
}


void RemoveFarAwayWatersMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_remove_far_away_water_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Removes waters from system." )
		.write_complex_type_to_schema( xsd );
}

std::string RemoveFarAwayWatersMoverCreator::keyname() const {
	return RemoveFarAwayWatersMover::mover_name();
}

protocols::moves::MoverOP
RemoveFarAwayWatersMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RemoveFarAwayWatersMover );
}

void RemoveFarAwayWatersMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RemoveFarAwayWatersMover::provide_xml_schema( xsd );
}



}  // namespace hydrate
}  // namespace protocols

