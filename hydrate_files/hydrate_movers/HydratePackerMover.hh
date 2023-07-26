// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/hydrate/HydratePackerMover.hh
/// @brief The HydratePacker mover protocol
/// @detailed
/// @author Joaquin Ambia, Jason K. Lai

#ifndef INCLUDED_protocols_hydrate_HydratePackerMover_HH
#define INCLUDED_protocols_hydrate_HydratePackerMover_HH

// Protocol headers
#include <protocols/hydrate/HydratePackerMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Core and Basic headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <iostream>
#include <string>

namespace protocols {
namespace hydrate {

/// @brief
class HydratePackerMover: public protocols::moves::Mover {
public:
	typedef core::pack::task::TaskFactoryOP TaskFactoryOP;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	// Default constructor
	HydratePackerMover();

	// Constructor that allows to specity a score function and the protein flexibility (region to pack & minimize )
	HydratePackerMover(
		core::scoring::ScoreFunctionOP scorefxn,
		std::string protein_flexibility = "all"
	);

	// Copy constructor
	HydratePackerMover(HydratePackerMover const & hyd);

	// Assignment operator
	HydratePackerMover & operator=(HydratePackerMover const & hyd);

	// Destructor
	virtual ~HydratePackerMover();

	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;
	
	virtual void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &
	) override;

	/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_score_function(
		TagCOP,
		basic::datacache::DataMap const &
	);
	
	/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
	virtual void parse_task_operations(
		TagCOP,
		basic::datacache::DataMap const &
	);

	static void register_options();

	// setters
	void score_function( ScoreFunctionOP sf );
	virtual void task_factory( TaskFactoryOP tf );

	/// @brief  Apply the corresponding move to <input_pose>.
	virtual void apply(core::pose::Pose & input_pose);
	
	static std::string mover_name();

	static utility::tag::XMLSchemaComplexTypeGeneratorOP complex_type_generator_for_hydrate_packer_mover( utility::tag::XMLSchemaDefinition & xsd );

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	



private:
	// Initialize data members.
	void init();

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(HydratePackerMover hyd_to, HydratePackerMover hyd_from);

	core::scoring::ScoreFunctionOP scorefxn_; // Score function used in the protocol

	bool water_hybrid_sf_;
	std::string resfile_;
	std::string hyfile_;
	bool hydrate_all_;
	utility::vector1< core::Size > nofasol_;
	std::string weights_;
	bool pretalaris2013_;
	core::Real hbond_threshold_;
	core::Real water_rotamers_cap_;
	std::string extra_res_fa_;
	std::string extra_res_cen_;
	std::string spanfile_;
	bool membed_init_;
	bool mhbond_depth_;
	core::Real cst_weight_;
	core::Real cst_fa_weight_;

	// rosettascripts options
	// water_hybrid_sf_				score::water_hybrid_sf true
	// resfile_               in::file::resfile ""
	// hyfile_								hydrate::hyfile ""
	// hydrate_all_           hydrate::hydrate_all false
	// nofasol_               hydrate::ignore_fa_sol_at_positions ""
	// weights_               score::weights ""
	// pretalaris2013_        mistakes::restore_pre_talaris_2013_behavior false
	// hbond_threshold_       hydrate::hbond_threshold -0.1
	// water_rotamers_cap_    hydrate::water_rotamers_cap 500
	// extra_res_fa_		      in::file::extra_res_fa ""
	// extra_res_cen_		      in::file::extra_res_cen ""
	// spanfile_					    in::file::spanfile ""
	// membed_init_				    membrane::Membed_init false
	// mhbond_depth_					membrane::Mhbond_depth false
	// cst_weight_						constraints::cst_weight 1.0
	// cst_fa_weight_					constraints::cst_fa_weight 1.0

	std::string protein_flexibility_;
	core::Real partial_hydrate_dew_; // Fraction of water molecules that will be packed during the first step using rotamers
	// with two optimized hb (double edge water, dew); defaults to 0.75
	bool short_residence_time_mode_; // Triggers a different algorithm to add water
	core::Real near_water_threshold_;  // Distance between a wat molecule and a residue to consider it flexible; yumeng
	bool minimize_bb_where_packing_;  // The minimizer will also minimize the backbone in the regions that were packed; yumeng
};  // class HydratePackerMover

}  // namespace hydrate
}  // namespace protocols

#endif  // INCLUDED_protocols_hydrate_HydratePackerMover_HH
